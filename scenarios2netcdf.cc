//#include <iostream>
#include <netcdfcpp.h>
#include "genimage.h"

#define TIMESOFYEAR 17

using namespace std;

class Datafile {
public:
  Datafile(char* filename) {
    this->filename = filename;
    data = 0;
    timeofyear = 0;
    slmask = 0;
  }

  Datafile() {
    filename = 0;
    data = 0;
    timeofyear = 0;
    slmask = 0;
  }

  ~Datafile() {
    if(data) {
      delete[] data;
    }
    if(slmask) {
      delete[] slmask;
    }
  }

  bool openFile(bool hasHeader = true) {
    char buf[1024];
    char buf2[1024];
    char * ptr = buf;
    int i;

    // Open data file
    infile = fopen(filename, "r");
  
    if(!infile) {
      return false;
    }
    if(hasHeader) {
      if(fgets(buf, 1024, infile)) {
	headerlen = strlen(buf);
	ptr = buf;
	
	// Assume first line is parameters -- valid.
	for(i = 0; i < 4; i++) {
	  ptr = seektoword(ptr);
	}
	
	// At this point, we're at the # of cols
	cols = readint(ptr, buf2);
	
	// Now seek to the next non-space
	ptr = seektoword(ptr);
	
	// Now load in # of rows
	rows = readint(ptr, buf);
	
	line_length = cols * VALUE_LENGTH + 3;
	
      } else {
	return false;
      }
    }
    if(!data) {
      data = new double[rows * cols * TIMESOFYEAR];
    }
    if(!slmask) {
      slmask = new int[rows * cols * TIMESOFYEAR];
    }
    return true;
  }

  bool closeFile() {
    fclose(infile);
  }

  bool loadSlmask() {
    char * dataline;

    dataline = new char[line_length];
    load_grid(infile, rows, cols, dataline, line_length, slmask, SLMASK_LENGTH);
    delete[] dataline;
    return true;
  }
  
  bool loadLats() {
    char * dataline;

    dataline = new char[line_length];
    load_grid(infile, rows, 1, dataline, line_length, data, LL_LENGTH);
    delete[] dataline;
    return true;
  }

  bool loadLongs() {
    char * dataline;

    dataline = new char[line_length];
    load_grid(infile, 1, cols, dataline, line_length, data, LL_LENGTH);
    delete[] dataline;
    return true;
  }

  bool loadData() {
    char * dataline;
    int i, timeofyear;
    double* dataptr = data;
    char buf[1024];

    dataline = new char[line_length];
    fseek(infile, 0, SEEK_SET);
    while(fgets(buf, 1024, infile)) {
      load_grid(infile, rows, cols, dataline, line_length, dataptr, VALUE_LENGTH);
      dataptr += (rows * cols);
    }
    delete[] dataline;
    return true;
  }

  double* data;
  int* slmask;
  char* filename;
  int cols, rows;
  int timeofyear;

private:
  int headerlen;
  int line_length;
  FILE* infile;
};

#define NUM_TS 4
#define NUM_VARS 20


int main(int argc, char** argv) {
  int i, j;
  Datafile d;
  NcFile* f;
  int timeofyear = 0;
  char buf[1024];
  char buf2[1024];
  char *c;
  NcDim* rows;
  NcDim* cols;
  NcDim* timesofyear;
  int firstRun = 1;

  f = new NcFile(argv[1], NcFile::Replace);

  while(fgets(buf, 1024, stdin)) {
    c = strchr(buf, '\n');
    *c = '\0';
    d.filename = buf;
    strcpy(buf2, buf);
    c = strrchr(buf2, '_');
    *c = '\0';
    c = strchr(buf2, '_') + 1;
    if(strstr(d.filename, "slmask")) {
      //printf("Load slmask\n");
      if(!d.openFile(false)) {
	printf("Missing file %s\n", d.filename);
      } else {
	NcVar* var1 = f->add_var(c, ncInt, rows, cols);
	d.loadSlmask();
	var1->put(d.slmask, d.rows, d.cols);
	d.closeFile();
      }
    } else if(strstr(d.filename, "lats")) {
      //printf("Load lats\n");
      if(!d.openFile(false)) {
	printf("Missing file %s\n", d.filename);
      } else {
	NcVar* var1 = f->add_var(c, ncDouble, rows);
	d.loadLats();
	var1->put(d.data, d.rows);
	d.closeFile();
      }
    } else if(strstr(d.filename, "longs")) {
      //printf("Load longs\n");
      if(!d.openFile(false)) {
	printf("Missing file %s\n", d.filename);
      } else {
	NcVar* var1 = f->add_var(c, ncDouble, cols);
	d.loadLongs();
	var1->put(d.data, d.cols);
	d.closeFile();
      }
    } else {
      //printf("Load datafile\n");
      if(!d.openFile()) {
	printf("Missing file %s\n", d.filename);
      } else {
	if(firstRun) {
	  timesofyear = f->add_dim("timesofyear", TIMESOFYEAR);
	  rows = f->add_dim("rows", d.rows);
	  cols = f->add_dim("columns", d.cols);
	  firstRun = 0;
	}
	NcVar* var1 = f->add_var(c, ncDouble, timesofyear, rows, cols);
	d.loadData();
	var1->put(d.data, TIMESOFYEAR, d.rows, d.cols);
	d.closeFile();
      }
    }
    f->sync();
  }


  delete f;
}
