//#include <iostream>
#include "common.h"

#define BUFSIZE 10240

void copy_from_list(vector<float>& l, NcFile* f, char* var) {
  NcVar* v = f->get_var(var);
  long* edges = v->edges();
  int recsize = get_recsize_and_edges(v, edges);
  float* data = new float[recsize];

  assert((unsigned int)recsize == l.size());

  v->get(data, edges);
  copy(l.begin(), l.end(), data);
  v->put(data, edges);

  delete[] data;
  delete[] edges;
}

void load_list_from_csv(vector<float>& l, FILE* f) {
  char buf[BUFSIZE];
  boost::char_separator<char> sep(",");
  
  while(fgets(buf, BUFSIZE, f)) {
    string line = buf;
    tokenizer<char_separator<char> > tok(line, sep);
    tokenizer<char_separator<char> >::const_iterator i;
    for(i = tok.begin(); i != tok.end(); ++i)
      l.push_back(atof((*i).c_str()));
  }
}

void reverse_width(vector<float>& l, int height, int width) {
  for(int i = 0; i < height; i++) {
    for(int j = 0; j < floor(width / 2); j++) {
      swap(l[(i * height) + j], l[(i * height) + (width - j - 1)]);
    }
  }
}

void rescan_data(vector<float>& l, int height, int width) {
  float* data = new float[height * width];
  int k = 0;
  for(int j = 0; j < width; j++) {
    for(int i = 0; i < height; i++) {
      data[((height - i - 1) * width) + j] = l[k];
      k++;
    }
  }

  copy(data, data + (height * width), l.begin());
  delete[] data;
}

void fix_longs(vector<float>& l) {
  for(vector<float>::iterator i = l.begin(); i != l.end(); i++)
    *i = fmod((*i + 180), 360) - 180;
}

#define WIDTH 80
#define HEIGHT 90

int main(int argc, char** argv) {
  char buf[BUFSIZE];

  // Try not to fall on your face, netcdf, when a dimension or variable is missing
  ncopts = NC_VERBOSE;

  if(argc < 3) {
    printf("Usage: splice_lats_longs <lat file> <long file>\n");
    exit(1);
  }

  FILE* latfile = fopen(argv[1], "r");
  FILE* lonfile = fopen(argv[2], "r");
  vector<float> newlats;
  vector<float> newlons;
  
  assert(latfile && lonfile);

  load_list_from_csv(newlats, latfile);
  load_list_from_csv(newlons, lonfile);
  fix_longs(newlons);

  rescan_data(newlats, HEIGHT, WIDTH);
  rescan_data(newlons, HEIGHT, WIDTH);
  //reverse_width(newlats, HEIGHT, WIDTH);
  //reverse_width(newlons, HEIGHT, WIDTH);

  fclose(latfile);
  fclose(lonfile);

  // Second step: Read the files in, one by one; extract the climatologies
  while(fgets(buf, BUFSIZE, stdin)) {
    *(strchr(buf, '\n')) = '\0';
    string filename = buf;
    NcFile* f = new NcFile(buf, NcFile::Write);

    if(!f->is_valid()) {
      printf("Failed to open file %s\n", buf);
      exit(1);
    }

    copy_from_list(newlats, f, "lat");
    copy_from_list(newlons, f, "lon");

    delete f;
  }
}
