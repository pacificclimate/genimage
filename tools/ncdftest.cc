//#include <iostream>
#include <netcdfcpp.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <list>

#include <boost/tokenizer.hpp>

using namespace std;
using namespace boost;

class FileRecord {
public:
  string filename;
  string calendar_start;
  string calendar_type;
  double start_day;
  double end_day;
  bool timeless;
  int start_month;
  int end_month;
  string var;
  string expt;
  string model;
  string run;
};

#define DAYS_400YRS 146097
#define DAYS_100YRS 36524
#define DAYS_4YRS 1461
#define DAYS_1LPYR 366
#define DAYS_1YR 365

enum FILEBITS{DOT,EXPT,VAR,MODEL,RUN,FILENAME};
const int noleap_days[] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 };
const int leap_days[] = { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 };
const int equal_days[] = { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360 };

int do_binary_search(int number, int max, const int array[]) {
  int min = 0;
  int offset = (min + max) / 2;
  max--;
  if(number < array[0] && number > array[max]) {
    return -1;
  }
  while(array[offset] > number || array[offset + 1] < number) {
    if(array[offset] < number) {
      min = (int)ceil((double)(max + min) / 2.0);
    } else {
      max = (int)floor((double)(max + min) / 2.0);
    }
    offset = (max + min) / 2;
    //printf("moo: [%i < %i < %i]\n", array[offset], number, array[offset + 1]);
  }
  return offset;
}

int get_gregorian_total_months(int days) {
  int year, month;
  
  int yrs_400 = days / DAYS_400YRS;
  days = days % DAYS_400YRS;
  int yrs_100 = days / DAYS_100YRS;
  days = days % DAYS_100YRS;
  int yrs_4 = days / DAYS_4YRS;
  days = days % DAYS_4YRS;
  int yrs = days / DAYS_1YR;
  days = days % DAYS_1YR;
  
  year = yrs_400 * 400 + yrs_100 * 100 + yrs_4 * 4 + yrs;
  
  if(yrs == 0 && (yrs_4 != 0 || (yrs_4 == 0 && yrs_100 == 0))) { 
    // Leap
    month = do_binary_search(days, 12, leap_days);
  } else {
    // No Leap
    month = do_binary_search(days, 12, noleap_days);
  }
  if(month < 0) {
    printf("Error: days value not within range!\n");
  }
  
  return year * 12 + month;
}

int get_365day_total_months(int start_month, int days) {
  int remainder = days % DAYS_1YR;
  int year = days / DAYS_1YR;
  int month = do_binary_search(remainder, 12, noleap_days);
  if(month < 0) {
    printf("Error: days value not within range!\n");
  }
  return start_month + 12 * year + month;
}

int get_360day_total_months(int start_month, int days) {
  int remainder = days % 360;
  int year = days / 360;
  int month = do_binary_search(remainder, 12, equal_days);
  if(month < 0) {
    printf("Error: days value not within range!\n");
  }
  return start_month + month + 12 * year;
}

void calc_month_offsets(FileRecord* fr, int start_year, int start_month, int start_day) {
  int begin_month = start_year * 12 + start_month - 1;
  if(fr->calendar_type == "365_day" || fr->calendar_type == "noleap") {
    // Calculate the offset in months that the offset in days expresses, for a calendar without leap years
    fr->start_month = get_365day_total_months(begin_month, (int)fr->start_day);
    fr->end_month = get_365day_total_months(begin_month, (int)fr->end_day);
  } else if(fr->calendar_type == "gregorian" || fr->calendar_type == "standard") {
    // Calculate the offset in months that the offset in days expresses, for a Gregorian calendar with leap years
    int leap_years = (start_year / 4) - (start_year / 100) + (start_year / 400);
    int noleap_years = start_year - leap_years;
    int days = leap_years * DAYS_1LPYR + noleap_years * DAYS_1YR + (start_day - 1);
    
    if(start_year % 4 == 0 && (start_year % 100 != 0 || start_year % 400 == 0)) {
      // Start year is a leap year
      days += leap_days[start_month - 1];
    } else {
      // Start year is not a leap year
      days += noleap_days[start_month - 1];
    }
    fr->start_month = get_gregorian_total_months(days + (int)fr->start_day);
    fr->end_month = get_gregorian_total_months(days + (int)fr->end_day);
  } else if(fr->calendar_type == "360_day") {
    // Calculate the offset in months that the offset in days expresses, for a calendar with 30-day months
    fr->start_month = get_360day_total_months(begin_month, (int)fr->start_day);
    fr->end_month = get_360day_total_months(begin_month, (int)fr->end_day);
  } else {
    printf("Error: Invalid calendar type\n");
  }
}

int main(int argc, char** argv) {
  NcFile* f;
  char buf[1024];
  list<FileRecord*> l;
  list<FileRecord*>::const_iterator li;

  while(fgets(buf, 1024, stdin)) {
    FileRecord* fr = new FileRecord();
    *(strchr(buf, '\n')) = '\0';
    f = new NcFile(buf, NcFile::ReadOnly);

    // Try not to fall on your face, netcdf, when a dimension or variable is missing
    ncopts = NC_VERBOSE;

    if(!f->is_valid()) {
      printf("Failed to open file %s\n", buf);
      delete f;
      continue;
    }

    fr->filename = buf;
    boost::char_separator<char> sep("/");
    tokenizer<char_separator<char> > tok(fr->filename, sep);
    int i = 0;
    for(tokenizer<char_separator<char> >::iterator beg = tok.begin(); beg != tok.end(); ++beg, ++i) {
      switch(i) {
      case VAR:
	fr->var = *beg;
	break;
      case EXPT:
	fr->expt = *beg;
	break;
      case MODEL:
	fr->model = *beg;
	break;
      case RUN:
	fr->run = *beg;
	break;
      }
    }

    if(f->get_dim("time")) {
      NcVar *t;
      fr->timeless = false;
      
      if((t = f->get_var("time"))) {
	int start_year, start_month, start_day;
	long* edges;
	int len = t->get_dim(0)->size();
	double* values = new double[len];
	t->get(values, len);

	fr->start_day = values[0];
	fr->end_day = values[len - 1];
	NcAtt* cal = t->get_att("calendar");
	NcAtt* units = t->get_att("units");

	fr->calendar_type = cal->as_string(0);
	fr->calendar_start = units->as_string(0);
	
	// Load in the base month
	if(sscanf(fr->calendar_start.c_str(), "days since %i-%i-%i", &start_year, &start_month, &start_day) != 3) {
	  printf("Failure to match start date\n");
	}

	// Compute the offsets in months
	calc_month_offsets(fr, start_year, start_month, start_day);

	delete[] values;
      }
    } else {
      fr->timeless = true;
    }

    delete f;
    l.push_back(fr);
  }


  for(li = l.begin(); li != l.end(); li++) {
    FileRecord* fr = *li;

    printf("%s,%s,%s,%s,%s,", fr->filename.c_str(), fr->model.c_str(), fr->expt.c_str(), fr->var.c_str(), fr->run.c_str());

    if(fr->timeless) {
      cout << "0,,,,,\n";
    } else {
      printf("1,%s,%s,%i,%f,%i,%f,\n", fr->calendar_type.c_str(), fr->calendar_start.c_str(), fr->start_month, fr->start_day, fr->end_month, fr->end_day);
    }

    delete fr;
  }
}
