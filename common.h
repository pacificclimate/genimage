#include <netcdfcpp.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <list>
#include <vector>
#include <sstream>

#include <boost/tokenizer.hpp>
#include <boost/shared_ptr.hpp>

using namespace std;
using namespace boost;

#define DAYS_400YRS 146097
#define DAYS_100YRS 36524
#define DAYS_4YRS 1461
#define DAYS_1LPYR 366
#define DAYS_1EQYR 360
#define DAYS_1YR 365

#define MAX_TOY 17

enum TIMESOFYEAR{JAN,FEB,MAR,APR,MAY,JUN,JUL,AUG,SEP,OCT,NOV,DEC,DJF,MAM,JJA,SON,ANN};

enum FILEBITS{DOT,EXPT,VAR,MODEL,RUN,FILENAME};
const int noleap_days[] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 };
const int leap_days[] = { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 };
const int equal_days[] = { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360 };

const int noleap_dpm[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
const int leap_dpm[] = { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
const int equal_dpm[] = { 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30 };

int do_binary_search(int number, int max, const int array[]) {
  int min = 0;
  max--;
  int offset = min + (max - min) / 2;
  while(min <= max) {
    if(array[offset] > number) {
      max = offset - 1;
    } else if(array[offset] < number) {
      min = offset + 1;
    } else {
      return offset;
    }
    offset = min + (max - min) / 2;
  }

  return -1;
}

int find_slot_in_range(int number, int max, const int array[]) {
  for(int i = 0; i < max - 1; i++) {
    if(array[i] <= number && array[i + 1] > number) {
      return i;
    }
  }
  return -1;
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
    month = find_slot_in_range(days, 13, leap_days);
  } else {
    // No Leap
    month = find_slot_in_range(days, 13, noleap_days);
  }
  if(month < 0) {
    printf("Error: days value not within range!\n");
  }
  
  return 12 * year + month;
}

int get_365day_total_months(int days) {
  int remainder = days % DAYS_1YR;
  int year = days / DAYS_1YR;
  int month = find_slot_in_range(remainder, 13, noleap_days);
  if(month < 0) {
    printf("Error: days value not within range!\n");
  }
  return 12 * year + month;
}

int get_360day_total_months(int days) {
  int remainder = days % 360;
  int year = days / 360;
  int month = find_slot_in_range(remainder, 13, equal_days);
  if(month < 0) {
    printf("Error: days value not within range!\n");
  }
  return 12 * year + month;
}

int get_total_months(string calendar_type, int days) {
  if(calendar_type == "gregorian") {
    return get_gregorian_total_months(days);
  } else if(calendar_type == "365_day") {
    return get_365day_total_months(days);
  } else if(calendar_type == "360_day") {
    return get_360day_total_months(days);
  } else {
    assert(false);
  }
  return -1;
}

// Convert a Gregorian year-month-day date into days since 0000-01-01
int date2days_greg(int start_year, int start_month, int start_day) {
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
  
  return days;
}

// Convert a 365-day year-month-day date into days since 0000-01-01
int date2days_365(int start_year, int start_month, int start_day) {
  return start_year * DAYS_1YR + noleap_days[start_month - 1] + (start_day - 1);
}

// Convert a 360-day year-month-day date into days since 0000-01-01
int date2days_360(int start_year, int start_month, int start_day) {
  return start_year * DAYS_1EQYR + equal_days[start_month - 1] + (start_day - 1);
}

// Convert a year-month-day date into days since 0000-01-01, given the supplied calendar type
int date2days(string calendar_type, int start_year, int start_month, int start_day) {
  if(calendar_type == "gregorian") {
    return date2days_greg(start_year, start_month, start_day);
  } else if(calendar_type == "365_day") {
    return date2days_365(start_year, start_month, start_day);
  } else if(calendar_type == "360_day") {
    return date2days_360(start_year, start_month, start_day);
  } else {
    assert(false);
  }
  return -1;
}

class FileRecord {
public:

  FileRecord(string filename, bool do_time_calcs = true): filename(filename), f(new NcFile(filename.c_str(), NcFile::ReadOnly)) {
    is_ok = f->is_valid();

    if(is_ok) {
      boost::char_separator<char> sep("/");
      tokenizer<char_separator<char> > tok(filename, sep);
      vector<string> tokens(tok.begin(), tok.end());

      var = tokens[VAR];
      expt = tokens[EXPT];
      model = tokens[MODEL];
      run = tokens[RUN];
      
      // Check if this is corrected data
      string fn = tokens[FILENAME];
      corrected = (fn.find("corrected") != string::npos);
      
      set_time_params(do_time_calcs);
    }
  }

  void set_time_params(bool do_time_calcs = true) {
    if(!is_ok)
      return;

    NcVar* t;
    if(f->get_dim("time")) {
      timeless = false;
      
      if((t = f->get_var("time"))) {
	int start_year, start_month, start_day;
	NcAtt* cal = t->get_att("calendar");
	NcAtt* units = t->get_att("units");

	if(cal && units) {
	  char* ctype = cal->as_string(0);
	  char* calendar_start = units->as_string(0);
	  calendar_type = ctype;
	  if(calendar_type == "standard")
	    calendar_type = "gregorian";
	  else if(calendar_type == "noleap")
	    calendar_type = "365_day";
	  delete[] ctype;

	  if(do_time_calcs) {
	    // Load in the base month
	    if(sscanf(calendar_start, "days since %i-%i-%i", &start_year, &start_month, &start_day) != 3) {
	      printf("Failure to match start date\n");
	      return;
	    }
	    this->start_day = date2days(calendar_type, start_year, start_month, start_day);
	  }
	  delete[] calendar_start;
	}
	delete cal;
	delete units;
      }
    } else {
      timeless = true;
    }
  }
  
  bool operator==(const FileRecord& f) {
    return (f.var == var) && (f.expt == expt) && (f.model == model) && (f.run == run);
  }

  bool operator!=(const FileRecord& f) {
    return !(*this == f);
  }

  string filename;

  string var;
  string expt;
  string model;
  string run;

  shared_ptr<NcFile> f;
  bool corrected;
  bool timeless;

  string calendar_type;
  int start_day;

  bool is_ok;
};

class DataRecord {
public:
  DataRecord(FileRecord& f, NcVar* v, int offset, int month): f(f) {
    this->v = v;
    this->offset = offset;
    this->month = month;
  }
  NcVar *v;
  int offset;
  int month;
  FileRecord& f;

  // Should check more but this is adequate
  bool operator==(const DataRecord& r) {
    return (month == r.month && offset == r.offset);
  }

  bool operator<(const DataRecord& r) {
    return (month < r.month);
  }
};

template <typename T> void copy_att(NcAtt* src, T* dst) {
  NcValues* values = src->values();
  int num_vals = src->num_vals();
  
  // So braindead
  // Copy the attribute's contents to a temporary and add the attribute to the new output variable
  switch(src->type()) {
  case ncByte:
  case ncChar:
    {
      char data[num_vals];
      for(int j = 0; j < num_vals; j++) {
	data[j] = values->as_char(j);
      }
      dst->add_att(src->name(), num_vals, data);
    }
    break;
  case ncShort:
    {
      short data[num_vals];
      for(int j = 0; j < num_vals; j++) {
	data[j] = values->as_short(j);
      }
      dst->add_att(src->name(), num_vals, data);
    }
    break;
  case ncLong:
    {
      long data[num_vals];
      for(int j = 0; j < num_vals; j++) {
	data[j] = values->as_int(j);
      }
      dst->add_att(src->name(), num_vals, data);
    }
    break;
  case ncFloat:
    {
      float data[num_vals];
      for(int j = 0; j < num_vals; j++) {
	data[j] = values->as_float(j);
      }
      dst->add_att(src->name(), num_vals, data);
    }
    break;
  case ncDouble:
    {
      double data[num_vals];
      for(int j = 0; j < num_vals; j++) {
	data[j] = values->as_double(j);
      }
      dst->add_att(src->name(), num_vals, data);
    }
    break;
  case ncNoType:
    assert(false);
    break;
  }

  delete values;
}

template <typename T> void copy_atts(T* src, T* dst) {
  // Copy attributes
  int num_atts = src->num_atts();
  for(int i = 0; i < num_atts; i++) {
    NcAtt* att = src->get_att(i);
    
    copy_att(att, dst);

    delete att;
  }
}

void populate_dimvec(NcVar* src, NcFile& dst, vector<NcDim*>& dims) {
  int num_dims = src->num_dims();
  for(int i = 0; i < num_dims; i++) {
    const string dimname = src->get_dim(i)->name();
    NcDim* d = dst.get_dim(dimname.c_str());
    assert(d);
    dims.push_back(d);
  }
}

// Fancy way:
// Recurse until 2 or less values; then copy that thing, adding it to the next layer... or somethin

// Simple way: Just copy the values
template <typename T> void copy_var_values(NcVar* src, NcVar* dst, T* temporary) {
  long* edges = src->edges();

  assert(src->get(temporary, edges));
  assert(dst->put(temporary, edges));

  delete[] edges;
}

NcVar* copy_var_structure(NcVar* src, NcFile& dst) {
  // Get dimensions to use
  vector<NcDim*> dims;
  populate_dimvec(src, dst, dims);
  
  // Create output variable
  return dst.add_var(src->name(), src->type(), dims.size(), (const NcDim**)&dims[0]);
}

template <typename T> NcVar* copy_var_t(NcVar* src, NcFile& dst, T* temporary = 0) {
  NcVar* dstvar = copy_var_structure(src, dst);
  copy_atts(src, dstvar);
  if(!temporary) {
    temporary = new T[src->num_vals()];
    copy_var_values<T>(src, dstvar, temporary);
    delete[] temporary;
  } else {
    copy_var_values<T>(src, dstvar, temporary);
  }

  return dstvar;
}

NcVar* copy_var(NcVar* src, NcFile& dst) {
  switch(src->type()) {
  case ncByte:
  case ncChar:
    return copy_var_t<char>(src, dst);
    break;
  case ncShort:
    return copy_var_t<short>(src, dst);
    break;
  case ncLong:
    return copy_var_t<long>(src, dst);
    break;
  case ncFloat:
    return copy_var_t<float>(src, dst);
    break;
  case ncDouble:
    return copy_var_t<double>(src, dst);
    break;
  case ncNoType:
    return 0;
  }
  return 0;
}

template<typename T> NcDim* copy_dim(NcFile& src, NcFile& dst, T dim) {
  NcDim* d = src.get_dim(dim);
  NcDim* out;
  if(d->is_unlimited()) {
    out = dst.add_dim(d->name());
  } else {
    out = dst.add_dim(d->name(), d->size());
  }
  return out;
}

void copy_dims(NcFile& src, NcFile& dst) {
  // Copy the dimensions
  for(int i = 0; i < src.num_dims(); i++) {
    assert(copy_dim(src, dst, i));
  }
}

void copy_double_to_float(NcVar* src, NcVar* dst, int bufsize, long* edges, double* indata, float* outdata) {
    assert(src->get(indata, edges));
    for(int i = 0; i < bufsize; i++) {
      outdata[i] = (float)indata[i];
    }
    assert(dst->put(outdata, edges));
}

void copy_double_to_double(NcVar* src, NcVar* dst, int bufsize, long* edges, double* buf) {
    assert(src->get(buf, edges));
    assert(dst->put(buf, edges));
}

void copy_float_to_float(NcVar* src, NcVar* dst, int bufsize, long* edges, float* buf) {
    assert(src->get(buf, edges));
    assert(dst->put(buf, edges));
}

template<class T> class Range {
public:
  Range(T min, T max): min(min), max(max) {
  }
  T min, max;
};

void add_to_grid(int size, float* input, float* accum) {
  for(int i = 0; i < size; i++) {
    accum[i] += input[i];
  }
}

void add_to_grid(int size, float* input, float* accum, float missing) {
  for(int i = 0; i < size; i++) {
    if(accum[i] != missing) {
      if(input[i] == missing) {
	accum[i] = missing;
      } else {
	accum[i] += input[i];
      }
    }
  }
}

void add_to_grid(int size, float input, float* accum) {
  for(int i = 0; i < size; i++) {
    accum[i] += input;
  }
}

void add_to_grid(int size, float input, float* accum, float missing) {
  for(int i = 0; i < size; i++) {
    if(accum[i] != missing) {
      if(input == missing) {
	accum[i] = missing;
      } else {
	accum[i] += input;
      }
    }
  }
}

template <typename T> void divide_grid_by_scalar(int size, float* input, T scalar) {
  for(int i = 0; i < size; i++) {
    input[i] /= scalar;
  }
}

template <typename T> void divide_grid_by_scalar(int size, float* input, T scalar, float missing) {
  for(int i = 0; i < size; i++) {
    if(input[i] != missing) {
      if(input[i] == missing) {
	input[i] = missing;
      } else {
	input[i] /= scalar;
      }
    }
  }
}

template <typename T> void multiply_grid_by_scalar(int size, float* input, T scalar) {
  for(int i = 0; i < size; i++) {
    input[i] *= scalar;
  }
}

template <typename T> void multiply_grid_by_scalar(int size, float* input, T scalar, float missing) {
  for(int i = 0; i < size; i++) {
    if(input[i] != missing) {
      if(input[i] == missing) {
	input[i] = missing;
      } else {
	input[i] *= scalar;
      }
    }
  }
}

int get_recsize_and_edges(NcVar* invar, long* edges) {
  // Construct the list of what to copy
  int recsize = 1;
  for(int i = 0; i < invar->num_dims(); i++) {
    NcDim* d = invar->get_dim(i);
    assert(d);
    string name = d->name();
    if(name == "time") {
      edges[i] = 1;
    }
    recsize *= edges[i];
  }
  return recsize;
}

