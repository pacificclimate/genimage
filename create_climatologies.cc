//#include <iostream>
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
  int old_offset = 0;
  int offset = (min + max) / 2;
  max--;
  if(number < array[0] || number > array[max]) {
    return -1;
  }
  while(old_offset != offset && (array[offset] > number || array[offset + 1] < number)) {
    if(array[offset] < number) {
      min = (int)ceil((double)(max + min) / 2.0);
    } else {
      max = (int)floor((double)(max + min) / 2.0);
    }
    old_offset = offset;
    offset = (max + min) / 2;
  }
  if(array[offset] > number || array[offset + 1] < number) {
    return -1;
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
    month = do_binary_search(days, 13, leap_days);
  } else {
    // No Leap
    month = do_binary_search(days, 13, noleap_days);
  }
  if(month < 0) {
    printf("Error: days value not within range!\n");
  }
  
  return 12 * year + month;
}

int get_365day_total_months(int days) {
  int remainder = days % DAYS_1YR;
  int year = days / DAYS_1YR;
  int month = do_binary_search(remainder, 13, noleap_days);
  if(month < 0) {
    printf("Error: days value not within range!\n");
  }
  return 12 * year + month;
}

int get_360day_total_months(int days) {
  int remainder = days % 360;
  int year = days / 360;
  int month = do_binary_search(remainder, 13, equal_days);
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

  FileRecord(string filename): filename(filename), f(new NcFile(filename.c_str(), NcFile::ReadOnly)) {
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
      
      set_time_params();
    }
  }

  void set_time_params() {
    if(!is_ok)
      return;

    NcVar* t;
    if(f->get_dim("time")) {
      timeless = false;
      
      if((t = f->get_var("time"))) {
	int start_year, start_month, start_day;
	NcAtt* cal = t->get_att("calendar");
	NcAtt* units = t->get_att("units");
	char* ctype = cal->as_string(0);
	char* calendar_start = units->as_string(0);

	calendar_type = ctype;
	if(calendar_type == "standard")
	  calendar_type = "gregorian";
	else if(calendar_type == "noleap")
	  calendar_type = "365_day";

	// Load in the base month
	if(sscanf(calendar_start, "days since %i-%i-%i", &start_year, &start_month, &start_day) != 3) {
	  printf("Failure to match start date\n");
	  return;
	}

	delete[] ctype;
	delete[] calendar_start;
	delete cal;
	delete units;
	this->start_day = date2days(calendar_type, start_year, start_month, start_day);
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

void copy_att(NcAtt* src, NcVar* dst) {
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
  }

  delete values;
}

void copy_atts(NcVar* src, NcVar* dst) {
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
  default:
    return 0;
  }
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

template <typename T> void divide_grid_by_scalar(int size, float* input, T scalar) {
  for(int i = 0; i < size; i++) {
    input[i] /= scalar;
  }
}

template <typename T> void multiply_grid_by_scalar(int size, float* input, T scalar) {
  for(int i = 0; i < size; i++) {
    input[i] *= scalar;
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

void create_climatology(FileRecord& f, string outpath, const Range<int>& r) {
  string path = outpath;
  struct stat s;
  list<string>::const_iterator p;
  list<string> pathbits;
  pathbits.push_back(f.expt);
  pathbits.push_back(f.var);
  pathbits.push_back(f.model);
  pathbits.push_back(f.run);

  for(p = pathbits.begin(); p != pathbits.end(); ++p) {
    path += "/" + *p;

    // If directory in path doesn't exist, create it.
    if(stat(path.c_str(), &s) != 0) {
      mkdir(path.c_str(), 0777);
    }
  }

  stringstream sst;
  sst << path << "/" << f.model << "-" << f.expt << "-" << f.var << "-" << f.run << "-" << (r.min / 12) << "-" << (r.max / 12) << ".nc";

  string ofile = path + "/" + sst.str();

  NcFile out(ofile.c_str(), NcFile::Replace);
  NcFile& in = *(f.f);

  assert(in.is_valid() && out.is_valid());

  out.set_fill(NcFile::NoFill);

  copy_dims(in, out);

  // Copy vars
  int num_vars = in.num_vars();

  for(int i = 0; i < num_vars; i++) {
    NcVar* v = in.get_var(i);
    assert(v);
    string name = v->name();
    if(name == "time" || name == "time_bnds" || name == f.var)
      continue;
    copy_var(v, out);
  }
  
  NcVar* invar = in.get_var(f.var.c_str());
  long* edges = invar->edges();
  int rec_size = get_recsize_and_edges(invar, edges);
  if(f.timeless) {
    // Just copy the damned thing
  } else {
    // Do the averaging
    
    // Get time
    NcVar* intime = in.get_var("time");
    long* tedges = intime->edges();
    int numtimes = intime->num_vals();
    int times[numtimes];
    intime->get(times, tedges);
    delete tedges;

    // Allocate data
    int days[MAX_TOY];
    int data_size = rec_size * MAX_TOY;
    float* data = new float[data_size];

    // Clear data
    for(int i = 0; i < data_size; i++) {
      data[i] = 0;
    }
    for(int i = 0; i < MAX_TOY; i++) {
      days[i] = 0;
    }

    // Loop over the data, starting from the start location
    int start_offset = do_binary_search(r.min, numtimes, times);
    assert(start_offset != -1);
    float* indata = new float[invar->rec_size()];
    for(int i = start_offset; times[i] <= r.max; i++) {
      invar->set_cur(i);
      invar->get(indata, edges);
      int days_in_month;
      int year = times[i] / 12;
      int month = times[i] % 12;

      // Get the right # of days for this month
      if(f.calendar_type == "gregorian") {
	if(year % 4 == 0 && (year % 100 != 0 || year % 400 == 0)) {
	  // leap year
	  days_in_month = leap_dpm[month];
	} else {
	  // not leap year
	  days_in_month = noleap_dpm[month];
	}
      } else if(f.calendar_type == "365_day") {
	days_in_month = noleap_dpm[month];
      } else if(f.calendar_type == "360_day") {
	days_in_month = equal_dpm[month];
      } else {
	assert(false);
      }

      // Multiply by # of days to give accumulated days of the mean in the month
      multiply_grid_by_scalar(rec_size, indata, days_in_month);

      // Do month
      add_to_grid(rec_size, indata, data + (rec_size * month));
      days[month] += days_in_month;
    }

    // Divide grids by # days
    // WARNING: DOES NOT HANDLE MISSING VALUES
    // FIXME BUG BUG BUG
    // FIXME MAYBE USE DAYS OF MONTHS?
    for(int i = JAN; i <= DEC; i++) {
      divide_grid_by_scalar(rec_size, data + (rec_size * i), days[i]);
      add_to_grid(rec_size, data + (rec_size * i), data + (rec_size * ANN));
    }
    add_to_grid(rec_size, data + (rec_size * DEC), data + (rec_size * DJF));
    add_to_grid(rec_size, data + (rec_size * JAN), data + (rec_size * DJF));
    add_to_grid(rec_size, data + (rec_size * FEB), data + (rec_size * DJF));
    add_to_grid(rec_size, data + (rec_size * MAR), data + (rec_size * MAM));
    add_to_grid(rec_size, data + (rec_size * APR), data + (rec_size * MAM));
    add_to_grid(rec_size, data + (rec_size * MAY), data + (rec_size * MAM));
    add_to_grid(rec_size, data + (rec_size * JUN), data + (rec_size * JJA));
    add_to_grid(rec_size, data + (rec_size * JUL), data + (rec_size * JJA));
    add_to_grid(rec_size, data + (rec_size * AUG), data + (rec_size * JJA));
    add_to_grid(rec_size, data + (rec_size * SEP), data + (rec_size * SON));
    add_to_grid(rec_size, data + (rec_size * OCT), data + (rec_size * SON));
    add_to_grid(rec_size, data + (rec_size * NOV), data + (rec_size * SON));

    divide_grid_by_scalar(rec_size, data + (rec_size * DJF), 3);
    divide_grid_by_scalar(rec_size, data + (rec_size * MAM), 3);
    divide_grid_by_scalar(rec_size, data + (rec_size * JJA), 3);
    divide_grid_by_scalar(rec_size, data + (rec_size * SON), 3);
    divide_grid_by_scalar(rec_size, data + (rec_size * ANN), 12);

    // Finally throw the data into the output
    // FIXME ACTUALLY DO THIS

    delete[] edges;
    delete[] indata;
    delete[] data;
  }
}

int main(int argc, char** argv) {
  char buf[1024];
  list<Range<int> > ranges;

  // Try not to fall on your face, netcdf, when a dimension or variable is missing
  ncopts = NC_VERBOSE;

  if(argc < 2) {
    printf("Usage: fix_missing <output_path> [<climatology start> <climatology end>...]");
  }

  string output_path = argv[1];

  // First step: Accumulate the climatologies (provided on command line after filename)
  for(int i = 2; i + 1 < argc; i += 2) {
    Range<int> r(atoi(argv[i]), atoi(argv[i + 1]));
    ranges.push_back(r);
  }

  // Second step: Read the files in, one by one; extract the climatologies
  while(fgets(buf, 1024, stdin)) {
    // Chomp a la perl
    *(strchr(buf, '\n')) = '\0';
    FileRecord fr(buf);

    if(!fr.is_ok) {
      printf("Failed to open file %s\n", buf);
      continue;
    }

    // Generate climatologies
    list<Range<int> >::const_iterator i;
    for(i = ranges.begin(); i != ranges.end(); i++) {
      create_climatology(fr, output_path, *i);
    }
  }
}
