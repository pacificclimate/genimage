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
#include <string.h>
#include <stdlib.h>


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

#define DAYS_1970 719528
#define MAX_TOY 17

enum TIMESOFYEAR{JAN,FEB,MAR,APR,MAY,JUN,JUL,AUG,SEP,OCT,NOV,DEC,DJF,MAM,JJA,SON,ANN};

enum PATHBITS{DOT,EXPT,VAR,MODEL,RUN,FILENAME};
enum PATHBITS_CMIP5{C5DOT,C5CENTER,C5MODEL,C5EXPT,C5FREQ,C5TYPE,C5JUNK,C5RUN,C5VER,C5VAR,C5FILENAME};
enum FILEBITS{FMODEL,FEXPT,FVAR,FRUN,YEARSTART,YEAREND};
const int noleap_days[] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 };
const int leap_days[] = { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 };
const int equal_days[] = { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360 };

const int noleap_dpm[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
const int leap_dpm[] = { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
const int equal_dpm[] = { 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30 };

vector<int> get_time_values(NcFile& f, int start_day);
float get_missing_value_float(NcVar* v);
double get_missing_value(NcVar* v);
int do_binary_search(int number, int max, const int array[]);
int find_slot_in_range(int number, int max, const int array[]);
int get_gregorian_total_months(int days);
int get_365day_total_months(int days);
int get_360day_total_months(int days);
int get_total_months(string calendar_type, int days);
int date2days_greg(int start_year, int start_month, int start_day);
int date2days_365(int start_year, int start_month, int start_day);
int date2days_360(int start_year, int start_month, int start_day);
int date2days(string calendar_type, int start_year, int start_month, int start_day);
void populate_dimvec(NcVar* src, NcFile& dst, vector<NcDim*>& dims);
void populate_dimvec(NcVar* src, vector<NcDim*>& dims);
NcVar* copy_var_structure(NcVar* src, NcFile& dst);
NcVar* copy_var(NcVar* src, NcFile& dst);
void copy_dims(NcFile& src, NcFile& dst);
void copy_double_to_float(NcVar* src, NcVar* dst, int bufsize, long* edges, double* indata, float* outdata);
void copy_double_to_double(NcVar* src, NcVar* dst, int bufsize, long* edges, double* buf);
void copy_float_to_float(NcVar* src, NcVar* dst, int bufsize, long* edges, float* buf);
void add_to_grid(int size, float* input, float* accum);
void add_to_grid(int size, float* input, float* accum, float missing);
void add_to_grid(int size, float input, float* accum, float missing);
int get_recsize_and_edges(NcVar* invar, long* edges);

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

template <typename T> void divide_grid_by_scalar(int size, T* input, T scalar) {
  for(int i = 0; i < size; i++) {
    input[i] /= scalar;
  }
}

template <typename T> void divide_grid_by_scalar(int size, T* input, T scalar, float missing) {
  for(int i = 0; i < size; i++) {
    if(input[i] != missing) {
      input[i] /= scalar;
    } else {
      input[i] = missing;
    }
  }
}

template <typename T> void multiply_grid_by_scalar(int size, T* input, T scalar) {
  for(int i = 0; i < size; i++) {
    input[i] *= scalar;
  }
}

template <typename T> void multiply_grid_by_scalar(int size, T* input, T scalar, T missing) {
  for(int i = 0; i < size; i++) {
    if(input[i] != missing) {
      input[i] *= scalar;
    } else {
      input[i] = missing;
    }
  }
}

template <typename T> void add_scalar_to_grid(int size, T* input, T scalar) {
  for(int i = 0; i < size; i++) {
    input[i] += scalar;
  }
}

template <typename T> void add_scalar_to_grid(int size, T* input, T scalar, T missing) {
  for(int i = 0; i < size; i++) {
    if(input[i] != missing) {
      input[i] += scalar;
    } else {
      input[i] = missing;
    }
  }
}

template<class T> class Range {
public:
  Range(T min, T max): min(min), max(max) {
  }
  T min, max;
};

template <typename T> void copy_to_float(NcVar* src, NcVar* dst, int bufsize, long* edges, T* indata, float* outdata, float scale_factor, float bias) {
  assert(src->get(indata, edges));
  for(int i = 0; i < bufsize; i++) {
    outdata[i] = ((float)indata[i] * scale_factor) - bias;
  }
  assert(dst->put(outdata, edges));
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

// Fancy way:
// Recurse until 2 or less values; then copy that thing, adding it to the next layer... or somethin

// Simple way: Just copy the values
template <typename T> void copy_var_values(NcVar* src, NcVar* dst, T* temporary) {
  long* edges = src->edges();

  assert(src->get(temporary, edges));
  assert(dst->put(temporary, edges));

  delete[] edges;
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

template <typename T> vector<int> get_time_values(NcVar* time_var, int start_day) {
  // Get and calculate size of object
  long* edges = time_var->edges();
  int dim_prod = 1;

  for(int i = 0; i < time_var->num_dims(); ++i)
    dim_prod *= edges[i];

  // Allocate storage, get data.
  vector<T> vals(dim_prod);
  vector<int> retvals(dim_prod);
  assert(time_var->get(&vals[0], edges));

  // Then run through and calculate the day offsets, putting them in the vector.
  for(int i = 0; i < dim_prod; ++i)
    retvals[i] = (int)vals[i] + start_day;

  delete[] edges;
  return retvals;
}

class FileRecord {
public:

 FileRecord(string filename, bool do_time_calcs = true, bool do_filename_date_parsing = false, bool cmip5_paths = false): filename(filename), f(new NcFile(filename.c_str(), NcFile::ReadOnly)) {
    is_ok = f->is_valid();

    hourly = false;

    if(is_ok) {
      boost::char_separator<char> sep("/");
      tokenizer<char_separator<char> > tok(filename, sep);
      vector<string> tokens(tok.begin(), tok.end());

      if(cmip5_paths) {
        var = tokens[C5VAR];
        expt = tokens[C5EXPT];
        model = tokens[C5MODEL];
        run = tokens[C5RUN];
        file = tokens[C5FILENAME];
      } else {
        var = tokens[VAR];
        expt = tokens[EXPT];
        model = tokens[MODEL];
        run = tokens[RUN];
        file = tokens[FILENAME];
      }
      
      // Check if this is corrected data
      corrected = (file.find("corrected") != string::npos);
      
      set_time_params(do_time_calcs, do_filename_date_parsing);
    }
  }

  void set_time_params(bool do_time_calcs = true, bool do_filename_date_parsing = false) {
    if(!is_ok)
      return;

    NcVar* t;
    if(f->get_dim("time")) {
      timeless = false;
      
      int start_year, start_month, start_day;
      if(do_filename_date_parsing) {
	int monyear = 0;
	const char* stuff = strrchr(strrchr(filename.c_str(), '/'), '_');
	calendar_type = "365_day";
	printf("Filename: %s\n", filename.c_str());
	if(sscanf(stuff, "_%i.nc", &monyear) != 1) {
	  printf("Failure to match date\n");
	  return;
	}
	start_day = 1;
	start_month = monyear % 100;
	start_year = monyear / 100;
	hourly = true;
	this->start_day = date2days(calendar_type, start_year, start_month, start_day);
      } else {
	if(do_time_calcs && (t = f->get_var("time"))) {
	  NcAtt* cal = t->get_att("calendar");
	  NcAtt* units = t->get_att("units");

	  if(cal && units) {
	    char* ctype = cal->as_string(0);
	    char* calendar_start = units->as_string(0);
	    calendar_type = ctype;
	    if(calendar_type == "standard" || calendar_type == "proleptic_gregorian")
	      calendar_type = "gregorian";
	    else if(calendar_type == "noleap")
	      calendar_type = "365_day";
	    delete[] ctype;
	    
	    if(do_time_calcs) {
	      // Load in the base month
	      if(sscanf(calendar_start, "days since %i-%i-%i", &start_year, &start_month, &start_day) != 3) {
		if(sscanf(calendar_start, "months since %i-%i", &start_year, &start_month) == 2) {
		  start_day = 1;
		} else {
		  printf("Failure to match start date\n");
		  delete cal;
		  delete units;
		  return;
		}
	      }
	      printf("start y-m-d: %d-%d-%d\n", start_year, start_month, start_day);
	      this->start_day = date2days(calendar_type, start_year, start_month, start_day);
	    }
	    delete[] calendar_start;
	    delete cal;
	    delete units;
	  }
	}
      }
    } else {
      timeless = true;
    }
  }

  bool close() {
    f.reset();
    return true;
  }
  
  bool open() {
    f.reset(new NcFile(filename.c_str(), NcFile::ReadOnly));
    is_ok = f->is_valid();
    return(is_ok);
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
  string file;

  shared_ptr<NcFile> f;
  bool corrected;
  bool timeless;
  bool hourly;

  string calendar_type;
  int start_day;

  bool is_ok;
};

class DataRecord {
public:
  DataRecord(FileRecord& f, int offset, int month): f(f) {
    this->offset = offset;
    this->month = month;
    this->v = NULL;
  }
 DataRecord(FileRecord& f, NcVar* v, int offset, int month): f(f) {
    this->v = v;
    this->offset = offset;
    this->month = month;
  }
  int offset;
  int month;
  FileRecord& f;
  NcVar* v;

  // Should check more but this is adequate
  bool operator==(const DataRecord& r) {
    return (month == r.month && offset == r.offset);
  }

  bool operator<(const DataRecord& r) {
    return (month < r.month);
  }
};

