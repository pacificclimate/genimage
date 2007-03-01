//#include <iostream>
#include <netcdfcpp.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <list>
#include <vector>

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

enum FILEBITS{DOT,EXPT,VAR,MODEL,RUN,FILENAME};
const int noleap_days[] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 };
const int leap_days[] = { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 };
const int equal_days[] = { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360 };

int do_binary_search(int number, int max, const int array[]) {
  int min = 0;
  int offset = (min + max) / 2;
  max--;
  if(number < array[0] || number > array[max]) {
    return -1;
  }
  while(array[offset] > number || array[offset + 1] < number) {
    if(array[offset] < number) {
      min = (int)ceil((double)(max + min) / 2.0);
    } else {
      max = (int)floor((double)(max + min) / 2.0);
    }
    offset = (max + min) / 2;
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

void remove_dups(list<DataRecord>& drlist) {
  drlist.sort();
  list<DataRecord>::iterator dit, old;
  dit = drlist.begin();
  ++dit;
  old = drlist.begin();
  for(; dit != drlist.end();) {
    if((*old).month == (*dit).month) {
      if((*old).f.corrected && !(*dit).f.corrected) {
	//printf("Corrected record found: 1!\n");
	dit = drlist.erase(dit);
      } else if(!(*old).f.corrected && (*dit).f.corrected) {
	//printf("Corrected record found: 2!\n");
	old = drlist.erase(old);
	++dit;
      } else {
	// Duplicate record without one being a "corrected" record. Error.
	printf("Duplicate record found!\n");
	assert(0 && "Duplicate record error");
      }
      continue;
    }
    ++dit;
    ++old;
  }
}

void populate_drlist(list<FileRecord>& l, list<DataRecord>& drlist) {
  list<FileRecord>::iterator li;
  for(li = l.begin(); li != l.end(); li++) {
    FileRecord& f = *li;
    NcVar* t;
    NcVar* v;
    assert((v = f.f->get_var(f.var.c_str())));
    if(f.timeless) {
      DataRecord dr(*li, v, 0, 0);
      drlist.push_back(dr);
    } else {
      assert((t = f.f->get_var("time")));
      int i;
      int len = t->get_dim(0)->size();
      double* days = new double[len];
      t->get(days, len);

      for(i = 0; i < len; i++) {
	//printf("Start day: %i, Day: %f\n", f->start_day, days[i]);
	DataRecord dr(*li, v, i, get_total_months(f.calendar_type, (int)(f.start_day + days[i])));

	drlist.push_back(dr);
      }
      delete[] days;
    }
  }
}

void copy_dims(NcFile& src, NcFile& dst) {
  // Copy the dimensions
  for(int i = 0; i < src.num_dims(); i++) {
    NcDim* d = src.get_dim(i);
    NcDim* out;
    if(d->is_unlimited()) {
      out = dst.add_dim(d->name());
    } else {
      out = dst.add_dim(d->name(), d->size());
    }
    assert(out);
  }
}

void copy_double_to_float(NcVar* src, NcVar* dst, int bufsize, long* edges, double* indata, float* outdata) {
    assert(src->get(indata, edges));
    for(int i = 0; i < bufsize; i++) {
      outdata[i] = (float)indata[i];
    }
    assert(dst->put(outdata, edges));
}

void copy_float_to_float(NcVar* src, NcVar* dst, int bufsize, long* edges, float* buf) {
    assert(src->get(buf, edges));
    assert(dst->put(buf, edges));
}

void emit_and_cleanup(list<FileRecord>& l) {
  string ofile;

  list<DataRecord> drlist;
  populate_drlist(l, drlist);
  remove_dups(drlist);

  FileRecord& f = *(l.begin());
  ofile = f.expt + "/" + f.var + "/" + f.model + "/" + f.run + "/" + f.model + "-" + f.expt + "-" + f.var + "-" + f.run + ".nc";

  NcFile out(ofile.c_str(), NcFile::Replace);
  out.set_fill(NcFile::NoFill);
  NcFile& in = *(f.f);

  assert(in.is_valid() && out.is_valid());

  copy_dims(in, out);

  // Copy some variables without any modifications
  int num_vars = in.num_vars();
  for(int i = 0; i < num_vars; i++) {
    NcVar* v = in.get_var(i);
    assert(v);
    string name = v->name();
    if(name == "time" || name == "time_bnds" || name == f.var)
      continue;
    copy_var(v, out);
  }

  // Special handling of time
  vector<NcDim*> dims;
  NcVar* intime = in.get_var("time");
  NcVar* outtime;
  if(intime) {
    populate_dimvec(intime, out, dims);
    outtime = out.add_var("time", ncInt, dims.size(), (const NcDim**)&dims[0]);
    assert(outtime);
    copy_atts(intime, outtime);
    NcAtt* timeatt = outtime->get_att("units");
    assert(timeatt);
    timeatt->rename("old_units");
    delete timeatt;
    outtime->add_att("units", "months since 0000-01");
  }

  // Create output variable and copy attributes
  dims.clear();
  NcVar* invar = in.get_var(f.var.c_str());
  assert(invar);
  populate_dimvec(invar, out, dims);
  NcVar* outvar = out.add_var(f.var.c_str(), ncFloat, dims.size(), (const NcDim**)&dims[0]);
  assert(outvar);
  copy_atts(invar, outvar);

  // Construct the list of what to copy
  long* edges = invar->edges();
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
  printf("recsize: %i, listsize: %i\n", recsize, drlist.size());
  
  // Run through the list of data bits
  float* fdata = new float[recsize];
  double* ddata = new double[recsize];
  list<DataRecord>::const_iterator ts;
  vector<int> months;

  // FIXME: Improve method of calculating output offset (use initial month?)
  int j = 0;
  for(ts = drlist.begin(); ts != drlist.end(); ++ts) {
    const DataRecord& d = *ts;
    months.push_back(d.month);

    // Add the data to the output variable
    assert(d.v->set_cur(d.offset));
    assert(outvar->set_cur(j));

    // Add the data to the output variable
    switch(d.v->type()) {
    case ncByte:
    case ncChar:
    case ncShort:
    case ncLong:
      assert(false);
      break;
    case ncFloat:
      copy_float_to_float(d.v, outvar, recsize, edges, fdata);
      break;
    case ncDouble:
      copy_double_to_float(d.v, outvar, recsize, edges, ddata, fdata);
      break;
    }
    j++;
  }
  delete[] fdata;
  delete[] ddata;
  delete[] edges;
  
  // Copy months
  if(intime) {
    assert(outtime->put(&months[0], months.size()));
  }

  printf("Output file: %s\n", ofile.c_str());
  l.clear();
  drlist.clear();
}

int main(int argc, char** argv) {
  char buf[1024];
  list<FileRecord> l;

  // Try not to fall on your face, netcdf, when a dimension or variable is missing
  ncopts = NC_VERBOSE;

  while(fgets(buf, 1024, stdin)) {
    // Chomp a la perl
    *(strchr(buf, '\n')) = '\0';
    FileRecord fr(buf);

    if(!fr.is_ok) {
      printf("Failed to open file %s\n", buf);
      continue;
    }

    if(l.size()) {
      const FileRecord& oldfr = l.back();
      if(fr != oldfr) {
	emit_and_cleanup(l);
      }
    }
    l.push_back(fr);
  }
  emit_and_cleanup(l);
}
