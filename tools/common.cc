#include "common.h"

vector<int> get_time_values(NcFile& f, int start_day) {
  NcVar* time_var = f.get_var("time");
  assert(time_var);

  switch(time_var->type()) {
  case ncShort:
    return get_time_values<short>(time_var, start_day);
    break;
  case ncLong:
    return get_time_values<long>(time_var, start_day);
    break;
  case ncFloat:
    return get_time_values<float>(time_var, start_day);
    break;
  case ncDouble:
    return get_time_values<double>(time_var, start_day);
    break;
  default:
    assert(false);
    {
      vector<int> vi;
      return vi;
    }
  }
}

double get_missing_value(NcVar* v) {
  double miss = 1e20;
  NcAtt* fillvalue_att = v->get_att("_FillValue");
  NcAtt* missing_att = v->get_att("missing_value");
  if(fillvalue_att) {
    miss = fillvalue_att->as_double(0);
    delete fillvalue_att;
  } else if(missing_att) {
    miss = missing_att->as_double(0);
    delete missing_att;
  }
  return(miss);
}

float get_missing_value_float(NcVar* v) {
  float miss = 1e20;
  NcAtt* fillvalue_att = v->get_att("_FillValue");
  NcAtt* missing_att = v->get_att("missing_value");
  if(fillvalue_att) {
    miss = fillvalue_att->as_float(0);
    delete fillvalue_att;
  } else if(missing_att) {
    miss = missing_att->as_float(0);
    delete missing_att;
  }
  if(missing_att)
    delete missing_att;
  return(miss);
}

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
  time_t seconds_since_1970 = (time_t)(days - DAYS_1970) * 86400;
  struct tm td;
  assert(gmtime_r(&seconds_since_1970, &td));

  //fprintf(stderr, "%i %i %i %i\n", td.tm_year, td.tm_mon, td.tm_mday, td.tm_hour);
  
  return 12 * (td.tm_year + 1900) + td.tm_mon;
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
  if(calendar_type == "gregorian" || calendar_type == "proleptic_gregorian") {
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
  struct tm intime = {
    0, 0, 0,
    start_day, start_month - 1, start_year - 1900,
    0, 0, 0,
    0, NULL
  };

  // Ensure that UTC is what we'll get...
  time_t time_secs_since_1970 = mktime(&intime);
  assert(time_secs_since_1970 != -1);
  
  return((int)(time_secs_since_1970 / 86400) + DAYS_1970);
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
  if(calendar_type == "gregorian" || calendar_type == "proleptic_gregorian") {
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

void populate_dimvec(NcVar* src, NcFile& dst, vector<NcDim*>& dims) {
  int num_dims = src->num_dims();
  for(int i = 0; i < num_dims; i++) {
    const string dimname = src->get_dim(i)->name();
    NcDim* d = dst.get_dim(dimname.c_str());
    assert(d);
    dims.push_back(d);
  }
}

void populate_dimvec(NcVar* src, vector<NcDim*>& dims) {
  int num_dims = src->num_dims();
  for(int i = 0; i < num_dims; i++) {
    NcDim* d = src->get_dim(i);
    assert(d);
    dims.push_back(d);
  }
}

NcVar* copy_var_structure(NcVar* src, NcFile& dst) {
  // Get dimensions to use
  vector<NcDim*> dims;
  populate_dimvec(src, dst, dims);
  
  // Create output variable
  return dst.add_var(src->name(), src->type(), dims.size(), (const NcDim**)&dims[0]);
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

void add_to_grid(int size, float* input, float* accum) {
  for(int i = 0; i < size; i++) {
    accum[i] += input[i];
  }
}

void add_to_grid(int size, float* input, float* accum, float missing) {
  for(int i = 0; i < size; i++) {
    if(accum[i] != missing) {
      if(input[i] != missing) {
	accum[i] += input[i];
      } else {
	accum[i] = missing;
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
      if(input != missing) {
	accum[i] += input;
      } else {
	accum[i] = missing;
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
