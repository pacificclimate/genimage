//#include <iostream>

#include "common.h"
#include <boost/lexical_cast.hpp>

vector<int> get_time_as_months_since_0000(const FileRecord& f, NcFile& in) {
  NcVar* intime = in.get_var("time");
  long* tedges = intime->edges();
  int numtimes = intime->num_vals();
  vector<int> times(numtimes);

  NcAtt* time_units_att = intime->get_att("units");
  assert(time_units_att);
  char* tu_temp = time_units_att->as_string(0);
  string time_units = tu_temp;
  delete[] tu_temp;
  delete time_units_att;
  
  // Handle climatologies with units of "days since x"
  if(time_units.compare(0, 4, "days") == 0) {
    vector<int> time_days = get_time_values(in, f.start_day);
    for(int i = 0; i < numtimes; ++i)
      times[i] = get_total_months(f.calendar_type, time_days[i]);
  } else {
    // Months...
    intime->get(&times[0], tedges);
  }
  delete[] tedges;
  return times;
}

string create_output_filename_make_path(const FileRecord& f, const Range<int>& r, const string outpath) {
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
      assert(mkdir(path.c_str(), 0777) == 0); 
    }
  }

  stringstream sst;
  sst << "/" << f.model << "-" << f.expt << "-" << f.var << "-" << f.run << "-" << (r.min / 12) << "-" << (r.max / 12) << ".nc";

  string ofile = path + "/" + sst.str();
  return ofile;
}

void copy_all_but_time_and_var(const FileRecord& f, NcFile& in, NcFile& out) {
  copy_dims(in, out);
  copy_atts(&in, &out);
  
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
}

void compute_climatology(NcVar* invar, const float missing, float* data, const vector<int>& times, const int start_offset, const int rec_size, const long* edges, const string calendar_type, const list<int>& omitlist, const Range<int>& r) {
  float* indata = new float[rec_size];
  int days[MAX_TOY];
  std::fill(&data[0], &data[rec_size * MAX_TOY], 0.0f);
  std::fill(&days[0], &days[MAX_TOY], 0);

  list<int>::const_iterator omitted = omitlist.begin();
  for(unsigned int i = start_offset; i < times.size() && times[i] <= r.max; i++) {
    assert(invar->set_cur(i));
    assert(invar->get(indata, edges));
    int days_in_month;
    int year = times[i] / 12;
    int month = times[i] % 12;
    
    if(omitted != omitlist.end() && times[i] == *omitted) {
      ++omitted;
      continue;
    }

    // Get the right # of days for this month
    if(calendar_type == "gregorian") {
      if(year % 4 == 0 && (year % 100 != 0 || year % 400 == 0)) {
	// leap year
	days_in_month = leap_dpm[month];
      } else {
	// not leap year
	days_in_month = noleap_dpm[month];
      }
    } else if(calendar_type == "365_day") {
      days_in_month = noleap_dpm[month];
    } else if(calendar_type == "360_day") {
      days_in_month = equal_dpm[month];
    } else {
      days_in_month = noleap_dpm[month];
    }

    // Multiply by # of days to give accumulated days of the mean in the month
    multiply_grid_by_scalar(rec_size, indata, (float)days_in_month, missing);
    
    // Do month
    add_to_grid(rec_size, indata, data + (rec_size * month), missing);
    days[month] += days_in_month;
  }
  
  // Divide grids by # days
  int ann_days = 0;
  for(int i = JAN; i <= DEC; i++) {
    add_to_grid(rec_size, data + (rec_size * i), data + (rec_size * ANN), missing);
    ann_days += days[i];
  }
  add_to_grid(rec_size, data + (rec_size * DEC), data + (rec_size * DJF), missing);
  add_to_grid(rec_size, data + (rec_size * JAN), data + (rec_size * DJF), missing);
  add_to_grid(rec_size, data + (rec_size * FEB), data + (rec_size * DJF), missing);
  add_to_grid(rec_size, data + (rec_size * MAR), data + (rec_size * MAM), missing);
  add_to_grid(rec_size, data + (rec_size * APR), data + (rec_size * MAM), missing);
  add_to_grid(rec_size, data + (rec_size * MAY), data + (rec_size * MAM), missing);
  add_to_grid(rec_size, data + (rec_size * JUN), data + (rec_size * JJA), missing);
  add_to_grid(rec_size, data + (rec_size * JUL), data + (rec_size * JJA), missing);
  add_to_grid(rec_size, data + (rec_size * AUG), data + (rec_size * JJA), missing);
  add_to_grid(rec_size, data + (rec_size * SEP), data + (rec_size * SON), missing);
  add_to_grid(rec_size, data + (rec_size * OCT), data + (rec_size * SON), missing);
  add_to_grid(rec_size, data + (rec_size * NOV), data + (rec_size * SON), missing);
  
  divide_grid_by_scalar(rec_size, data + (rec_size * DJF), (float)(days[DEC] + days[JAN] + days[FEB]), missing);
  divide_grid_by_scalar(rec_size, data + (rec_size * MAM), (float)(days[MAR] + days[APR] + days[MAY]), missing);
  divide_grid_by_scalar(rec_size, data + (rec_size * JJA), (float)(days[JUN] + days[JUL] + days[AUG]), missing);
  divide_grid_by_scalar(rec_size, data + (rec_size * SON), (float)(days[SEP] + days[OCT] + days[NOV]), missing);
  divide_grid_by_scalar(rec_size, data + (rec_size * ANN), (float)ann_days, missing);
  
  for(int i = JAN; i <= DEC; i++) {
    divide_grid_by_scalar(rec_size, data + (rec_size * i), (float)days[i], missing);
  }
  delete[] indata;
}

void add_cf_compliant_time_dim(NcFile& out, const FileRecord& f, NcVar* outvar, const Range<int>& r) {
  // Add valid time information here as per CF 1.6
  NcVar* time_var = out.get_var("time");
  assert(time_var);
  time_var->add_att("calendar", f.calendar_type.c_str());
  int total_months_start_day = get_total_months(f.calendar_type, f.start_day);
  char* buf = new char[100];
  sprintf(buf, "days since %i-%i-%i", total_months_start_day / 12, total_months_start_day % 12 + 1, 1);
  time_var->add_att("units", buf);
  delete[] buf;

  string clim_bnds_name = "climatology_bounds";
  string bnds_name = "bnds";
  // FIXME: This should reflect the actual cell_methods if there are any.
  outvar->add_att("cell_methods", "time: mean within days time: mean over years");
  time_var->add_att("climatology", clim_bnds_name.c_str());
  NcDim* time_dim = out.get_dim("time");
  NcDim* bnds = 0;
  for(int i = 0; i < out.num_dims(); ++i) {
    if(bnds_name == out.get_dim(i)->name()) {
      bnds = out.get_dim(i);
      break;
    }
  }
  if(!bnds)
    bnds = out.add_dim(bnds_name.c_str(), 2);
  
  NcVar* climatology_bounds = out.add_var(clim_bnds_name.c_str(), ncDouble, time_dim, bnds);
  
  int start_year = r.min / 12;
  int end_year = r.max / 12;
  int mid_year = (start_year + end_year) / 2;
  //                    Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec DJF MAM JJA SON ANN
  int start_month[] = { 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 12,  3,  6,  9,  1 };
  int mid_month[]   = { 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,  1,  4,  7, 10,  7 };
  int end_month[]   = { 2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,  1,  3,  6,  9, 12,  1 };
  int mid_day[]     = {16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 15, 15, 15, 15,  1 };
  
  int time_days[17];
  int time_bnds[17 * 2];

  for(int i = 0; i < 17; ++i) {
    time_days[i] = date2days(f.calendar_type, mid_year, mid_month[i], mid_day[i]) - f.start_day;
    time_bnds[i * 2] = date2days(f.calendar_type, start_year, start_month[i], 1) - f.start_day;
    time_bnds[i * 2 + 1] = date2days(f.calendar_type, (i == 17 ? end_year + 1 : end_year), end_month[i], 1) - f.start_day;
  }
  
  time_var->put(time_days, 17);
  climatology_bounds->put(time_bnds, 17, 2);
}

void create_climatology(FileRecord& f, string outpath, const Range<int>& r, list<int>& omitlist) {
  int start_offset = 0;
  vector<int> times;
  NcFile& in = *(f.f);
  assert(in.is_valid());
  // Get time
  if(!f.timeless) {
    times = get_time_as_months_since_0000(f, in);

    start_offset = do_binary_search(r.min, times.size(), &times[0]);
    printf("Start offset: %i\n", start_offset);
    if(start_offset == -1 || (unsigned int)(start_offset + (r.max - r.min)) >= times.size()) {
      printf("Specified offsets are not within range of data set: %i to %i, start is %i, max is %i\n", r.min, r.max, times[0], times[times.size() - 1]);
      return;
    }
  }

  string ofile = create_output_filename_make_path(f, r, outpath);
  NcFile out(ofile.c_str(), NcFile::Replace);
  assert(out.is_valid());
  printf("Output file: %s\n", ofile.c_str());
  out.set_fill(NcFile::NoFill);

  if(!f.timeless)
    printf("Start: %i, End: %i, Max: %i\n", start_offset, start_offset + (r.max - r.min + 1), (int)times.size() - 1);

  copy_all_but_time_and_var(f, in, out);

  NcVar* invar = in.get_var(f.var.c_str());

  NcVar* outvar = 0;
  long* edges = invar->edges();
  int rec_size = get_recsize_and_edges(invar, edges);
  float* data = new float[rec_size * MAX_TOY];

  if(f.timeless) {
    float* indata = new float[rec_size];
    NcDim* lat = out.get_dim("lat");
    NcDim* lon = out.get_dim("lon");
    NcDim* toy = out.add_dim("time");
    NcVar* time = out.add_var("time", ncDouble, toy);
    assert(time);
    f.calendar_type = "gregorian";
    f.start_day = date2days(f.calendar_type, r.min / 12, r.min % 12 + 1, 1);
    outvar = out.add_var(invar->name(), invar->type(), toy, lat, lon);
    copy_atts(invar, outvar);
    invar->get(indata, edges);
    // Just copy the damned thing
    for(int i = 0; i < MAX_TOY; i++)
      std::copy(indata, &indata[rec_size], &data[i * rec_size]);
    delete[] indata;
  } else {
    NcVar* time = out.add_var("time", ncDouble, out.get_dim("time"));
    assert(time);

    vector<NcDim*> dims;
    populate_dimvec(invar, out, dims);
    
    outvar = out.add_var(invar->name(), invar->type(), dims.size(), (const NcDim**)&dims[0]);
    copy_atts(invar, outvar);

    float missing = get_missing_value_float(outvar);

    compute_climatology(invar, missing, data, times, start_offset, rec_size, edges, f.calendar_type, omitlist, r);
  }

  add_cf_compliant_time_dim(out, f, outvar, r);

  // Finally, throw the data into the output
  long* out_edges = outvar->edges();
  outvar->put(data, out_edges);

  delete[] out_edges;
  delete[] edges;
  delete[] data;
}

int main(int argc, char** argv) {
  char buf[10240];
  bool cmip5_paths = false;

  // Try not to fall on your face, netcdf, when a dimension or variable is missing
  NcError n(NcError::silent_nonfatal);

  if(argc < 2) {
    printf("Usage: create_climatologies <output_path> [<cmip5_paths>]\n");
    exit(1);
  }

  string output_path = argv[1];

  if(argc == 3)
    cmip5_paths = true;

  // Second step: Read the files in, one by one; extract the climatologies
  while(fgets(buf, 10240, stdin)) {
    // Chomp a la perl
    *(strchr(buf, '\n')) = '\0';

    if(*buf == '#') {
      printf("Comment!\n");
      continue;
    }

    list<int> omitlist;
    list<Range<int> > ranges;
    list<Range<int> >::const_iterator i;
    string input_line = buf;

    // Split up input line into filename and omit bits
    boost::char_separator<char> commasep(",");
    boost::char_separator<char> colonsep(":");
    tokenizer<char_separator<char> > tok(input_line, commasep);
    tokenizer<char_separator<char> >::const_iterator bits = tok.begin();

    assert(tok.begin() != tok.end());

    string filename = *bits;
    
    for(++bits; bits != tok.end(); ++bits) {
      if((*bits).find(":") != string::npos) {
	tokenizer<char_separator<char> > rtok(*bits, colonsep);
	tokenizer<char_separator<char> >::const_iterator rbits = rtok.begin();
	int min, max;
	try {
	  min = lexical_cast<int>(*rbits);
	  ++rbits;
	  assert(rbits != rtok.end());
	  max = lexical_cast<int>(*rbits);
	  Range<int> r(min, max);
	  ranges.push_back(r);
	} catch(bad_lexical_cast& b) {
	  assert(false);
	}
      } else {
	omitlist.push_back(atoi((*bits).c_str()));
      }
    }
    FileRecord fr(filename.c_str(), true, false, cmip5_paths);

    if(!fr.is_ok) {
      printf("Failed to open file %s\n", buf);
      continue;
    }

    // Generate climatologies
    for(i = ranges.begin(); i != ranges.end(); i++) {
      create_climatology(fr, output_path, *i, omitlist);
    }
  }
}
