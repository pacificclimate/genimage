//#include <iostream>

#include "common.h"

void create_climatology(FileRecord& f, string outpath, const Range<int>& r, list<int>& omitlist) {
  int numtimes, start_offset;
  NcFile& in = *(f.f);
  assert(in.is_valid());
  int* times = 0;
  // Get time
  if(!f.timeless) {
    NcVar* intime = in.get_var("time");
    long* tedges = intime->edges();
    numtimes = intime->num_vals();
    times = new int[numtimes];
    intime->get(times, tedges);
    delete[] tedges;
    
    start_offset = do_binary_search(r.min, numtimes, times);
    if(start_offset == -1 || start_offset + (r.max - r.min + 1) >= numtimes) {
      printf("Specified offsets are not within range of data set\n");
      return;
    }
  }

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
  sst << "/" << f.model << "-" << f.expt << "-" << f.var << "-" << f.run << "-" << (r.min / 12) << "-" << (r.max / 12) << ".nc";

  string ofile = path + "/" + sst.str();

  NcFile out(ofile.c_str(), NcFile::Replace);
  
  printf("Output file: %s\n", ofile.c_str());

  if(!f.timeless)
    printf("Start: %i, End: %i, Max: %i\n", start_offset, start_offset + (r.max - r.min + 1), numtimes - 1);

  assert(out.is_valid());

  out.set_fill(NcFile::NoFill);

  copy_dims(in, out);
  copy_atts(&in, &out);
  NcDim* toy = out.add_dim("timeofyear", MAX_TOY);
  NcDim* lat = out.get_dim("lat");
  NcDim* lon = out.get_dim("lon");

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
  NcVar* outvar = out.add_var(invar->name(), invar->type(), toy, lat, lon);
  copy_atts(invar, outvar);
  long* edges = invar->edges();
  int rec_size = get_recsize_and_edges(invar, edges);
  int data_size = rec_size * MAX_TOY;
  float* data = new float[data_size];
  float* indata = new float[rec_size];
  if(f.timeless) {
    invar->get(indata, edges);
    // Just copy the damned thing
    for(int i = 0; i < MAX_TOY; i++) {
      for(int j = 0; j < rec_size; j++) {
	data[(i * rec_size) + j] = indata[j];
      }
    }
  } else {
    // Do the averaging
    int days[MAX_TOY];

    // Clear data
    for(int i = 0; i < data_size; i++) {
      data[i] = 0.0;
    }
    for(int i = 0; i < MAX_TOY; i++) {
      days[i] = 0;
    }

    list<int>::const_iterator omitted = omitlist.begin();
    for(int i = start_offset; i < numtimes && times[i] <= r.max; i++) {
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
  }
  
  delete[] edges;

  // Finally, throw the data into the output
  edges = outvar->edges();
  outvar->put(data, edges);

  delete[] edges;
  delete[] indata;
  delete[] data;
  if(times)
    delete[] times;
}

int main(int argc, char** argv) {
  char buf[1024];
  list<Range<int> > ranges;

  // Try not to fall on your face, netcdf, when a dimension or variable is missing
  ncopts = NC_VERBOSE;

  if(argc < 2) {
    printf("Usage: create_climatologies <output_path> [<climatology start> <climatology end>...]");
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
    list<int> omitlist;
    *(strchr(buf, '\n')) = '\0';
    string filename = buf;

    // Split up input line into filename and omit bits
    boost::char_separator<char> sep(",");
    tokenizer<char_separator<char> > tok(filename, sep);
    tokenizer<char_separator<char> >::const_iterator bits = tok.begin();
    filename = *bits;
    for(++bits; bits != tok.end(); ++bits) {
      
      omitlist.push_back(atoi((*bits).c_str()));
    }
    FileRecord fr(filename.c_str(), false);

    if(!fr.is_ok) {
      printf("Failed to open file %s\n", buf);
      continue;
    }

    // Generate climatologies
    list<Range<int> >::const_iterator i;
    for(i = ranges.begin(); i != ranges.end(); i++) {
      create_climatology(fr, output_path, *i, omitlist);
    }
  }
}
