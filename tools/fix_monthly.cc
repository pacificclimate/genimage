#include "common.h"

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
	printf("Duplicate record found; file %s and file %s\n", (*dit).f.filename.c_str(), (*old).f.filename.c_str());
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
    if(!f.open()) {
      printf("Failed to reopen file %s\n", f.filename.c_str());
      assert(false);
    }
    assert((v = f.f->get_var(f.var.c_str())));
    if(f.timeless) {
      DataRecord dr(*li, 0, 0);
      drlist.push_back(dr);
    } else {
      assert((t = f.f->get_var("time")));
      int i;
      int len = t->get_dim(0)->size();
      double* days = new double[len];
      t->get(days, len);

      if(f.hourly) 
	for(i = 0; i < len; i++) 
	  days[i] /= 24;

      for(i = 0; i < len; i++) {
	//printf("Start day: %i, Day: %f\n", f->start_day, days[i]);
	DataRecord dr(*li, i, get_total_months(f.calendar_type, (int)(f.start_day + days[i])));

	drlist.push_back(dr);
      }
      delete[] days;
    }
    assert(f.close());
  }
}

void emit_and_cleanup(list<FileRecord>& l, string outpath) {
  string path = outpath;
  struct stat s;
  list<DataRecord> drlist;

  FileRecord& f = *(l.begin());

  populate_drlist(l, drlist);
  remove_dups(drlist);

  string ofile = f.model + "-" + f.expt + "-" + f.var + "-" + f.run + ".nc";

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

  ofile = path + "/" + ofile;

  f.open();

  NcFile out(ofile.c_str(), NcFile::Replace);
  NcFile& in = *(f.f);

  assert(in.is_valid() && out.is_valid());

  out.set_fill(NcFile::NoFill);

  copy_dims(in, out);
  copy_atts(&in, &out);

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
  NcVar* outtime = 0;
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
  int recsize = get_recsize_and_edges(invar, edges);
  printf("recsize: %i, listsize: %li\n", recsize, drlist.size());
  
  // Run through the list of data bits
  float* fdata = new float[recsize];
  double* ddata = new double[recsize];
  list<DataRecord>::const_iterator ts;
  vector<int> months;

  f.close();

  // FIXME: Improve method of calculating output offset (use initial month?)
  int j = 0;

  for(ts = drlist.begin(); ts != drlist.end(); ++ts) {
    const DataRecord& d = *ts;
    NcVar* v;
    d.f.open();
    assert((v = d.f.f->get_var(f.var.c_str())));

    float bias = 0;
    float scale_factor = 1;
    
    NcAtt* var_add_offset = v->get_att("add_offset");
    NcAtt* var_scale_factor = v->get_att("scale_factor");

    if(var_add_offset) {
      bias += -var_add_offset->as_float(0);
    }

    if(var_scale_factor) {
      scale_factor *= var_scale_factor->as_float(0);
    }

    months.push_back(d.month);

    // Add the data to the output variable
    assert(v->set_cur(d.offset));
    assert(outvar->set_cur(j));

    // Add the data to the output variable
    switch(v->type()) {
    case ncByte:
    case ncChar:
    case ncShort:
      copy_to_float(v, outvar, recsize, edges, (short*)ddata, fdata, scale_factor, bias);
      break;
    case ncLong:
      copy_to_float(v, outvar, recsize, edges, (int*)ddata, fdata, scale_factor, bias);
      break;
    case ncFloat:
      copy_float_to_float(v, outvar, recsize, edges, fdata);
      break;
    case ncDouble:
      copy_double_to_float(v, outvar, recsize, edges, ddata, fdata);
      break;
    default:
      assert(false);
    }
    d.f.close();
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
  bool filename_date_parsing = false;

  // Try not to fall on your face, netcdf, when a dimension or variable is missing
  NcError n(NcError::silent_nonfatal);

  if(argc < 2) {
    printf("Usage: fix_missing <output_path> [<filename date parsing boolean>]\n");
  }

  if(argc > 2) {
    filename_date_parsing = (atoi(argv[2]) > 0);
  }

  string output_path = argv[1];

  while(fgets(buf, 1024, stdin)) {
    // Chomp a la perl
    *(strchr(buf, '\n')) = '\0';
    FileRecord fr(buf, true, filename_date_parsing);

    if(!fr.is_ok) {
      printf("Failed to open file %s\n", buf);
      continue;
    }

    if(l.size()) {
      const FileRecord& oldfr = l.back();
      if(fr != oldfr) {
	emit_and_cleanup(l, output_path);
      }
    }
    assert(fr.close());
    l.push_back(fr);
  }
  emit_and_cleanup(l, output_path);
}
