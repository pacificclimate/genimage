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
    default:
      assert(false);
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

  if(argc < 2) {
    printf("Usage: fix_missing <output_path>");
  }

  string output_path = argv[1];

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
	emit_and_cleanup(l, output_path);
      }
    }
    l.push_back(fr);
  }
  emit_and_cleanup(l, output_path);
}
