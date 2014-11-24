#include "common.h"

enum FNOFF{FMODEL,FEXPT,FVAR,FRUN,STARTYEAR,ENDYEAR};

int main(int argc, char** argv) {
  char buf[1024];
  list<FileRecord> l;

  // Try not to fall on your face, netcdf, when a dimension or variable is missing
  ncopts = NC_VERBOSE;

  if(argc < 3) {
    printf("Usage: extract_level <height> <output_path> [<output variable name>]\n");
  }

  bool climatology = false;
  string specified_outvar = "";
  double height = atof(argv[1]);
  string output_path = argv[2];
  vector<NcDim*> dims;
  vector<NcDim*> outdims;
  struct stat s;
  if(argc > 3) {
    specified_outvar = argv[3];
  }

  while(fgets(buf, 1024, stdin)) {
    // Chomp a la perl
    *(strchr(buf, '\n')) = '\0';
    FileRecord fr(buf, false);

    if(!fr.is_ok) {
      printf("Failed to open file %s\n", buf);
      continue;
    }

    // Parse out filename
    boost::char_separator<char> sep("-");
    string f = fr.file.substr(0, fr.file.rfind("."));
    tokenizer<char_separator<char> > tok(f, sep);
    vector<string> tokens(tok.begin(), tok.end());

    climatology = (tokens.size() >= ENDYEAR);

    // Format an output string (zg500, hur1000, etc)
    stringstream foo;
    if(specified_outvar.length() > 0) {
      foo << specified_outvar;
    } else {
      foo << fr.var << (height / 100);
    }
    const string output_var = foo.str();

    // Format output filename and create bits
    string path = output_path;
    string ofile = fr.model + "-" + fr.expt + "-" + output_var + "-" + fr.run;
    if(climatology) {
      ofile += "-" + tokens[STARTYEAR] + "-" + tokens[ENDYEAR];
    }
    ofile += ".nc";

    list<string>::const_iterator p;
    list<string> pathbits;
    pathbits.push_back(fr.expt);
    pathbits.push_back(output_var);
    pathbits.push_back(fr.model);
    pathbits.push_back(fr.run);

    for(p = pathbits.begin(); p != pathbits.end(); ++p) {
      path += "/" + *p;

      // If directory in path doesn't exist, create it.
      if(stat(path.c_str(), &s) != 0) {
	mkdir(path.c_str(), 0777);
      }
    }

    ofile = path + "/" + ofile;

    printf("Input file: %s\n", buf);
    printf("Output file: %s\n", ofile.c_str());

    NcFile out(ofile.c_str(), NcFile::Replace);
    NcFile& in = *(fr.f);

    assert(in.is_valid() && out.is_valid());
    
    out.set_fill(NcFile::NoFill);
    
    // Copy some dimensions (all except plev)
    for(int i = 0; i < in.num_dims(); i++) {
      NcDim* d = in.get_dim(i);
      assert(d);
      string name = d->name();
      if(name == "plev") 
	continue;
      assert(copy_dim(in, out, i));
    }

    // Copy all global attributes
    copy_atts(&in, &out);

    // Copy all vars verbatim except plev and fr.var
    int num_vars = in.num_vars();
    for(int i = 0; i < num_vars; i++) {
      NcVar* v = in.get_var(i);
      assert(v);
      string name = v->name();
      if(name == fr.var || name == "plev" || name == "lev")
	continue;
      copy_var(v, out);
    }

    // Scan plev for the pressure level we want (fail/return if we don't find it)
    NcVar* plev = in.get_var("plev");
    if(!plev) {
      plev = in.get_var("lev");
    }
    if(!plev) {
      printf("%s: Couldn't find plev layer\n", buf);
      continue;
    }
    long* edges = plev->edges();
    double* levels = new double[plev->num_vals()];
    int idx;
    assert(plev->get(levels, edges));
    for(idx = 0; idx < plev->num_vals(); idx++) {
      if(levels[idx] == height) {
	break;
      }
    }
    delete[] levels;
    delete[] edges;

    if(idx == plev->num_vals()) {
      printf("%s: Couldn't find height\n", buf);
      continue;
    }
    
    // Get the old var's dimension list
    dims.clear();
    outdims.clear();
    NcVar* invar = in.get_var(fr.var.c_str());
    assert(invar);
    populate_dimvec(invar, in, dims);

    // Get edges, modify the plev param in edges to 1, and create dim vector for new var
    edges = invar->edges();
    long cur[dims.size()];
    int timeidx = -1;
    for(unsigned int i = 0; i < dims.size(); i++) {
      string name = dims[i]->name();
      cur[i] = 0;
      if(name == "time") {
	edges[i] = 1;
	timeidx = i;
      }
      if(name == "plev" || name == "lev") {
	edges[i] = 1;
	cur[i] = idx;
      } else {
	outdims.push_back(out.get_dim(dims[i]->name()));
      }
    }

    assert(climatology || timeidx != -1);

    // Create the new variable
    NcVar *outvar = out.add_var(output_var.c_str(), invar->type(), outdims.size(), (const NcDim**)&outdims[0]);

    // Copy the attributes
    copy_atts(invar, outvar);

    assert(outvar);

    // Calculate the size of the data
    int size = 1;
    for(unsigned int i = 0; i < dims.size(); i++)
      size *= edges[i];

    float* dat = new float[size];
    long* outedges = outvar->edges();
    assert(get_recsize_and_edges(outvar, outedges) == size);
    if(climatology) {
      // Set the offset to the offset in plev of the pressure level
      invar->set_cur(cur);
      
      // Get the old data using the modified edge vector
      invar->get(dat, edges);
      
      // Copy the data to the new variable
      outvar->put(dat, outedges);
    } else {
      long outcur[outdims.size()];
      int outtimeidx = -1;
      for(unsigned int i = 0; i < outdims.size(); i++) {
	outcur[i] = 0;
	string name = dims[i]->name();
	if(name == "time") {
	  outtimeidx = i;
	}
      }
      assert(outtimeidx != -1);
      for(unsigned int i = 0; i < dims[timeidx]->size(); i++) {
	cur[timeidx] = i;
	outcur[outtimeidx] = i;

	// Set the offset to the correct time period and pressure level
	invar->set_cur(cur);
	outvar->set_cur(outcur);
	
	// Get the old data using the modified edge vector, and copy it into the new var
	invar->get(dat, edges);
	outvar->put(dat, outedges);
      }
    }

    delete[] dat;
    delete[] outedges;
    delete[] edges;
  }
}
