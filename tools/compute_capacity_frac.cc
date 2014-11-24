#include "common.h"

enum FNOFF{FMODEL,FEXPT,FVAR,FRUN,STARTYEAR,ENDYEAR};

int main(int argc, char** argv) {
  char buf[1024];
  list<FileRecord> l;

  // Try not to fall on your face, netcdf, when a dimension or variable is missing
  ncopts = NC_VERBOSE;

  if(argc < 2) {
    printf("Usage: compute_capacity_frac <output path>\n");
    printf("\tInput: List of files, 2 per line: mrso on left, mrsofc on right\n");
  }

  struct stat s;
  const string specified_outvar = "mrsocf";
  const string output_path = argv[1];

  while(fgets(buf, 1024, stdin)) {
    // Chomp a la perl
    *(strchr(buf, '\n')) = '\0';

    // Split on comma
    char* filename = buf;
    char* cap_filename = strchr(buf, ',');
    *cap_filename = '\0';
    cap_filename++;
    
    // Create file records
    FileRecord cap_fr(cap_filename, false);
    FileRecord fr(filename, false);
    if(!cap_fr.is_ok) {
      printf("Failed to open file %s\n", buf);
      continue;
    }
    if(!fr.is_ok) {
      printf("Failed to open file %s\n", buf);
      continue;
    }
    const string output_var = "soil";

    // Format output filename and create bits
    string path = output_path;
    string ofile = fr.model + "-" + fr.expt + "-" + output_var + "-" + fr.run + ".nc";

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
    copy_dims(in, out);

    // Copy all global attributes
    copy_atts(&in, &out);

    // Copy all vars verbatim except fr.var
    int num_vars = in.num_vars();
    for(int i = 0; i < num_vars; i++) {
      NcVar* v = in.get_var(i);
      assert(v);
      string name = v->name();
      if(name == fr.var)
	continue;
      copy_var(v, out);
    }

    // Get mrsofc
    NcVar* mrsofc_var = cap_fr.f->get_var("mrsofc");
    assert(mrsofc_var);
    long* edges = mrsofc_var->edges();
    int recsize = get_recsize_and_edges(mrsofc_var, edges);
    float* mrsofc = new float[recsize];
    assert(mrsofc_var->get(mrsofc, edges));
    delete[] edges;
    
    // Create the new variable
    NcVar* invar = in.get_var(fr.var.c_str());
    NcVar* outvar = copy_var_structure(invar, out);
    outvar->rename(output_var.c_str());
    copy_atts(invar, outvar);
    assert(outvar);

    // Get ready to copy/modify data
    edges = outvar->edges();
    recsize = get_recsize_and_edges(outvar, edges);
    assert(recsize > 0);

    // Allocate memory
    float* dat = new float[recsize];
    
    // Get time dimesion so we know how far to go
    NcDim* time = out.get_dim("time");

    // Get missing value
    NcAtt* mv;
    float missing = 1e20f;
    if((mv = outvar->get_att("missing_value"))) {
      missing = mv->as_float(0);
    }

    for(unsigned int i = 0; i < time->size(); i++) {
      // Set the offset to the correct time period and pressure level
      assert(invar->set_cur(i));
      assert(outvar->set_cur(i));
      
      // Get the old data
      invar->get(dat, edges);

      // Process
      divide_grid_by_grid(recsize, dat, mrsofc, missing);
      for(int j = 0; j < recsize; j++) {
	if(dat[j] > 2) {
	  dat[j] = missing;
	}
      }

      // Put the new data
      outvar->put(dat, edges);
    }

    delete[] mrsofc;
    delete[] dat;
    delete[] edges;
  }
}
