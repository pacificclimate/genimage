#include "common.h"
#include "../core/derived_vars.h"
#include <map>

#define NUM_COLS 3

enum COLNAMES { DV_NAME=0, TEMP_FILE, PREC_FILE };

/* Expected input format: 3 columns of text, comma separated, on stdin.
    - Output variable name.
    - NetCDF file containing temperature data.
    - NetCDF file containing precipitation data.
 */

string get_filename_create_path(string output_path, FileRecord& temp, string new_var_name, string temp_varname="tas") {
  // Format output filename and create bits
  struct stat s;
  string path = output_path;
  string ofile = "";
  list<string>::const_iterator p;
  list<string> pathbits;
  pathbits.push_back(temp.expt);
  pathbits.push_back(new_var_name);
  pathbits.push_back(temp.model);
  pathbits.push_back(temp.run);
  
  for(p = pathbits.begin(); p != pathbits.end(); ++p) {
    path += "/" + *p;
    
    // If directory in path doesn't exist, create it.
    if(stat(path.c_str(), &s) != 0) {
      mkdir(path.c_str(), 0777);
    }
  }
  
  size_t varname_offset = temp.file.find(temp_varname);
  assert(varname_offset != string::npos);

  ofile = path + "/" + temp.file.replace(varname_offset, temp_varname.length(), new_var_name);
  return(ofile);
}

void copy_all_but_named_var(NcFile& in, NcFile& out, string named_var) {
  // Copy some dimensions (all except plev)
  copy_dims(in, out);
  
  // Copy all global attributes
  copy_atts(&in, &out);
  
  // Copy all vars verbatim except the main variable in the file.
  const int num_vars = in.num_vars();
  for(int i = 0; i < num_vars; i++) {
    NcVar* v = in.get_var(i);
    assert(v);
    string name = v->name();
    if(name == named_var)
      continue;
    copy_var(v, out);
  }
}

bool var_edges_same(NcVar* v1, NcVar* v2) {
    assert(v1->num_dims() == v2->num_dims());
    
    long* v1_edges = v1->edges();
    long* v2_edges = v2->edges();

    assert(v1_edges && v2_edges);

    for(int i = 0; i < v1->num_dims(); ++i) {
      if(v1_edges[i] != v2_edges[i]) {
	delete[] v1_edges;
	delete[] v2_edges;
	return(false);
      } 
    }

    delete[] v1_edges;
    delete[] v2_edges;
    return(true);
}


void add_derived_var_to_grid(NcVar* temp_var, NcVar* prec_var, string derived_var_name, int mon, double temp_missing, double prec_missing, double out_missing, double* out_grid, unsigned int grid_size, long starts[3], long counts[3]) {
  double* prec_grid = new double[grid_size];
  double* temp_grid = new double[grid_size];

  assert(temp_var->set_cur(starts));
  
  assert(temp_var->get(temp_grid, counts));
  add_scalar_to_grid(grid_size, temp_grid, -273.15, temp_missing);
  
  if(derived_var_name == "dl00") {
    // Degree-days less than 0
    degree_days_pcic(temp_grid, out_grid, temp_missing, out_missing, 0, LT, mon, grid_size);
  } else if(derived_var_name == "dg05") {
    // Degree-days greater than 5
    degree_days_pcic(temp_grid, out_grid, temp_missing, out_missing, 5, GT, mon, grid_size);
  } else if(derived_var_name == "dl18") {
    // Degree-days less than 18
    degree_days_pcic(temp_grid, out_grid, temp_missing, out_missing, 18, LT, mon, grid_size);
  } else if(derived_var_name == "dg18") {
    // Degree-days greater than 18
    degree_days_pcic(temp_grid, out_grid, temp_missing, out_missing, 18, GT, mon, grid_size);
  } else if(derived_var_name == "nffd") {
    // Number of frost-free days
    number_frost_free_days(temp_grid, out_grid, temp_missing, out_missing, mon, grid_size);
  } else if(derived_var_name == "pass") {
    // Prec as snow
    assert(prec_var->set_cur(starts));
    assert(prec_var->get(prec_grid, counts));
    multiply_grid_by_scalar(grid_size, prec_grid, 86400.0, prec_missing);
    precip_as_snow(temp_grid, prec_grid, out_grid, temp_missing, prec_missing, out_missing, mon, grid_size);
  }

  delete[] prec_grid;
  delete[] temp_grid;
}

void add_derived_var_attributes(NcVar* outvar, string derived_var_name) {
  if(derived_var_name == "dl00") {
    // Degree-days less than 0
    outvar->add_att("standard_name", "degree_days_lt_0");
    outvar->add_att("long_name", "Degree-days less than 0 degrees C");
    outvar->add_att("units", "degree-days");
  } else if(derived_var_name == "dg05") {
    // Degree-days greater than 5
    outvar->add_att("standard_name", "degree_days_gt_5");
    outvar->add_att("long_name", "Degree-days greater than 5 degrees C");
    outvar->add_att("units", "degree-days");
  } else if(derived_var_name == "dl18") {
    // Degree-days less than 18
    outvar->add_att("standard_name", "degree_days_lt_18");
    outvar->add_att("long_name", "Degree-days less than 18 degrees C");
    outvar->add_att("units", "degree-days");
  } else if(derived_var_name == "dg18") {
    // Degree-days greater than 18
    outvar->add_att("standard_name", "degree_days_gt_18");
    outvar->add_att("long_name", "Degree-days greater than 18 degrees C");
    outvar->add_att("units", "degree-days");
  } else if(derived_var_name == "nffd") {
    // Number of frost-free days
    outvar->add_att("standard_name", "number_frost_free_days");
    outvar->add_att("long_name", "Number of frost-free days");
    outvar->add_att("units", "days");
  } else if(derived_var_name == "pass") {
    // Prec as snow
    outvar->add_att("standard_name", "snowfall_amount");
    outvar->add_att("long_name", "Snowfall amount");
    outvar->add_att("units", "kg m-2");
  }
}

int get_month(string calendar_type, int time_val) {
  return get_total_months(calendar_type, time_val) % 12;
}


int main(int argc, char** argv) {
  list<FileRecord> l;
  const string temp_var_name = "tas", prec_var_name = "pr";
  unsetenv("TZ");

  NcError n(NcError::silent_nonfatal);

  if(argc < 2) {
    printf("Usage: compute_derived_vars <output path>\n");
    printf("\tInput: List of stuff: derived_var, temp_file, prec_file\n");
  }

  const string output_path = argv[1];

  std::string filedata((std::istreambuf_iterator<char>(cin)), std::istreambuf_iterator<char>());
  
  boost::tokenizer<boost::escaped_list_separator<char> > tok(filedata, boost::escaped_list_separator<char>("\\", "\n", "\""));
  boost::tokenizer<boost::escaped_list_separator<char> >::iterator beg;

  for(beg = tok.begin(); beg != tok.end(); ++beg) {
    // Split up the comma separated data into individual sub-components
    vector<string> bits(NUM_COLS);
    boost::char_separator<char> commasep(",");
    tokenizer<char_separator<char> > tok_line(*beg, commasep);
    if(tok_line.begin() == tok_line.end()) 
      break;

    tokenizer<char_separator<char> >::const_iterator tok_line_iter = tok_line.begin();
    for(int col = 0; col < NUM_COLS && tok_line_iter != tok_line.end(); ++col, ++tok_line_iter) {
      bits[col] = *tok_line_iter;
    }

    cout << "Variable: " << bits[DV_NAME] << ", Temp file: " << bits[TEMP_FILE] << ", Prec file: " << bits[PREC_FILE] << endl;

    // Open files, make sure they're valid.
    FileRecord temp(bits[TEMP_FILE]);
    FileRecord prec(bits[PREC_FILE]);
    assert(temp.is_ok && prec.is_ok);

    // Get variables, make sure they exist.
    NcVar* temp_var = temp.f->get_var(temp_var_name.c_str());
    NcVar* prec_var = prec.f->get_var(prec_var_name.c_str());
    assert(temp_var);
    assert(prec_var);

    double temp_missing = get_missing_value(temp_var);
    double prec_missing = get_missing_value(prec_var);
    double out_missing = 1e20;

    // Test that dimensions are same for prec and temp...
    assert(var_edges_same(temp_var, prec_var));

    // Create and set up output file.
    string ofile = get_filename_create_path(output_path, temp, bits[DV_NAME]);
    NcFile out(ofile.c_str(), NcFile::Replace);
    assert(out.is_valid());
    out.set_fill(NcFile::NoFill);
    copy_all_but_named_var(*(temp.f), out, temp_var_name);
    vector<NcDim*> dims;
    populate_dimvec(temp_var, out, dims);
    NcVar* outvar = out.add_var(bits[DV_NAME].c_str(), temp_var->type(), dims.size(), (const NcDim**)&dims[0]);
    outvar->add_att("missing_value", out_missing);
    outvar->add_att("_FillValue", out_missing);

    // Set up dimensions for use in reading and writing files.
    long* edges = outvar->edges();
    long counts[3] = { 1, dims[1]->size(), dims[2]->size() };
    long starts[3] = { 0, 0, 0 };
    cout << "Dim sizes: " << dims[0]->size() << " " << dims[1]->size() << " " << dims[2]->size() << endl;
    int grid_size = dims[2]->size() * dims[1]->size();

    // Prepare temporary results storage...
    double* out_grid = new double[grid_size];

    // Add required metadata attributes.
    add_derived_var_attributes(outvar, bits[DV_NAME]);

    cout << "Calendar: " << temp.calendar_type << endl;

    // Only process data with an unlimited dim.
    if(dims[0]->is_unlimited()) {
      vector<int> time_vals = get_time_values(*(temp.f), temp.start_day);
      NcVar* time_var = out.get_var("time");
      NcAtt* climatology_att = time_var->get_att("climatology");
      if(climatology_att) {
	// If we have a climatology attribute, we should respect what it says and process the data
	// according to the metadata stored. This code should do that.
	
	// Identify var and fetch climatology bounds data
	char* clim_bnds_name = climatology_att->as_string(0);
	NcVar* clim_bnds_var = out.get_var(clim_bnds_name);
	vector<int> clim_bnds = get_time_values<double>(clim_bnds_var, temp.start_day);
	delete[] clim_bnds_name;
	delete climatology_att;

	// Establish mapping structures
	map<int,int> mon2index;
        vector<vector<int> > clim_months(dims[0]->size());

	// Use the climatology information to determine what months should be included in the sums
	for(int i = 0; i < dims[0]->size(); ++i) {
	  int start_mon = get_month(temp.calendar_type, clim_bnds[i * 2]);
	  int mid_mon = get_month(temp.calendar_type, time_vals[i]);
	  int end_mon = get_month(temp.calendar_type, clim_bnds[i * 2 + 1]);

	  fprintf(stderr, "start_mon: %i, mid_mon: %i, end_mon: %i\n", start_mon, mid_mon, end_mon);

	  if(start_mon == mid_mon && ((mid_mon + 1) % 12) == end_mon) {
	    // Single month
	    clim_months[i].push_back(start_mon);
	    mon2index[start_mon] = i;
	  } else {
	    // Tricky code. We need to handle climatologies without wrap-around, climatologies with wrap-around, 
	    // and annual climatologies. This handles wrap-around with the modulus, annual by testing only after 
	    // changing the value of cur_mon, and normal without any difficulty.
	    int term_month = end_mon % 12;
	    int cur_mon = start_mon; 
	    fprintf(stderr, "months:");
	    do {
	      clim_months[i].push_back(cur_mon);
	      fprintf(stderr, " %i", cur_mon);
	      cur_mon = (cur_mon + 1) % 12;
	    } while(cur_mon != term_month);
	    fprintf(stderr, "\n");
	  }
	}

	// Finally, compute the indices according to the climatology and store them.
	int timestep = 0;
	for(vector<vector<int> >::const_iterator ts_iter = clim_months.begin(); ts_iter != clim_months.end(); ++ts_iter) {
	  std::fill(out_grid, &out_grid[grid_size], 0.0);
	  for(vector<int>::const_iterator mon_iter = (*ts_iter).begin(); mon_iter != (*ts_iter).end(); ++mon_iter) {
	    assert(mon2index.find(*mon_iter) != mon2index.end());
	    int mon_idx = mon2index[*mon_iter];
	    starts[0] = mon_idx;
	    add_derived_var_to_grid(temp_var, prec_var, bits[DV_NAME], *mon_iter, temp_missing, prec_missing, out_missing, out_grid, grid_size, starts, counts);
	  }
	  starts[0] = timestep++;
	  assert(outvar->set_cur(starts));
	  assert(outvar->put(out_grid, counts));
	}

      } else {
	// Just process it month-by-month.
	for(int i = 0; i < dims[0]->size(); ++i) {
	  starts[0] = i;
	  std::fill(out_grid, &out_grid[grid_size], 0.0);
	  
	  // Calculate month
	  int mon = get_month(temp.calendar_type, time_vals[i]);
	  
	  // Call sub.
	  add_derived_var_to_grid(temp_var, prec_var, bits[DV_NAME], mon, temp_missing, prec_missing, out_missing, out_grid, grid_size, starts, counts);
	  
	  assert(outvar->set_cur(starts));
	  assert(outvar->put(out_grid, counts));
	}
      }
    } else {
      // Not going to handle things without an unlimited dim.
      assert(false);
    }

    temp.f->close();
    prec.f->close();
    out.close();

    delete[] edges;
    delete[] out_grid;
  }
}
