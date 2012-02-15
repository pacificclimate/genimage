#include "common.h"
#include "../genimage/derived_vars.h"

#define NUM_COLS 3

enum COLNAMES { DV_NAME=0, TEMP_FILE, PREC_FILE };

/* Expected input format: 3 columns of text, comma separated, on stdin.
    - Output variable name.
    - NetCDF file containing temperature data.
    - NetCDF file containing precipitation data.
 */

string get_filename_create_path(string output_path, FileRecord& temp, string new_var_name) {
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
    ofile += "-" + *p;
    
    // If directory in path doesn't exist, create it.
    if(stat(path.c_str(), &s) != 0) {
      mkdir(path.c_str(), 0777);
    }
  }
  
  ofile = path + "/" + ofile.substr(1) + ".nc";
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

int main(int argc, char** argv) {
  list<FileRecord> l;
  const string temp_var_name = "tas", prec_var_name = "pr";

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

    string ofile = get_filename_create_path(output_path, temp, bits[DV_NAME]);

    NcFile out(ofile.c_str(), NcFile::Replace);
    assert(out.is_valid());
    out.set_fill(NcFile::NoFill);

    copy_all_but_named_var(*(temp.f), out, temp_var_name);

    vector<NcDim*> dims;
    populate_dimvec(temp_var, out, dims);
    NcVar* outvar = out.add_var(bits[DV_NAME].c_str(), temp_var->type(), dims.size(), (const NcDim**)&dims[0]);

    long* edges = outvar->edges();
    long counts[3] = { 1, dims[1]->size(), dims[2]->size() };
    long starts[3] = { 0, 0, 0 };

    cout << "Dim sizes: " << dims[0]->size() << " " << dims[1]->size() << " " << dims[2]->size() << endl;

    int grid_size = dims[2]->size() * dims[1]->size();

    double* prec_grid = new double[grid_size];
    double* temp_grid = new double[grid_size];
    double* out_grid = new double[grid_size];

    if(bits[DV_NAME] == "dl00") {
      // Degree-days less than 0
      outvar->add_att("standard_name", "degree_days_lt_0");
      outvar->add_att("long_name", "Degree-days less than 0 degrees C");
      outvar->add_att("units", "degree-days");
    } else if(bits[DV_NAME] == "dg05") {
      // Degree-days greater than 5
      outvar->add_att("standard_name", "degree_days_gt_5");
      outvar->add_att("long_name", "Degree-days greater than 5 degrees C");
      outvar->add_att("units", "degree-days");
    } else if(bits[DV_NAME] == "dl18") {
      // Degree-days less than 18
      outvar->add_att("standard_name", "degree_days_lt_18");
      outvar->add_att("long_name", "Degree-days less than 18 degrees C");
      outvar->add_att("units", "degree-days");
    } else if(bits[DV_NAME] == "dg18") {
      // Degree-days greater than 18
      outvar->add_att("standard_name", "degree_days_gt_18");
      outvar->add_att("long_name", "Degree-days greater than 18 degrees C");
      outvar->add_att("units", "degree-days");
    } else if(bits[DV_NAME] == "nffd") {
      // Number of frost-free days
      outvar->add_att("standard_name", "number_frost_free_days");
      outvar->add_att("long_name", "Number of frost-free days");
      outvar->add_att("units", "days");
    } else if(bits[DV_NAME] == "pass") {
      // Prec as snow
      outvar->add_att("standard_name", "snowfall_amount");
      outvar->add_att("long_name", "Snowfall amount");
      outvar->add_att("units", "kg m-2");
    }

    cout << "Calendar: " << temp.calendar_type << endl;

    if(dims[0]->is_unlimited()) {
      vector<int> time_vals = get_time_values(*(temp.f), temp.start_day);
      
      // Just process it month-by-month.
      for(int i = 0; i < dims[0]->size(); ++i) {
	starts[0] = i;
	temp_var->set_cur(starts);
	outvar->set_cur(starts);
	assert(temp_var->get(temp_grid, counts));
        add_scalar_to_grid(grid_size, temp_grid, -273.15, temp_missing);
	
	// Calculate month
	int total_months = get_total_months(temp.calendar_type, time_vals[i]);
	int mon = total_months % 12;
	
	std::fill(out_grid, &out_grid[grid_size], 0.0);
	if(bits[DV_NAME] == "dl00") {
	  // Degree-days less than 0
	  degree_days_pcic(temp_grid, out_grid, temp_missing, out_missing, 0, LT, mon, grid_size);
	} else if(bits[DV_NAME] == "dg05") {
	  // Degree-days greater than 5
	  degree_days_pcic(temp_grid, out_grid, temp_missing, out_missing, 5, GT, mon, grid_size);
	} else if(bits[DV_NAME] == "dl18") {
	  // Degree-days less than 18
	  degree_days_pcic(temp_grid, out_grid, temp_missing, out_missing, 18, LT, mon, grid_size);
	} else if(bits[DV_NAME] == "dg18") {
	  // Degree-days greater than 18
	  degree_days_pcic(temp_grid, out_grid, temp_missing, out_missing, 18, GT, mon, grid_size);
	} else if(bits[DV_NAME] == "nffd") {
	  // Number of frost-free days
	  number_frost_free_days(temp_grid, out_grid, temp_missing, out_missing, mon, grid_size);
	} else if(bits[DV_NAME] == "pass") {
	  // Prec as snow
	  prec_var->set_cur(starts);
	  assert(prec_var->get(prec_grid, counts));
	  multiply_grid_by_scalar(grid_size, prec_grid, 86400.0, prec_missing);
	  precip_as_snow(temp_grid, prec_grid, out_grid, temp_missing, prec_missing, out_missing, mon, grid_size);
	}
	assert(outvar->put(out_grid, counts));
      }
    } else {
      // Special handling for seasons and annual
      assert(false);
    }

    temp.f->close();
    prec.f->close();
    out.close();

    delete[] edges;
    delete[] temp_grid;
    delete[] prec_grid;
    delete[] out_grid;
  }
}
