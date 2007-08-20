//#include <iostream>

#include "common.h"

int main(int argc, char** argv) {
  list<Range<int> > ranges;
  bool percent_calcs = false;

  // Try not to fall on your face, netcdf, when a dimension or variable is missing
  ncopts = NC_VERBOSE;

  if(argc < 3) {
    printf("Usage: create_anomaly <normal data> <prediction data (will be modified)> [<percent>]\n");
    exit(1);
  }

  float missing = 1e20f;

  string normal_file = argv[1];
  string prediction_file = argv[2];

  NcFile pred(prediction_file.c_str(), NcFile::Write);
  FileRecord norm(normal_file, false);

  assert(norm.is_ok && pred.is_valid());

  // Get units, change units to "change"
  NcVar* pred_var = pred.get_var(norm.var.c_str());
  NcVar* norm_var = norm.f->get_var(norm.var.c_str());
  
  assert(pred_var && norm_var);

  NcAtt* missing_value = pred_var->get_att("missing_value");
  if(missing_value) {
    missing = missing_value->as_float(0);
    delete missing_value;
  }

  NcAtt* units = pred_var->get_att("units");
  assert(units);
  string unitstr;
  {
    char* bs = units->as_string(0);
    unitstr = bs;
    delete[] bs;
  }
  unitstr += " change";
  units->remove();
  delete units;
  pred_var->add_att("units", unitstr.c_str());

  // For each record, subtract the normals from the climatology and write the data back
  long* pred_edges = pred_var->edges();
  int pred_rec_size = get_recsize_and_edges(pred_var, pred_edges);
  float* pred_data = new float[pred_rec_size];
  assert(pred_var->get(pred_data, pred_edges));
      
  long* norm_edges = norm_var->edges();
  int norm_rec_size = get_recsize_and_edges(norm_var, norm_edges);
  float* norm_data = new float[norm_rec_size];
  assert(norm_var->get(norm_data, norm_edges));

  printf("Norm rec size: %i, Pred rec size: %i\n", norm_rec_size, pred_rec_size);

  assert(norm_rec_size == pred_rec_size);

  for(int i = 0; i < pred_rec_size; i++) {
    if(pred_data[i] != missing && norm_data[i] != missing) {
      pred_data[i] -= norm_data[i];
    } else {
      pred_data[i] = missing;
    }
  }

  assert(pred_var->put(pred_data, pred_edges));

  delete[] pred_edges;
  delete[] norm_edges;
  delete[] pred_data;
  delete[] norm_data;
}
