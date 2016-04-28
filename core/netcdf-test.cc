#include <netcdfcpp.h>
#include "ConfigFile.h"
#include "genimage.h"

#define TIMESOFYEAR 17

using namespace std;

int main(int argc, char** argv) {
  NcFile* f;
  NcDim* rows;
  NcDim* cols;
  NcVar* slmask;
  NcVar* lats;
  NcVar* longs;
  NcVar* data;

  if (argc < 4) {
    printf("Usage: netcdf-test <outfile.dat> <variable> <timeofyear>\n");
    exit(1);
  }

  f = new NcFile(argv[1]);

  rows = f->get_dim("rows");
  cols = f->get_dim("columns");
  slmask = f->get_var("slmask");
  lats = f->get_var("lats");
  longs = f->get_var("longs");

  data = f->get_var(argv[2]);

  double* values = new double[rows->size() * cols->size()];
  data->set_cur(atoi(argv[3]));
  if (data->get(values, 1, rows->size(), cols->size())) {
    printf("Loaded values...\n");
  }

  for (int i = 0; i < rows->size(); i++) {
    for (int j = 0; j < cols->size(); j++) {
      printf(" % 0.5E", values[(i * cols->size()) + j]);
    }
    printf("\n");
  }

  delete[] values;
  delete f;
}
