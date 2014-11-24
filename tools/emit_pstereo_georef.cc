//#include <iostream>

#include "common.h"

int main(int argc, char** argv) {
  // Try not to fall on your face, netcdf, when a dimension or variable is missing
  ncopts = NC_VERBOSE;

  if(argc < 2) {
    printf("Usage: emit_georef <climatology file> [<bias> <scale factor>]");
    exit(1);
  }

  float scale_factor = 1;
  float bias = 0;

  if(argc > 2) {
    bias = atof(argv[2]);
    if(argc > 3) {
      scale_factor = atof(argv[3]);
    }
  }

  float missing = 1e20f;
  string in_file = argv[1];
  FileRecord f(in_file, false);
  assert(f.is_ok);

  NcVar* var = f.f->get_var(f.var.c_str());

  NcVar* lats = f.f->get_var("lat");
  NcVar* longs = f.f->get_var("lon");
  assert(var && lats && longs);

  NcDim* ydim = f.f->get_dim("yc");
  NcDim* xdim = f.f->get_dim("xc");
  NcDim* toy = var->get_dim(0);

  assert(xdim && ydim && toy);

  NcAtt* missing_value = var->get_att("missing_value");
  if(missing_value) {
    missing = missing_value->as_float(0);
    delete missing_value;
  }

  const int recsize = xdim->size() * ydim->size();
  float* data = new float[recsize * toy->size()];
  float* latdata = new float[recsize];
  float* longdata = new float[recsize];

  long* edges = var->edges();
  int var_recsize = get_recsize_and_edges(var, edges);

  assert(var_recsize == toy->size() * recsize);

  assert(lats->get(latdata, ydim->size(), xdim->size()));
  assert(longs->get(longdata, ydim->size(), xdim->size()));
  assert(var->get(data, edges));

  printf("Lat,Long,Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec,DJF,MAM,JJA,SON,ANN\n");

  for(int i = 0; i < recsize; i++) {
    printf("%f,%f", latdata[i], fmod((longdata[i] + 180.0f), 360.0f) - 180.0f);
    for(int k = 0; k < toy->size(); k++) {
      const int offset = i + k * recsize;
      printf(",%f", (data[offset] * scale_factor) - bias);
    }
    printf("\n");
  }

  delete[] latdata;
  delete[] longdata;
  delete[] data;
}
