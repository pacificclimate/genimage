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

  assert(var->num_dims() == 3);

  NcDim* latdim = lats->get_dim(0);
  NcDim* longdim = longs->get_dim(0);
  NcDim* toy = var->get_dim(0);

  assert(latdim && longdim && toy);

  NcAtt* missing_value = var->get_att("missing_value");
  if(missing_value) {
    missing = missing_value->as_float(0);
    delete missing_value;
  }

  const int recsize = longdim->size() * latdim->size();
  float* data = new float[recsize * toy->size()];
  double* latdata = new double[latdim->size()];
  double* longdata = new double[longdim->size()];

  assert(lats->get(latdata, latdim->size()));
  assert(longs->get(longdata, longdim->size()));
  assert(var->get(data, toy->size(), latdim->size(), longdim->size()));

  printf("Lat,Long,Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec,DJF,MAM,JJA,SON,ANN\n");

  for(int i = 0; i < latdim->size(); i++) {
    const int latoff = i * longdim->size();
    for(int j = 0; j < longdim->size(); j++) {
      const int longoff = latoff + j;
      printf("%f,%f", latdata[i], fmod((longdata[j] + 180.0f), 360.0f) - 180.0f);
      for(int k = 0; k < toy->size(); k++) {
	const int offset = longoff + k * recsize;
    	printf(",%f", (data[offset] - bias) * scale_factor);
      }
      printf("\n");
    }
  }

  delete[] latdata;
  delete[] longdata;
  delete[] data;
}
