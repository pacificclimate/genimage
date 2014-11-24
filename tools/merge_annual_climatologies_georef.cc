//#include <iostream>

#include "common.h"

int main(int argc, char** argv) {
  list<Range<int> > ranges;

  // Try not to fall on your face, netcdf, when a dimension or variable is missing
  ncopts = NC_VERBOSE;

  if(argc < 3) {
    printf("Usage: create_climatologies <output file> [<input file> ... ]\n");
    exit(1);
  }

  float missing = 1e20f;
  FILE* outfile = fopen(argv[1], "w");

  assert(outfile);

  int rec_size = -1;
  int old_rec_size = -1;
  int lat_rec_size, long_rec_size;
  vector<float*> dat;
  double *lats, *longs;
  for(int i = 0; i < argc - 2; i++) {
    string infile = argv[i + 2];
    FileRecord clim(infile, false);
    assert(clim.is_ok);

    // Get units, change units to "change"
    NcVar* var = clim.f->get_var(clim.var.c_str());
    assert(var);
    NcAtt* missing_value = var->get_att("missing_value");
    if(missing_value) {
      missing = missing_value->as_float(0);
      delete missing_value;
    }


    // For each record, subtract the normals from the climatology and write the data back
    long* edges = var->edges();
    rec_size = get_recsize_and_edges(var, edges);
    
    if(old_rec_size > 0) {
      assert(rec_size == old_rec_size);
    } else {
      NcVar* vlats = clim.f->get_var("lat");
      NcVar* vlongs = clim.f->get_var("lon");
      assert(vlats && vlongs);

      long* lat_edges = vlats->edges();
      lat_rec_size = get_recsize_and_edges(vlats, lat_edges);
      lats = new double[lat_rec_size];

      long* long_edges = vlongs->edges();
      long_rec_size = get_recsize_and_edges(vlongs, long_edges);
      longs = new double[long_rec_size];

      assert(vlats->get(lats, lat_edges));
      assert(vlongs->get(longs, long_edges));

      assert(lat_rec_size * long_rec_size * 17 == rec_size);

      delete[] long_edges;
      delete[] lat_edges;

    }
    float* data = new float[rec_size];
    assert(var->get(data, edges));
    dat.push_back(data);

    old_rec_size = rec_size;
    delete[] edges;
  }

  for(int i = 0; i < lat_rec_size; i++) {
    for(int j = 0; j < long_rec_size; j++) {
      const int offset = (i * long_rec_size) + j;
      fprintf(outfile, "%f,%f", lats[i], ((longs[j] + 180) % 360) - 180);

      vector<float*>::const_iterator k;
      for(k = dat.begin(); k != dat.end(); ++k) {
	fprintf(outfile, ",%f", (*k)[offset]);
      }
      fprintf(outfile, "\n");
    }
  }

  vector<float*>::const_iterator k;
  for(k = dat.begin(); k != dat.end(); ++k) {
    delete[] *k;
  }

  delete[] lats;
  delete[] longs;
  fclose(outfile);
}
