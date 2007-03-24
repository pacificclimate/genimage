//#include <iostream>
#include "common.h"

#define ZERO_C_IN_K 273.15

int main() {
  char buf[1024];

  // Try not to fall on your face, netcdf, when a dimension or variable is missing
  ncopts = NC_VERBOSE;

  // Second step: Read the files in, one by one; extract the climatologies
  while(fgets(buf, 1024, stdin)) {
    *(strchr(buf, '\n')) = '\0';
    string filename = buf;
    float missing = 1e20f;
    NcFile* f = new NcFile(buf, NcFile::Write);
    boost::char_separator<char> sep("/");
    tokenizer<char_separator<char> > tok(filename, sep);
    vector<string> tokens(tok.begin(), tok.end());

    if(!f->is_valid()) {
      printf("Failed to open file %s\n", buf);
      exit(1);
    }

    NcAtt* missing_value = f->get_att("missing_value");
    if(missing_value) {
      missing = missing_value->as_float(0);
      delete missing_value;
    }

    NcVar* v = f->get_var(tokens[VAR].c_str());
    NcDim* time = f->get_dim("time");
    long* edges = v->edges();
    int recsize = get_recsize_and_edges(v, edges);
    float* data = new float[recsize];

    for(int i = 0; i < time->size(); i++) {
      v->set_cur(i);
      v->get(data, edges);
      add_to_grid(recsize, -ZERO_C_IN_K, data, missing);
      v->put(data, edges);
    }

    delete[] data;
    delete[] edges;
    delete f;
  }
}
