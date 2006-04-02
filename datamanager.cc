#include "datamanager.h"
#include "genimage.h"
#include <iostream>

using namespace std;

bool DataManager::get_drawmask(int* values, const int* slmask) {
  int datasize = data_size();
  switch(config.plot_over_ocean) {
  case 2:
    for(int i = 0; i < datasize; i++) {
      values[i] = !slmask[i];
    }
    break;
  case 1:
    for(int i = 0; i < datasize; i++) {
      values[i] = 1;
    }
    break;
  case 0:
    for(int i = 0; i < datasize; i++) {
      values[i] = slmask[i];
    }
    break;
  default:
    cerr << "Invalid plot_over_ocean value!\n";
    break;
  }
  return true;
}

bool DataManager::get_datamask(int* values, const int* slmask, const double* grid_lats, const double* grid_longs) {
  const int rows = numrows();
  const int cols = numcols();
  const int datasize = data_size();
  double* r_grid = new double[datasize];
  for(int i = 0; i < datasize; i++) {
    r_grid[i] = 0;
  }

  get_drawmask(values, slmask);

  if(config.numpoints > 2) {
    draw_polygon(rows, cols, r_grid, grid_lats, grid_longs, config.points, config.numpoints);

    // Mask off all boxes which are outside the poly in the data (stipple) mask
    // Need to allow most-covered grid box to be selected if no others available
    for(int i = 0; i < rows; i++) {
      for(int j = 0; j < cols; j++) {
	double area = (grid_lats[i] - grid_lats[i + 1]) * (grid_longs[j + 1] - grid_longs[j]);
	
	// If coverage is less than 10%, throw this box out
	values[(i * cols) + j] &= (r_grid[(i * cols) + j] / area > 0.10);
      }
    }
  }

  delete[] r_grid;

  return true;
}

gdImagePtr DataManager::get_basemap() {
  stringstream s;
  s << config.map_dir << model << "_" << config.zoom_factor << "_" << config.resolution << "_" << region << ".png";

  FILE* infile;
  infile = fopen(s.str().c_str(), "rb");
  if(!infile) {
    cerr << "Couldn't read map file!\n";
    exit(1);
  }
  gdImagePtr image = gdImageCreateFromPng(infile);
  if(!image) {
    cerr << "Couldn't load map file!\n";
    exit(1);
  }
  fclose(infile);
  return image;
}

bool DataManager::get_slmask(int* values) {
  if(!f)
    return false;
  
  NcVar* slmask = f->get_var("slmask");
  slmask->get(values, numrows(), numcols());
  return true;
}

int DataManager::slmask_size() const {
  if(!f) 
    return -1;
  return numrows() * numcols();
}

bool DataManager::get_longs(double* values) {
  if(!f)
    return false;
  
  NcVar* longs = f->get_var("longs");
  longs->get(values, numcols());
  return true;
}

int DataManager::longs_size() const {
  if(!f) 
    return -1;
  return numcols();
}

bool DataManager::get_lats(double* values) {
  if(!f)
    return false;
  
  NcVar* lats = f->get_var("lats");
  lats->get(values, numrows());
  return true;
}

int DataManager::lats_size() const {
  if(!f) 
    return -1;
  return numrows();
}

bool DataManager::get_gridlongs(double* values, double* longs) {
  bool delete_longs = false;
  if(!f)
    return false;
  
  int cols = numcols();
  if(!longs) {
    longs = new double[longs_size()];
    get_longs(longs);
    delete_longs = true;
  }

  values[0] = longs[0] - ((longs[1] - longs[0]) / 2);
  for(int i = 1; i < cols; i++) {
    values[i] = longs[i] - ((longs[i] - longs[i - 1]) / 2);
  }
  values[cols] = longs[cols - 1] + ((longs[cols - 1] - longs[cols - 2]) / 2);
  
  maxlong = longs[cols - 1] + ((longs[cols - 1] - longs[cols - 2]) / 2);
  minlong = longs[0] - ((longs[1] - longs[0]) / 2);
  difflong = maxlong - minlong;
  
  if(delete_longs) {
    delete[] longs;
  }
  return true;
}

int DataManager::gridlongs_size() const {
  if(!f) 
    return -1;
  return numcols() + 1;
}

bool DataManager::get_gridlats(double* values, double* lats) {
  bool delete_lats = false;
  if(!f)
    return false;

  int rows = numrows();
  if(!lats) {
    lats = new double[lats_size()];
    get_lats(lats);
    delete_lats = true;
  }

  values[0] = lats[0] - ((lats[1] - lats[0]) / 2);
  for(int i = 1; i < rows; i++) {
    values[i] = lats[i] - ((lats[i] - lats[i - 1]) / 2);
  }
  values[rows] = lats[rows - 1] - ((lats[rows - 2] - lats[rows - 1]) / 2);
  
  maxlat = lats[0] - ((lats[1] - lats[0]) / 2);
  minlat = lats[rows - 1] + ((lats[rows - 1] - lats[rows - 2]) / 2);
  difflat = maxlat - minlat;

  if(delete_lats) {
    delete[] lats;
  }
  return true;
}

int DataManager::gridlats_size() const {
  if(!f) 
    return -1;
  return numrows() + 1;
}

bool DataManager::get_data(double* values) {
  if(!f)
    return false;
  
  string varname = scenario + "_" + timeslice + "_" + variable;
  NcVar* data = f->get_var((char*)varname.c_str());
  data->set_cur(timeofyear);
  data->get(values, 1, numrows(), numcols());
  return true;
}

int DataManager::data_size() const {
  if(!f) 
    return -1;
  return numrows() * numcols();
}

bool DataManager::get_basedata(double* values) {
  if(!f)
    return false;
  
  string varname = scenario + "_1961_1990_" + variable;
  NcVar* data = f->get_var((char*)varname.c_str());
  data->set_cur(timeofyear);
  data->get(values, 1, numrows(), numcols());
  return true;
}

void DataManager::open_datafile() {
  if(f) {
    delete f;
  }
  string filename = config.data_dir + model + "_" + region + ".dat";
  //cout << filename << endl;
  f = new NcFile((char*)filename.c_str());
  if(!f->is_valid()) {
    cout << "File not valid... " << endl;
    exit(3);
  }
  rows = f->get_dim("rows");
  cols = f->get_dim("columns");
}
