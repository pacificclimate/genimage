#ifndef __GENIMAGE_DATAMANAGER_H
#define __GENIMAGE_DATAMANAGER_H

#include <gd.h>
#include <netcdfcpp.h>
#include "scattervars.h"
#include "displayer.h"
#include "config.h"
#include <string>

class DataManager {
 public:
  DataManager(Config& c): config(c) {
    f = 0;
  }
  ~DataManager() {
    if(f) 
      delete f;
  }
    
  DataManager(DataManager& d): config(d.config) {
    model = d.model;
    expt = d.expt;
    timeslice = d.timeslice;
    region = d.region;
    variable = d.variable;
    timeofyear = d.timeofyear;
    f = 0;
  }
  
  void loadScatterVars(ScatterVars* s, bool loadXVar = false) {
    model = s->model;
    expt = s->expt;
    timeslice = s->timeslice;
    if(loadXVar) {
      variable = s->xvariable;
    } else {
      variable = s->yvariable;
    }
  }

  // Calculates a mask, given the selected region and options, of which 
  // data values are inside the region
  bool get_datamask(int* values, const int* slmask, const double* grid_lats, const double* grid_longs);

  // Returns a mask of what squares should be drawn
  bool get_drawmask(int* values, const int* slmask);

  gdImagePtr get_basemap();
  bool get_slmask(int* values);
  int slmask_size() const;
  bool get_longs(double* values);
  int longs_size() const;
  bool get_lats(double* values);
  int lats_size() const;
  bool get_gridlongs(double* values, double* longs = 0);
  int gridlongs_size() const;
  bool get_gridlats(double* values, double* lats = 0);
  int gridlats_size() const;
  bool get_data(double* values);
  int data_size() const;
  bool get_basedata(double* values);

  int numrows() const { return rows->size(); }
  int numcols() const { return cols->size(); }

  void open_datafile();

  std::string model;
  std::string expt;
  std::string timeslice;
  std::string region;
  std::string variable;

  int timeofyear;
  Config& config;

  double minlat, maxlat, minlong, maxlong;
  double difflat, difflong;

private:
  NcFile* f;
  NcDim* rows;
  NcDim* cols;

};

#endif
