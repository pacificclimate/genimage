#ifndef __GENIMAGE_SCATTERVARS_H
#define __GENIMAGE_SCATTERVARS_H

#include "range.h"
#include "legendtoken.h"
#include <iostream>
#include <string>

class MetaData {
public:
  // Always in native units (ie, %)
  double base_wmean;
  double wmean;
  double median;
  double wmedian;

  // Always in real units (ie, mm/month or whatever)
  double variance;
  double stddev;

  // In km^2
  double area;

  // Number of boxes in region
  int num_boxes;

  // Always in native units (covers only region selected)
  LatLonRange llrange;

  // Spatial extent of data
  Range longrange;
  Range latrange;

  // Always in native units (covers entire data set)
  Range datarange;

  // Parameters of gridbox selected
  int gridbox_x;
  int gridbox_y;
  double gridbox_data;
  double gridbox_lat;
  double gridbox_long;

  std::string proj4_string;
};

class ScatterVars {
public:
  ScatterVars(std::string model, std::string expt, std::string timeslice, std::string xvariable, std::string yvariable = "") {
    this->model = model;
    this->expt = expt;
    this->timeslice = timeslice;
    this->xvariable = xvariable;
    this->yvariable = yvariable;
  }
  void setCoord(double lon, double lat) {
    this->lon = lon;
    this->lat = lat;
  }
  void setXData(double x) {
    datx = x;
  }
  void setYData(double y) {
    daty = y;
  }
  std::string model;
  std::string expt;
  std::string timeslice;
  std::string xvariable;
  std::string yvariable;
  double lon;
  double lat;
  double datx;
  double daty;
  LegendToken* symbol;
};

#endif
