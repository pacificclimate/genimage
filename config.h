#ifndef __GENIMAGE_CONFIG_H
#define __GENIMAGE_CONFIG_H

#include "point.h"
#include "genimage.h"
#include <list>
#include <string>
#include <vector>

using namespace std;

enum SCATTER_TYPE { ST_POINT, ST_REGION };

class Config {
public:
  Config() {
    scatter_type = ST_POINT;
    percentiles = 0;
    zoom_factor = 1;
    plot_type = 0;
    resolution = 0;
    pct_change_data = 0;
    plot_over_ocean = 0;
    point_present = 0;
    xvariable = yvariable = "";
    gcminfofile = fontfile = data_dir = map_dir = "";
    fringe_size = 0;
    center.x = 0;
    center.y = 0;
    upperleft.x = -180;
    upperleft.y = 90;
    lowerright.x = 180;
    lowerright.y = -90;
    region = "world";
    map_projection = "";
    missing = 1e20;
    base_period = "1961_1990";
    use_anomaly = DEFAULT;
    model = model2 = expt = expt2 = timeslice = timeslice2 = baseline_model = baseline_expt = "";
  }
  // Specified in config file
  std::string gcminfofile;
  std::string fontfile;
  std::string data_dir;
  std::string map_dir;

  // Specified on command line
  std::string outfile;
  int zoom_factor;
  int plot_type;
  int pct_change_data;
  int plot_over_ocean;
  int percentiles;
  int point_present;
  int resolution;
  ANOM_TYPE use_anomaly;
  Point upperleft;
  Point lowerright;
  Point center;
  Point dpoint;
  vector<Point> points;
  int scatter_type;
  double fringe_size;
  double missing;
  double threshold;

  list<int> scenario_set;

  // Variables for axes
  string xvariable, yvariable, region, map_projection;
  string model, model2, expt, expt2, timeslice, timeslice2;
  string baseline_model, baseline_expt;
  string base_period;
  int timeofyear, timeofyear2;
};

#endif
