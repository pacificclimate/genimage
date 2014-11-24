#ifndef __GENIMAGE_DATAMANAGER_H
#define __GENIMAGE_DATAMANAGER_H

#include <gd.h>
#include <netcdfcpp.h>
#include "scattervars.h"
#include "support.h"
#include "config.h"
#include "genimage.h"
#include <string>
#include <algorithm>
#include <list>
#include <vector>
#include <map>
#include <proj_api.h>
#include <assert.h>
#include <boost/shared_ptr.hpp>

template<typename T>
T clip_to_range(T minval, T maxval, T val) {
  return min(max(minval, val), maxval);
}

double* get_centers_from_grid(const double* grid, int size);

class Window {
 public:
  Window() { top = bottom = left = right = 0; } 
 Window(double top, double left, double bottom, double right): top(top), left(left), bottom(bottom), right(right) {}
  Window* intersect(Window w) {
    return new Window(clip_to_range(w.bottom, w.top, top), clip_to_range(w.left, w.right, left), 
		      clip_to_range(w.bottom, w.top, bottom), clip_to_range(w.left, w.right, right));
  }
  double top;
  double left;
  double bottom;
  double right;
};

class Tile {
 public:
 Tile(Window w, gdImagePtr tile): w(w), tile(tile) {}
  Window w;
  gdImagePtr tile;
};

class DataSpec {
 public:
  DataSpec() { model = expt = timeslice = variable = ""; timeofyear = percent_change = 0; anom = DEFAULT; }
 DataSpec(ScatterVars& s, int timeofyear, bool loadXVar, ANOM_TYPE anom=DEFAULT, int percent_change=0): model(s.model), expt(s.expt), timeslice(s.timeslice), variable(loadXVar ? s.xvariable : s.yvariable), timeofyear(timeofyear), anom(anom), percent_change(percent_change) { }
 DataSpec(const string model, const string expt, const string timeslice, const string variable, const int timeofyear, ANOM_TYPE anom=DEFAULT, int percent_change=0): model(model), expt(expt), timeslice(timeslice), variable(variable), timeofyear(timeofyear), anom(anom), percent_change(percent_change) { }
  string model, expt, timeslice, variable;
  int timeofyear;
  ANOM_TYPE anom;
  int percent_change;
};

// This needs to be implemented outside DataManager
//virtual bool get_difference_data(DataObject<double&> diffdata, DataManager& dm2);


template <typename T> 
class DataGrid {
 public:
  DataGrid() { assert(false); }
  DataGrid(const DataSpec& r): r(r) { _x_size = _y_size = -1; _missing = 0; _projection = _proj4_string = ""; }
  void set_data(int x_size, int y_size, T missing, T* values, double* x_grid, double* y_grid, string proj4_string, string projection, bool anomaly, string base_period) { _values.reset(values); _x_grid.reset(x_grid); _y_grid.reset(y_grid); _x_size = x_size; _y_size = y_size; _missing = missing; _proj4_string = proj4_string; _projection = projection; _anomaly = anomaly; _base_period = base_period; }
  void set_data(int x_size, int y_size, T missing, boost::shared_ptr<T> values, boost::shared_ptr<double> x_grid, boost::shared_ptr<double> y_grid, string proj4_string, string projection, bool anomaly, string base_period) { _values = values; _x_grid = x_grid; _y_grid = y_grid; _x_size = x_size; _y_size = y_size; _missing = missing; _proj4_string = proj4_string; _projection = projection; _anomaly = anomaly; _base_period = base_period; }

  boost::shared_ptr<T> values() const { return(_values); }
  boost::shared_ptr<double> x_grid() const { return(_x_grid); }
  boost::shared_ptr<double> y_grid() const { return(_y_grid); }
  int grid_size() const { return(_x_size * _y_size); }
  int xgrid_size() const { return(_x_size * 2); }
  int ygrid_size() const { return(_y_size * 2); }
  int x_size() const { return(_x_size); }
  int y_size() const { return(_y_size); }
  string projection() const { return(_projection); }
  string proj4_string() const { return(_proj4_string); }
  string base_period() const { return(_base_period); }
  bool anomaly() const { return(_anomaly); }
  T missing() const { return(_missing); }
  bool has_data() const { return(_values != 0); }

  // HATE: This should be const, but C++ is a mofo
  DataSpec r;

 private:
  boost::shared_ptr<T> _values;
  boost::shared_ptr<double> _x_grid;
  boost::shared_ptr<double> _y_grid;
  int _x_size, _y_size;
  T _missing;
  string _proj4_string;
  string _projection;
  bool _anomaly;
  string _base_period;
};

enum SLOT{DATA_SLOT=0,BASELINE_SLOT,SLMASK_SLOT,LATS_SLOT,LONGS_SLOT,LAT_BNDS_SLOT,LON_BNDS_SLOT,XC_SLOT,YC_SLOT,PROJ_MAPPING_SLOT,MAX_SLOTS};

class VarsAndFiles {
 public:
  VarsAndFiles() { }
 VarsAndFiles(NcFile *f): f(f) {}
  shared_ptr<NcFile> f;
  map<string, NcVar*> var_cache;
};

class DataProvider {
 public:
 DataProvider(Window& w, Config& c): w(w), c(c) { current_filename = ""; used_for_bilinear = false; }
 DataProvider(Config& c): c(c) { current_filename = ""; used_for_bilinear = false; }
  ~DataProvider() {
  }

  void set_bilin_flag() { used_for_bilinear = true; }
  void set_window(Window w) { this->w = w; }
  int x_size(const DataSpec& s) { open_file(s); NcDim* cols_dim = f->get_dim("columns"); assert(cols_dim); return(cols_dim->size()); }
  int y_size(const DataSpec& s) { open_file(s); NcDim* rows_dim = f->get_dim("rows"); assert(rows_dim); return(rows_dim->size()); }
  int sub_x_size(const DataSpec& s) { return sub_x_size(s, l); }
  int sub_y_size(const DataSpec& s) { return sub_y_size(s, l); }
  int sub_x_size(const DataSpec& s, const list<Window>& l);
  int sub_y_size(const DataSpec& s, const list<Window>& l);
  int xgrid_size(const DataSpec& s) { return(x_size(s) * 2); }
  int ygrid_size(const DataSpec& s) { return(y_size(s) * 2); }
  bool get_windowed_data(const DataSpec& s, SLOT slot, double* values);
  bool get_windowed_mask_data(const DataSpec& s, SLOT slot, int* values);

  template <typename T> T get_missing(const DataSpec& s, const SLOT slot) {
    open_file(s);
    NcVar* v = get_ncdf_var(s, slot);
    assert(v);

    // Get missing value
    T missing;
    NcAtt* missing_att = v->get_att("missing_value");
    if(missing_att) {
      switch(missing_att->type()) {
      case NC_INT:
	missing = missing_att->as_int(0);
	break;
      case NC_FLOAT:
        missing = missing_att->as_float(0);
	break;
      case NC_DOUBLE:
        missing = missing_att->as_double(0);
	break;
      default:
	assert(false);
	break;
      }
      delete missing_att;
    } else {
      missing = 200000000;
    }
    return missing;
  }
  bool is_anomaly(const DataSpec& s, const SLOT slot) {
    NcVar* v = get_ncdf_var(s, slot);
    bool anomaly;
    NcAtt* units_att = v->get_att("units");
    if(units_att) {
      char* units_char = units_att->as_string(0);
      string units(units_char);
      delete[] units_char;
      anomaly = (units.find("change") != string::npos);
      delete units_att;
    } else {
      anomaly = !(s.expt.substr(0, 3) == "ABS") && (s.timeslice == "2080" || s.timeslice == "2050" || s.timeslice == "2020");
    }
    return anomaly;
  }

  string get_projection(const DataSpec& s) {
    NcVar* v = get_ncdf_var(s, SLMASK_SLOT);
    assert(v);

    // Get proj4 string
    string projection = "equidistant_cylindrical";

    NcAtt* grid_mapping = v->get_att("grid_mapping");
    if(grid_mapping) {
      // Get grid mapping information
      char* gm_desc = grid_mapping->as_string(0);
      projection = gm_desc;
      delete[] gm_desc;
      delete grid_mapping;
    }
    return projection;
  }
  string get_proj4_string(const DataSpec& s) {
    string projection = get_projection(s);

    if(projection == "equidistant_cylindrical") {
      return "";
    }

    NcVar* v = get_ncdf_var(s, SLMASK_SLOT);
    assert(v);
    NcVar* mapping = get_ncdf_var(s, PROJ_MAPPING_SLOT);
    assert(mapping);
    NcAtt* p4s = mapping->get_att("proj4_string");
    assert(p4s);
    char* p4s_str = p4s->as_string(0);
    
    string proj4_string = p4s_str;
    
    delete[] p4s_str;
    delete p4s;
    return proj4_string;
  }

  double* get_x_grid(const DataSpec& s) { return get_x_grid(s, l); }
  double* get_y_grid(const DataSpec& s) { return get_y_grid(s, l); }
  double* get_x_grid(const DataSpec& s, const list<Window>& l);
  double* get_y_grid(const DataSpec& s, const list<Window>& l);

 private:
  string get_filename(const DataSpec& s) { string ret; ret.reserve(100); return(ret.append(c.data_dir).append(s.model).append(".dat")); }
  bool get_slice(double* values, NcVar* data, const int timeofyear, const Window& w);
  NcVar* get_ncdf_var(const DataSpec& s, SLOT slot=DATA_SLOT);
  bool open_file(const DataSpec& s) {
    string filename = get_filename(s);
    if(file_cache.find(filename) == file_cache.end()) {
      fprintf(stderr, "Opening file %s\n", filename.c_str());
      NcFile *f2 = new NcFile(filename.c_str());
      assert(f2); assert(f2->is_valid());
      file_cache[filename] = VarsAndFiles(f2);
    }

    if(filename != current_filename) {
      fprintf(stderr, "Changing from file '%s' to file '%s'\n", current_filename.c_str(), filename.c_str());
      current_filename = filename;
      f = file_cache[current_filename].f;
      l = get_subsets(s);
    }

    return true;
  }
  Window w;
  Config& c;
  map<string, VarsAndFiles> file_cache;
  shared_ptr<NcFile> f;
  string current_filename;
  list<Window> l;
  bool is_same_grid(const DataSpec& s) { return(get_filename(s) == current_filename); }
  bool used_for_bilinear;

  // Returns a list of windows into the data, shifted to eliminate the need for stitching/shifting later
  list<Window> get_subsets(const DataSpec& s);
};

class DataManager {
 public:
 DataManager(Config& c): config(c), dp(DataProvider(c)), bp(DataProvider(c)) {
    w = get_display_window();
    dp.set_window(w);
    bp.set_window(w);
    // FIXME: Add a check here for config.baseline_model without config.baseline_expt or vice versa.
  }
  virtual ~DataManager() {
  }
  
  // Returns the window to display data in
  Window get_display_window();

  // Gets the polygon's vertexes in the projection of the specified data grid
  vector<Point> get_projected_points(const DataSpec& o);

  // Calculates a mask, given the selected region and options, of which 
  // data values are inside the region
  virtual DataGrid<int> get_datamask(const DataSpec& datamask, const double threshold);

  // Calculates a grid of areas
  virtual DataGrid<double> get_areagrid(const DataSpec& area);

  // Returns a mask of what squares should be drawn
  virtual DataGrid<int> get_drawmask(const DataSpec& drawmask);

  // Returns the requested sea-land mask
  virtual DataGrid<int> get_slmask(const DataSpec& slmask);

  // Returns the requested data
  virtual DataGrid<double> get_data(const DataSpec& data);

  // Gets basemaps and/or bits of basemaps
  virtual gdImagePtr get_basemap();


  Config& config;
  
 protected:
  virtual gdImagePtr get_basemap_image(Window& w);
  bool use_alt_baseline() { return (config.baseline_model != "" && config.baseline_expt != ""); }

  template<typename T>
    bool set_up_grid(DataGrid<T>& s, SLOT slot=DATA_SLOT, bool use_alt_baseline = false) {
    int xsize, ysize;
    double *x_grid, *y_grid;
    string proj4_string;
    string projection;
    if(use_alt_baseline) {
      dp.set_bilin_flag();
      DataSpec s2(config.baseline_model, config.baseline_expt, s.r.timeslice, s.r.variable, s.r.timeofyear, s.r.anom, s.r.percent_change);
      xsize = bp.sub_x_size(s2);
      ysize = bp.sub_y_size(s2);
      x_grid = bp.get_x_grid(s2);
      y_grid = bp.get_y_grid(s2);
      proj4_string = bp.get_proj4_string(s2);
      projection = bp.get_projection(s2);
    } else {
      xsize = dp.sub_x_size(s.r);
      ysize = dp.sub_y_size(s.r);
      x_grid = dp.get_x_grid(s.r);
      y_grid = dp.get_y_grid(s.r);
      proj4_string = dp.get_proj4_string(s.r);
      projection = dp.get_projection(s.r);
    }
    T* values = new T[xsize * ysize];

    // Get missing value
    T missing = dp.get_missing<T>(s.r, slot);
    bool anomaly = dp.is_anomaly(s.r, slot);

    s.set_data(xsize, ysize, missing, values, x_grid, y_grid, proj4_string, projection, anomaly, config.base_period);
    return(true);
  }

  Window w;
  DataProvider dp;
  DataProvider bp;
};

#endif
