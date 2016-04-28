#include "genimage.h"
#include "canvas.h"
#include "displayer.h"
#include "interpdatamanager.h"
#include "subset_derived_dm.h"
#include "datamanager.h"
#include "config.h"
#include "point.h"
#include "ConfigFile.h"

#include <assert.h>
#include <map>
#include <algorithm>
#include <boost/tokenizer.hpp>
#include <gdal/gdal_priv.h>
#include <gdal/ogr_spatialref.h>

/* Equatorial Cylindrical Equidistant Projection
    * lat = ((-y + center_y) / height) * 180
    * lon = ((x - center_x) / width) * 360
    * x = ((lon / 360) * width) + center_x
    * y = ((-lat / 180) * height) + center_y
*/

/*
  THOUGHTS:
  - Also need some other bits (in the native coordinate system of the data set)
    - List of X coordinates
    - List of Y coordinates
    - List of X coordinates of edges of grid boxes
    - List of Y coordinates of edges of grid boxes
  - All of these need to be subsetted for the selected region / zoom level

  - Use projection parameters to project coordinates of polygon from lat-lon to native coordinate space
    - After that, can use standard routines

  - None of this deals with projecting GCMs into Albers or whatnot
    - The easy way is to translate every line into a curve and translate the drawing, rather than doing anything fancier; worry about it later
    - Can handle masking by also clipping region mask by region being displayed
    - Fundamentals: Don't worry about it

*/

class RegionDescription {
public:
  Point upperleft, lowerright;
  string projection;
  string proj4_string;
};

// NOTE: CONVERT TO auto_ptr ET AL: http://ootips.org/yonat/4dev/smart-pointers.html
using namespace std;

void parseArgs(Config& c, Displayer& disp, int argc, char** argv, map<string, RegionDescription > regions) {
  const int ARG_TIMESLICE2 =  1000; // JMS: Out of letters.
  const int ARG_TIMEOFYEAR2 = 1001;
  const int ARG_EXPT2 = 1002;
  const int ARG_MODEL2 = 1003;
  const int ARG_FRINGESIZE = 1004;
  const int ARG_CENTERPOINT = 1005;
  const int ARG_WARNING = 1006;
  const int ARG_ANOMALY = 1007;
  const int ARG_BASELINE_MODEL = 1008;
  const int ARG_BASELINE_EXPT = 1009;

  int val;
  bool have_center = false;
  Point temp;
  vector<Point>::const_iterator pt_iter;
  char* optstring = (char *)"a:b:c:e:f:g:h:i:j:k:l:m:n:o:p:r:s:t:u:v:w:x:y:z:";
  static struct option options[] = {
    { "xrange-min", 1, 0, 'a'},
    { "xvariable", 1, 0, 'b'},
    { "colour-map", 1, 0, 'c'},
    { "yvariable", 1, 0, 'd'},
    { "xrange-max", 1, 0, 'e'},
    { "identify-text", 1, 0, 'f'},
    { "ocean-plot", 1, 0, 'g'},
    { "box-threshold", 1, 0, 'i' },
    { "poly-point", 1, 0, 'j'},
    { "y-axis-text", 1, 0, 'k'},
    { "x-axis-text", 1, 0, 'l'},
    { "model", 1, 0, 'm'},
    { "model2", 1, 0, ARG_MODEL2},
    { "region", 1, 0, 'n'},
    { "output-file", 1, 0, 'o'},
    { "plot-type", 1, 0, 'p'},
    { "scenario-set", 1, 0, 'q'},
    { "resolution", 1, 0, 'r'},
    { "expt", 1, 0, 's'},
    { "expt2", 1, 0, ARG_EXPT2},
    { "timeslice", 1, 0, 't' },
    { "timeslice2", 1, 0, ARG_TIMESLICE2 },
    { "yrange-min", 1, 0, 'u'},
    { "yrange-max", 1, 0, 'v'},
    { "timeofyear", 1, 0, 'y'},
    { "timeofyear2", 1, 0, ARG_TIMEOFYEAR2},
    { "zoom-factor", 1, 0, 'z'},
    { "show-grid", 0, &disp.grid, 1},
    { "colour-map-inverse", 0, &disp.colour_map_rev, 1},
    { "dynamic-range", 0, &disp.range_dynamic, 1},
    { "percentiles", 0, &c.percentiles, 1},
    { "percent-change-calculations", 0, &c.pct_change_data, 1},
    { "fringe-size", 1, 0, ARG_FRINGESIZE},
    { "center-point", 1, 0, ARG_CENTERPOINT},
    { "show-warning", 1, 0, ARG_WARNING},
    { "use-anomalies", 1, 0, ARG_ANOMALY},
    { "baseline-model", 1, 0, ARG_BASELINE_MODEL },
    { "baseline-expt", 1, 0, ARG_BASELINE_EXPT },
    { "no-region-vertices", 0, &disp.region_vertices, 0 },
    { 0, 0, 0, 0 }
  };
  while ((val = getopt_long(argc, argv, optstring, options, NULL)) != -1) {
    switch (val) {
    case 0:
      break;
    case 'a':
      // xrange-min
      disp.xrange_min = atof(optarg);
      break;
    case 'b':
      // X variable
      c.xvariable = optarg;
      break;
    case 'c':
      // colour-map
      disp.colour_map = atoi(optarg);
      break;
    case 'd':
      // Y variable
      c.yvariable = optarg;
      break;
    case 'e':
      // xrange-max
      disp.xrange_max = atof(optarg);
      break;
    case 'f':
      // identify-text
      disp.identify_text = optarg;
      break;
    case 'g':
      // ocean-plot
      c.plot_over_ocean = atoi(optarg);
      break;
    case 'i':
      // box-threshold
      c.threshold = atof(optarg);
      break;
    case 'j':
      // poly-point
      temp = parse_data_point(optarg);
      if (temp.valid()) {
        // Check for duplicate points in the list of points; if so, reject them.
        bool duplicate = false;
        for (pt_iter = c.points.begin(); pt_iter != c.points.end() && !duplicate; ++pt_iter) {
          duplicate = (fabs((*pt_iter).x - temp.x) < 1e-12 && fabs((*pt_iter).y - temp.y) < 1e-12);
        }
        if (duplicate)
          fprintf(stderr, "Rejecting duplicate point %s\n", optarg);
        else
          c.points.push_back(temp);
      }
      break;
    case 'k':
      // Y axis text
      disp.yaxis_text = optarg;
      break;
    case 'l':
      // X axis text
      disp.xaxis_text = optarg;
      break;
    case 'm':
      // model
      c.model = optarg;
      break;
    case ARG_MODEL2:
      c.model2 = optarg;
      break;
    case 'n':
      // region
      c.region = optarg;
      break;
    case 'o':
      // output-file
      c.outfile = optarg;
      break;
    case 'p':
      // plot-type
      c.plot_type = atoi(optarg);
      break;
    case 'q':
      // scenario-set
      if (1) {
        list<int> li;
        std::string dat(optarg);
        boost::tokenizer<boost::escaped_list_separator<char> > tok(dat);
        boost::tokenizer<boost::escaped_list_separator<char> >::iterator beg = tok.begin();

        for (; beg != tok.end(); ++beg) {
          li.push_back(atoi((*beg).c_str()));
        }
        c.scenario_set = li;
      }
      break;
    case 'r':
      // resolution
      c.resolution = atoi(optarg);
      break;
    case 's':
      // expt
      c.expt = optarg;
      break;
    case ARG_EXPT2:
      c.expt2 = optarg;
      break;
    case 't':
      c.timeslice = optarg;
      break;
    case ARG_TIMESLICE2:
      c.timeslice2 = optarg;
      break;
    case 'u':
      // yrange-min
      disp.yrange_min = atof(optarg);
      break;
    case 'v':
      // yrange-max
      disp.yrange_max = atof(optarg);
      break;
    case 'w':
      // plot-width
      disp.plot_width = atoi(optarg);
      break;
    case 'y':
      // timeofyear
      c.timeofyear = atoi(optarg);
      break;
    case ARG_TIMEOFYEAR2:
      c.timeofyear2 = atoi(optarg);
      break;
    case 'z':
      // zoom-factor
      c.zoom_factor = atoi(optarg);
      break;
    case ARG_FRINGESIZE:
      c.fringe_size = atof(optarg);
      break;
    case ARG_CENTERPOINT:
      assert(parse_data_point(optarg, c.center));
      have_center = true;
      break;
    case ARG_WARNING:
      disp.warning = optarg;
      break;
    case ARG_BASELINE_MODEL:
      c.baseline_model = optarg;
      break;
    case ARG_BASELINE_EXPT:
      c.baseline_expt = optarg;
      break;
    case ARG_ANOMALY:
      switch (atoi(optarg)) {
      case DEFAULT:
        c.use_anomaly = DEFAULT;
        break;
      case ANOMALY:
        c.use_anomaly = ANOMALY;
        break;
      case ABSOLUTE:
        c.use_anomaly = ABSOLUTE;
        break;
      default:
        fprintf(stderr, "--use-anomaly passed invalid argument\n");
        exit(1);
        break;
      }
      break;
    }
  }

  // Set bounding box for region
  assert(regions.find(c.region) != regions.end());
  RegionDescription d = regions[c.region];
  c.upperleft = d.upperleft;
  c.lowerright = d.lowerright;
  c.map_projection = d.projection;
  if (!have_center) {
    c.center = (c.upperleft + c.lowerright) / 2;
  } else {
    if (c.map_projection != "equidistant_cylindrical") {
      fprintf(stderr, "proj4 string for map: %s\n", d.proj4_string.c_str());
      projPJ pj = pj_init_plus(d.proj4_string.c_str());
      assert(pj);
      projUV uv;
      uv.u = c.center.x * DEG_TO_RAD;
      uv.v = c.center.y * DEG_TO_RAD;
      projUV val = pj_fwd(uv, pj);
      c.center = Point(val.u, val.v);
    }
  }
  fprintf(stderr, "window: %f, %f, %f, %f\n", c.upperleft.y, c.lowerright.y, c.upperleft.x, c.lowerright.x);
}

// WARNING: Must calculate means to get ranges and statistics!
MetaData calcMetaData(DataManager& dm, const DataSpec& s, bool calcRanges = true, bool calcStats = true, bool calcMeans = true) {
  MetaData m;

  vector<WPoint > wpoints;
  vector<WPoint > wpoints_pct;
  vector<WPoint >::const_iterator p_iter;
  double base_data_sum = 0;
  double total_squared_weight = 0;
  double data_sum = 0;
  int i, j;

  //  fprintf(stderr, "Begin allocating memory\n");

  const DataGrid<int> mask = dm.get_datamask(s, dm.config.threshold);
  const DataGrid<double> areagrid = dm.get_areagrid(s);
  const DataGrid<double> data = dm.get_data(s);
  const double* dv = data.values().get();
  const double* agv = areagrid.values().get();
  const int* mv = mask.values().get();
  const double missing = data.missing();
  const double* x_grid = data.x_grid().get();
  const double* y_grid = data.y_grid().get();

  m.proj4_string = data.proj4_string();

  if (calcRanges) {
    m.datarange = Range(data.values().get(), data.grid_size(), Range::RANGE_DATA, data.missing());
    m.longrange = Range(data.x_grid().get(), data.xgrid_size());
    m.latrange = Range(data.y_grid().get(), data.ygrid_size());
  }

  if (calcMeans) {
    wpoints.reserve(data.grid_size());

    // Compute sums and areas
    double marea = 0;
    if (s.percent_change && data.anomaly()) {
      wpoints_pct.reserve(data.grid_size());
      const DataGrid<double> base_data = dm.get_data(DataSpec(s.model, s.expt, dm.config.base_period, s.variable, s.timeofyear, ABSOLUTE, dm.config.pct_change_data));
      const double* bdv = base_data.values().get();
      for (i = data.y_size(); i--; ) {
        const int r_off = (i * data.x_size());
        const double* ag_ptr = &agv[r_off];
        const double* data_ptr = &dv[r_off];
        const double* bdata_ptr = &bdv[r_off];
        const int* mask_ptr = &mv[r_off];
        for (j = data.x_size(); j--; ) {
          if (mask_ptr[j] == 1 && data_ptr[j] != missing && bdata_ptr[j] != missing) {
            const double cdat = data_ptr[j];
            const double bdat = bdata_ptr[j];
            const double area = ag_ptr[j];
            const double val = bdat + bdat * (cdat / 100);
            base_data_sum += area * bdat;
            data_sum += area * val;
            total_squared_weight += area * area;
            marea += area;
            if (calcStats) {
              wpoints_pct.push_back(WPoint(cdat, area));
              wpoints.push_back(WPoint(val, area));
            }
            if (calcRanges)
              // FIXME: This is imperfect (but do we care?)
              m.llrange.addPoint(cdat, (x_grid[j * 2] + x_grid[j * 2 + 1]) / 2, (y_grid[i * 2] + y_grid[i * 2 + 1]) / 2);
          }
        }
      }
    } else {
      for (i = data.y_size(); i--; ) {
        const int r_off = (i * data.x_size());
        const double* ag_ptr = &agv[r_off];
        const double* data_ptr = &dv[r_off];
        const int* mask_ptr = &mv[r_off];
        for (j = data.x_size(); j--; ) {
          if (mask_ptr[j] == 1 && data_ptr[j] != missing) {
            const double area = ag_ptr[j];
            const double val = data_ptr[j];
            data_sum += area * val;
            total_squared_weight += area * area;
            marea += area;
            if (calcStats)
              wpoints.push_back(WPoint(val, area));
            if (calcRanges)
              // FIXME: This is imperfect (but do we care?)
              m.llrange.addPoint(val, (x_grid[j * 2] + x_grid[j * 2 + 1]) / 2, (y_grid[i * 2] + y_grid[i * 2 + 1]) / 2);
          }
        }
      }
    }

    m.area = marea;

    // Calculate weighted means
    m.wmean = data_sum / m.area;

    if (calcStats) {
      // Go through the list of pts and weights getting variance components
      double mvariance = 0;
      const double wmean = m.wmean;
      for (p_iter = wpoints.begin(); p_iter != wpoints.end(); ++p_iter) {
        const WPoint& wp = (*p_iter);
        const double residual = wp.p - wmean;
        mvariance += residual * residual * wp.w;
      }

      // Calculate variance
      // See http://pygsl.sourceforge.net/reference/pygsl/node36.html
      mvariance *= m.area / (squared(m.area) - total_squared_weight);
      m.variance = mvariance;

      // Calculate standard deviation
      m.stddev = sqrt(m.variance);

      // If the data is % change data, switch back to % here
      if (dm.config.pct_change_data)
        wpoints = wpoints_pct;

      m.num_boxes = (int)wpoints.size();

      // Work out the median
      if (data.projection() == "equidistant_cylindrical") {
        std::sort(wpoints.begin(), wpoints.end());
        double desired_weight = m.area / 2;
        WPoint curr_pt(0, 0);
        WPoint last_pt(0, 0);
        double wsum = 0;
        for (p_iter = wpoints.begin(); p_iter != wpoints.end() && wsum < desired_weight; ++p_iter) {
          last_pt = curr_pt;
          curr_pt = *p_iter;
          wsum += curr_pt.w;
        }

        // Interpolate and find best fit median
        double wdiff = curr_pt.w;
        m.wmedian = ((desired_weight - (wsum - wdiff)) / wdiff) * last_pt.p + ((wsum - desired_weight) / wdiff) * curr_pt.p;


        // FIXME: REPLACE WITH nth_element (or simple index)
        desired_weight = (double)wpoints.size() / 2;
        wsum = 0;
        for (p_iter = wpoints.begin(); p_iter != wpoints.end() && wsum <= desired_weight; ++p_iter) {
          last_pt = curr_pt;
          curr_pt = *p_iter;
          wsum++;
        }

        // Finish the median calculation
        if (wsum - desired_weight < 1) {
          m.median = curr_pt.p;
        } else {
          // Average of 2
          m.median = (last_pt.p + curr_pt.p) / 2;
        }
      } else {
        if (wpoints.size() % 2 == 1) {
          vector<WPoint >::iterator elem1 = wpoints.begin() + ((wpoints.end() - wpoints.begin()) / 2);
          vector<WPoint >::iterator elem2 = wpoints.begin() + ((wpoints.end() - wpoints.begin()) / 2) + 1;
          std::nth_element(wpoints.begin(), elem1, wpoints.end());
          std::nth_element(wpoints.begin(), elem2, wpoints.end());
          m.median = m.wmedian = ((*elem1).p + (*elem2).p) / 2;
        } else {
          vector<WPoint >::iterator elem = wpoints.begin() + (wpoints.end() - wpoints.begin()) / 2;
          std::nth_element(wpoints.begin(), elem, wpoints.end());
          m.median = m.wmedian = (*elem).p;
        }
      }
    }

    // Correct values to %
    if (s.percent_change && data.anomaly()) {
      m.base_wmean = base_data_sum / m.area;
      if (m.base_wmean == 0)
        m.wmean = 0;
      else
        m.wmean = ((m.wmean - m.base_wmean) / m.base_wmean) * 100;
    }
  }

  m.gridbox_data = missing;
  m.gridbox_x =  m.gridbox_y = 0;
  m.gridbox_long =  m.gridbox_lat = -999;
  if (dm.config.points.size() == 1) {
    // FIXME: This doesn't account for translated points. Should be updated
    // Find the column it's in
    int i, j;
    if (find_in_grid_var(dm.config.points[0].x, data.x_size(), x_grid, i, true) &&
        find_in_grid_var(dm.config.points[0].y, data.y_size(), y_grid, j, (y_grid[0] < y_grid[1]))) {
      m.gridbox_x = i;
      m.gridbox_y = j;
      m.gridbox_data = data.values().get()[(j * data.x_size()) + i];

      // FIXME: This can be in error
      m.gridbox_long = (x_grid[i * 2] + x_grid[i * 2 + 1]) / 2;
      m.gridbox_lat = (y_grid[j * 2] + y_grid[j * 2 + 1]) / 2;
    }
  }

  return m;
}

void showMap(Displayer& disp, const DataGrid<double>& data, const DataGrid<int>& mask, const DataGrid<int>& drawmask, gdImagePtr basemap, const vector<Point>& proj_points, const Window& display_win, const string map_projection, const string outfile, const Config& c) {
  Range datarange(data.values().get(), data.grid_size(), Range::RANGE_DATA, data.missing());
  Range xrange(display_win.left, display_win.right);
  Range yrange(display_win.bottom, display_win.top);
  Legend* legend = disp.getLegend(datarange);

  fprintf(stderr, "Data: min: %f, max: %f\n", datarange.min(), datarange.max());
  fprintf(stderr, "Legend: min: %f, max: %f\n", legend->range.min(), legend->range.max());

  disp.setOffsets(basemap);
  disp.createCanvas(c.fontfile);
  disp.drawMap(basemap, *legend, drawmask, mask, data, display_win);

  if (map_projection == "equidistant_cylindrical") {
    disp.drawTicks(xrange, yrange);
  } else {
    disp.fillTickAreas();
  }
  disp.drawScale(*legend);
  disp.clearIdentifyArea();
  disp.drawIdentifyText();

  // FIXME: Kind of a hack
  if (c.resolution != 3) {
    disp.drawCreditText();
  }
  disp.drawPolygon(proj_points, xrange, yrange, disp.region_vertices == 1);
  disp.fillMapGaps();

  disp.writePng(outfile);

  delete legend;
}

void handleMap(Displayer& disp, DataManager& dm) {
  const DataSpec s(dm.config.model, dm.config.expt, dm.config.timeslice, dm.config.xvariable, dm.config.timeofyear, dm.config.use_anomaly, dm.config.pct_change_data);
  const DataGrid<double> data = dm.get_data(s);
  const DataGrid<int> mask = dm.get_datamask(s, dm.config.threshold);
  const DataGrid<int> drawmask = dm.get_drawmask(s);
  gdImagePtr basemap = dm.get_basemap();
  vector<Point> proj_points = dm.get_projected_points(s);
  const Window display_win = dm.get_display_window();
  showMap(disp, data, mask, drawmask, basemap, proj_points, display_win, dm.config.map_projection, dm.config.outfile, dm.config);
  gdImageDestroy(basemap);
}

// If user supplied a timeslice2, or a timeofyear2 we can do a diff plot
bool canDoDifferencePlot(DataManager& dm) {
  // if specified and different
  return (!dm.config.timeslice2.empty() && dm.config.timeslice != dm.config.timeslice2) ||
         (dm.config.timeofyear2 != -1 && dm.config.timeofyear != dm.config.timeofyear2) ||
         (dm.config.model != dm.config.model2) ||
         (dm.config.expt != dm.config.expt2);
}

// Won't work quite right for world, but oh well
void handleDiffMap(Displayer& disp, DataManager& dm) {
  if (!canDoDifferencePlot(dm)) {
    fprintf(stderr, "Insufficient variables to do difference plot\n");
    return;
  }

  const DataSpec s1(dm.config.model, dm.config.expt, dm.config.timeslice, dm.config.xvariable, dm.config.timeofyear, dm.config.use_anomaly, dm.config.pct_change_data);
  const DataSpec s2(dm.config.model2, dm.config.expt2, dm.config.timeslice2, dm.config.xvariable, dm.config.timeofyear2, dm.config.use_anomaly, dm.config.pct_change_data);
  const DataGrid<double> data1 = dm.get_data(s1);
  const DataGrid<double> data2 = dm.get_data(s2);
  const DataGrid<int> mask = dm.get_datamask(s1, dm.config.threshold);
  const DataGrid<int> drawmask = dm.get_drawmask(s1);
  gdImagePtr basemap = dm.get_basemap();
  vector<Point> proj_points = dm.get_projected_points(s1);
  const Window display_win = dm.get_display_window();

  // BADNESS
  if (dm.config.model == dm.config.model2) {
    // Difference the data, storing result in data1
    const int datasize = data1.grid_size();
    double* values1 = data1.values().get();
    double* values2 = data2.values().get();
    for (int i = 0; i < datasize; i++)
      values1[i] = values2[i] - values1[i];

    showMap(disp, data1, mask, drawmask, basemap, proj_points, display_win, dm.config.map_projection, dm.config.outfile, dm.config);
  } else {
    Range data1_x_range(data1.x_grid().get(), data1.xgrid_size());
    Range data1_y_range(data1.y_grid().get(), data1.ygrid_size());
    Range data2_x_range(data2.x_grid().get(), data2.xgrid_size());
    Range data2_y_range(data2.y_grid().get(), data2.ygrid_size());
    InterpRect rect;
    rect.dst_boxlong = 2.0;
    rect.dst_boxlat =  2.0;
    rect.src_minlong = max(data1_x_range.min(), data2_x_range.min());
    rect.src_minlat =  max(data1_y_range.min(), data2_y_range.min());
    rect.src_maxlong = min(data1_x_range.max(), data2_x_range.max());
    rect.src_maxlat =  min(data1_y_range.max(), data2_y_range.max());

    const DataGrid<double> data1_interp = interpolate_grid(data1, rect);
    const DataGrid<double> data2_interp = interpolate_grid(data2, rect);
    const DataGrid<int> mask_interp = interpolate_grid(mask, rect);
    const DataGrid<int> drawmask_interp = interpolate_grid(mask, rect);

    // Difference the data, storing result in data1
    const int datasize = data1_interp.grid_size();
    double* values1 = data1_interp.values().get();
    double* values2 = data2_interp.values().get();
    for (int i = 0; i < datasize; i++)
      values1[i] = values2[i] - values1[i];

    showMap(disp, data1_interp, mask_interp, drawmask_interp, basemap, proj_points, display_win, dm.config.map_projection, dm.config.outfile, dm.config);
  }
  gdImageDestroy(basemap);
}

void format_lat(char* output, const char* format, float lat) {
  if (lat == 0) {
    sprintf(output, format, lat, "");
  } else if (lat < 0) {
    sprintf(output, format, -lat, "S");
  } else if (lat > 0) {
    sprintf(output, format, lat, "N");
  }
}

void format_lon(char* output, const char* format, float lon) {
  if (lon == 0) {
    sprintf(output, format, lon, "");
  } else if (lon < 0) {
    sprintf(output, format, -lon, "W");
  } else if (lon > 0) {
    sprintf(output, format, lon, "E");
  }
}

void handleScenariosData(DataManager& dm) {
  // Vars needed for header
  std::string region, month;
  const std::string months[17] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "DJF", "MAM", "JJA", "SON", "ANN" };

  // Open
  FILE* outfd;
  if (dm.config.outfile.length() && dm.config.outfile != "-") {
    outfd = fopen(dm.config.outfile.c_str(), "w");
    if (!outfd) {
      fprintf(stderr, "Couldn't open file %s\n", dm.config.outfile.c_str());
      exit(1);
    }
  } else {
    outfd = stdout;
  }

  if (dm.config.region == "cdn")
    region = "Canada";
  else if (dm.config.region == "nam")
    region = "North America";
  else if (dm.config.region == "world")
    region = "World";
  else
    region = dm.config.region;

  // For all months, seasons, and ann
  for (int timeofyear = 0; timeofyear < 17; ++timeofyear) {
    const DataSpec s(dm.config.model, dm.config.expt, dm.config.timeslice, dm.config.xvariable, timeofyear, dm.config.use_anomaly, dm.config.pct_change_data);
    const DataGrid<double> data = dm.get_data(s);
    const int rows = data.y_size();
    const int cols = data.x_size();
    Range x_range(data.x_grid().get(), data.xgrid_size());
    Range y_range(data.y_grid().get(), data.ygrid_size());

    // Dump header (fortran style)
    const char* fmt = "% 8.4f%s"; // Lat/long fortran style format
    char buf[128];

    fprintf(outfd, " (%s) %s %s % 10d% 10d %s:   ", s.timeslice.c_str(), months[timeofyear].c_str(), s.variable.c_str(), cols, rows, region.c_str());

    format_lat(buf, fmt, y_range.max());
    fprintf(outfd, "%s to ", buf);

    format_lat(buf, fmt, y_range.min());
    fprintf(outfd, "%s ", buf);

    format_lon(buf, fmt, x_range.max());
    fprintf(outfd, "%s to ", buf);

    format_lon(buf, fmt, x_range.min());
    fprintf(outfd, "%s\n", buf);

    for (int j = 0; j < rows; ++j) {
      for (int i = 0; i < cols; ++i) {
        fprintf(outfd, " % .5E", data.values().get()[(j * cols) + i]);
      }
      fprintf(outfd, "\n");
    }
  }

  // cleanup
  fclose(outfd);
}

// How to look up value to be dumped
enum DUMP_TYPE { LATS, LONGS, DATA };

// Helper func for writing lats/longs/slmask, as this is common between the two
template <class T> void dumpVals(const DataGrid<T>& g, const char* format, enum DUMP_TYPE type, const Config& config) {
// Open out file
  FILE* outfd;
  if (config.outfile.length() && config.outfile != "-") {
    outfd = fopen(config.outfile.c_str(), "w");
    if (!outfd) {
      fprintf(stderr, "Couldn't open file %s\n", config.outfile.c_str());
      exit(1);
    }
  } else {
    outfd = stdout;
  }

  // Dump the data
  const int rows = g.y_size();
  const int cols = g.x_size();
  switch (type) {
  case LATS:
  {
    double* vals = get_centers_from_grid(g.y_grid().get(), g.ygrid_size());
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++)
        fprintf(outfd, format, vals[i]);
      fprintf(outfd, "\n");
    }
    delete[] vals;
  }
  break;
  case LONGS:
  {
    double* vals = get_centers_from_grid(g.x_grid().get(), g.xgrid_size());
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++)
        fprintf(outfd, format, vals[j]);
      fprintf(outfd, "\n");
    }
    delete[] vals;
  }
  break;
  case DATA:
  {
    T* vals = g.values().get();
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++)
        fprintf(outfd, format, vals[(i * cols) + j]);
      fprintf(outfd, "\n");
    }
  }
  break;
  }

  // Cleanup
  if (config.outfile.length() && config.outfile != "-")
    fclose(outfd);
}

void handleLatsData(DataManager& dm) {
  const DataSpec s(dm.config.model, dm.config.expt, dm.config.timeslice, dm.config.xvariable, dm.config.timeofyear);
  const DataGrid<double> data = dm.get_data(s);
  dumpVals(data, " % 9.3f", LATS, dm.config);
}

void handleLongsData(DataManager& dm) {
  const DataSpec s(dm.config.model, dm.config.expt, dm.config.timeslice, dm.config.xvariable, dm.config.timeofyear);
  const DataGrid<double> data = dm.get_data(s);
  dumpVals(data, " % 9.3f", LONGS, dm.config);
}

void handleSlmaskData(DataManager& dm) {
  const DataSpec s(dm.config.model, dm.config.expt, dm.config.timeslice, dm.config.xvariable, dm.config.timeofyear);
  const DataGrid<int> slmask = dm.get_slmask(s);
  dumpVals(slmask, " %d", DATA, dm.config);
}

void handleText(Displayer& disp, DataManager& dm) {
  FILE* outfd;
  if (dm.config.outfile.length() && dm.config.outfile != "-") {
    outfd = fopen(dm.config.outfile.c_str(), "w");
    if (!outfd) {
      fprintf(stderr, "Couldn't open file %s\n", dm.config.outfile.c_str());
      exit(1);
    }
  } else {
    outfd = stdout;
  }

  const DataSpec s(dm.config.model, dm.config.expt, dm.config.timeslice, dm.config.xvariable, dm.config.timeofyear, dm.config.use_anomaly, dm.config.pct_change_data);
  MetaData m = calcMetaData(dm, s);
  Legend* legend = disp.getLegend(m.datarange);
  gdImagePtr basemap = dm.get_basemap();
  disp.setOffsets(basemap);

  fprintf(outfd, "Lat range: (%0.2f)-(%0.2f)\n", dm.config.lowerright.y, dm.config.upperleft.y);
  fprintf(outfd, "Lon range: (%0.2f)-(%0.2f)\n", dm.config.upperleft.x, dm.config.lowerright.x);
  fprintf(outfd, "Data range: (%0.6f)-(%0.6f)\n", m.datarange.min(), m.datarange.max());
  fprintf(outfd, "Legend colour range: (%f)-(%f)\n", legend->range.min(), legend->range.max());
  fprintf(outfd, "Size: (%i)-(%i)\n", disp.img_width, disp.img_height);
  fprintf(outfd, "Map size: (%i)-(%i)\n", disp.plot_width, disp.plot_height);
  fprintf(outfd, "Map offset: (%i)-(%i)\n", disp.plot_offset_x, disp.plot_offset_y);
  disp.setScatterOffsets();
  fprintf(outfd, "Scatter size: (%i)-(%i)\n", disp.img_width, disp.img_height);

  fprintf(outfd, "Selection area weighted mean: (%0.6f)\n", m.wmean);
  fprintf(outfd, "Selection area weighted median: (%0.6f)\n", m.wmedian);
  fprintf(outfd, "Selection area weighted standard deviation: (%0.6f)\n", m.stddev);
  fprintf(outfd, "Selection median: (%0.6f)\n", m.median);
  fprintf(outfd, "Selection data min (d, lon, lat): (%0.6f)-(%0.2f)-(%0.2f)\n", m.llrange.getMin(), m.llrange.getMinLong(), m.llrange.getMinLat());
  fprintf(outfd, "Selection data max (d, lon, lat): (%0.6f)-(%0.2f)-(%0.2f)\n", m.llrange.getMax(), m.llrange.getMaxLong(), m.llrange.getMaxLat());
  fprintf(outfd, "Selection area: (%.0f) km<sup>2</sup>\n", m.area);
  fprintf(outfd, "Selection num grid boxes: (%i)\n", m.num_boxes);

  fprintf(outfd, "Data proj4 string: (%s)\n", m.proj4_string.c_str());

  // If the user wants to know about the data box that a specific lat/lon
  // is in, report back
  if (dm.config.points.size() == 1) {
    fprintf(outfd, "Grid box: (%i)-(%i)\n", m.gridbox_x, m.gridbox_y);
    fprintf(outfd, "Data point (d, lon, lat): (%0.6f)-(%0.2f)-(%0.2f)\n", m.gridbox_data, m.gridbox_long, m.gridbox_lat);
  }

  if (dm.config.outfile.length() && dm.config.outfile != "-") {
    fclose(outfd);
  }

  delete legend;
  gdImageDestroy(basemap);
}

void handleMask(Displayer& disp, DataManager& dm) {
  const DataSpec s(dm.config.model, dm.config.expt, dm.config.timeslice, dm.config.xvariable, dm.config.timeofyear);
  const DataGrid<int> datamask = dm.get_datamask(s, dm.config.threshold);

  FILE* outfd;
  if (dm.config.outfile.length() && dm.config.outfile != "-") {
    outfd = fopen(dm.config.outfile.c_str(), "w");
    if (!outfd) {
      fprintf(stderr, "Couldn't open file %s\n", dm.config.outfile.c_str());
      exit(1);
    }
  } else {
    outfd = stdout;
  }

  int rows = datamask.y_size();
  int cols = datamask.x_size();
  const int* mask = datamask.values().get();
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++)
      fprintf(outfd, "% i", mask[(i * cols) + j]);
    fprintf(outfd, "\r\n");
  }

  if (dm.config.outfile.length() && dm.config.outfile != "-") {
    fclose(outfd);
  }
}

void handleGeoref(Displayer& disp, DataManager& dm) {
  const DataSpec s(dm.config.model, dm.config.expt, dm.config.timeslice, dm.config.xvariable, dm.config.timeofyear);
  const DataGrid<int> datamask = dm.get_datamask(s, dm.config.threshold);

  FILE* outfd;
  if (dm.config.outfile.length() && dm.config.outfile != "-") {
    outfd = fopen(dm.config.outfile.c_str(), "w");
    if (!outfd) {
      fprintf(stderr, "Couldn't open file %s\n", dm.config.outfile.c_str());
      exit(1);
    }
  } else {
    outfd = stdout;
  }

  // Snatch up the data
  list<DataGrid<double> > fulldata;
  for (int i = 0; i < NUM_TIMESLICES; i++) {
    const DataSpec s2(dm.config.model, dm.config.expt, dm.config.timeslice, dm.config.xvariable, i, dm.config.use_anomaly, dm.config.pct_change_data);
    fulldata.push_back(dm.get_data(s2));
  }

  int rows = datamask.y_size();
  int cols = datamask.x_size();
  const int* mask = datamask.values().get();
  const double* lats = get_centers_from_grid(datamask.y_grid().get(), datamask.y_size());
  const double* longs = get_centers_from_grid(datamask.x_grid().get(), datamask.x_size());
  fprintf(outfd, "Y,X,Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec,DJF,MAM,JJA,SON,ANN\r\n");
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      const int offset = (i * cols) + j;
      if (mask[offset]) {
        fprintf(outfd, "%9.4f, %9.4f", lats[i], longs[j]);
        for (list<DataGrid<double> >::const_iterator q = fulldata.begin(); q != fulldata.end(); ++q) {
          fprintf(outfd, ",%E", (*q).values().get()[offset]);
        }
        fprintf(outfd, "\r\n");
      }
    }
  }

  if (dm.config.outfile.length() && dm.config.outfile != "-") {
    fclose(outfd);
  }

  delete[] lats;
  delete[] longs;
}

void handlePlotInfo(Displayer& disp, DataManager& dm) {
  const DataSpec s(dm.config.model, dm.config.expt, dm.config.timeslice, dm.config.xvariable, dm.config.timeofyear, dm.config.use_anomaly, dm.config.pct_change_data);
  const DataGrid<double> data = dm.get_data(s);

  FILE* outfd;
  if (dm.config.outfile.length() && dm.config.outfile != "-") {
    outfd = fopen(dm.config.outfile.c_str(), "w");
    if (!outfd) {
      fprintf(stderr, "Couldn't open file %s\n", dm.config.outfile.c_str());
      exit(1);
    }
  } else {
    outfd = stdout;
  }

  gdImagePtr basemap = dm.get_basemap();
  disp.setOffsets(basemap);

  Range xrange(data.x_grid().get(), data.xgrid_size());
  Range yrange(data.y_grid().get(), data.ygrid_size());

  fprintf(outfd, "Lat range: (%0.2f)-(%0.2f)\n", yrange.min(), yrange.max());
  fprintf(outfd, "Lon range: (%0.2f)-(%0.2f)\n", xrange.min(), xrange.max());
  fprintf(outfd, "Size: (%i)-(%i)\n", disp.img_width, disp.img_height);
  fprintf(outfd, "Map size: (%i)-(%i)\n", disp.plot_width, disp.plot_height);
  fprintf(outfd, "Map offset: (%i)-(%i)\n", disp.plot_offset_x, disp.plot_offset_y);

  if (dm.config.outfile.length() && dm.config.outfile != "-") {
    fclose(outfd);
  }

  gdImageDestroy(basemap);
}

void handleRegionOnly(Displayer& disp, DataManager& dm) {
  const DataSpec s(dm.config.model, dm.config.expt, dm.config.timeslice, dm.config.xvariable, dm.config.timeofyear);

  FILE* outfd;
  if (dm.config.outfile.length() && dm.config.outfile != "-") {
    outfd = fopen(dm.config.outfile.c_str(), "w");
    if (!outfd) {
      fprintf(stderr, "Couldn't open file %s\n", dm.config.outfile.c_str());
      exit(1);
    }
  } else {
    outfd = stdout;
  }

  gdImagePtr basemap = dm.get_basemap();
  vector<Point> proj_points = dm.get_projected_points(s);

  Window display_win = dm.get_display_window();
  Range xrange(display_win.left, display_win.right);
  Range yrange(display_win.bottom, display_win.top);

  disp.setOffsets(basemap);
  disp.createCanvas(dm.config.fontfile);
  disp.copyMap(basemap);

  if (dm.config.map_projection == "equidistant_cylindrical") {
    disp.drawTicks(xrange, yrange);
  } else {
    disp.fillTickAreas();
  }

  disp.fillMapGaps();
  disp.fillGapsRegionMap();
  disp.clearIdentifyArea();
  disp.drawCreditText();
  disp.drawPolygon(proj_points, xrange, yrange, disp.region_vertices == 1);

  disp.writePng(dm.config.outfile);

  if (dm.config.outfile.length() && dm.config.outfile != "-") {
    fclose(outfd);
  }

  gdImageDestroy(basemap);
}

inline bool is_available(string token) {
  return !(token == "o" || token == "0");
}

void readGCMInfo(DataManager& dm, list<vector<string> >& list) {
  ifstream f(dm.config.gcminfofile.c_str());
  if (!f.is_open()) {
    fprintf(stderr, "Could not read gcminfo file!\n");
    exit(3);
  }

  fprintf(stderr, "readGCMInfo\n");

  std::string filedata((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());

  boost::tokenizer<boost::escaped_list_separator<char> > tok(filedata, boost::escaped_list_separator<char>("\\", ",\n", "\""));
  boost::tokenizer<boost::escaped_list_separator<char> >::iterator beg = tok.begin();

  int (*lowercase)(int) = tolower;
  while (beg != tok.end()) {
    int count;
    vector<string> bits(GCMINFO_COLS);
    // Read a line of vars
    string model = *beg;
    std::transform(model.begin(), model.end(), model.begin(), lowercase);
    model.erase(0, model.find_first_not_of("\n\r"));
    bits[0] = model;
    ++beg;
    for (count = 1; count < GCMINFO_COLS && beg != tok.end(); ++beg, ++count)
      bits[count] = *beg;
    if (count == GCMINFO_COLS) {
      list.push_back(bits);
    }
  }

  fprintf(stderr, "Done readGCMInfo\n");

}

// FIXME: This is broken (because we can't add another model without another colour)
map<const string, int> getModelColours() {
  map<const string, int> colours;

  colours["climatebc"] = BLACK;
  colours["climatepp"] = BLACK;
  colours["bccr_bcm20"] = LTRED;
  colours["ccsrnies"] = LTORANGE;
  colours["cgcm1"] = LTGRAY;
  colours["cgcm2"] = DKGRAY;
  colours["cccma_cgcm3"] = LTBLUE;
  colours["cccma_cgcm3_t63"] = BLACK;
  colours["cnrm_cm3"] = LTYELLOW;
  colours["csiromk2b"] = LTYELLOWGRN;
  colours["csiro_mk30"] = DKYELLOWGRN;
  colours["echam4"] = LTGREENYLW;
  colours["mpi_echam5"] = DKGREENYLW;
  colours["hadcm2"] = BRTGREEN;
  colours["hadcm3"] = LTGREEN;
  colours["gfdlr15"] = BRTBLUEGRN;
  colours["gfdlr30"] = LTBLUEGRN;
  colours["gfdl_cm20"] = DKBLUEGRN;
  colours["gfdl_cm21"] = DKRBLUEGRN;
  colours["giss_aom"] = LTGCYAN;
  colours["giss_eh"] = DKGCYAN;
  colours["giss_er"] = DKRGCYAN;
  colours["iap_fgoals10g"] = LTBCYAN;
  colours["inmcm30"] = LTCBLUE;
  colours["ipsl_cm4"] = DKRGRAY;
  colours["miroc32_medres"] = LTPURPLE;
  colours["miroc32_hires"] = DKPURPLE;
  colours["miub_echog"] = LTPMAGENTA;
  colours["mri_cgcm232a"] = LTMAGENTA;
  colours["ncarpcm"] = LTRMAGENTA;
  colours["ncar_pcm1"] = DKRMAGENTA;
  colours["ncar_ccsm30"] = DKRRMAGENTA;
  colours["ukmo_hadcm3"] = DKGREEN;
  colours["ukmo_hadgem1"] = DKRGREEN;

  return colours;
}

// FIXME: This is broken (Can't add new experiments without pain)
map<const string, enum SYMBOL> getExptSymbols() {
  map<const string, enum SYMBOL> symbols;

  symbols["HIST"] = CIRCLE;
  symbols["A1"] = LTRIANGLE;
  symbols["A2"] = UTRIANGLE;
  symbols["B1"] = RTRIANGLE;
  symbols["ABS_A1"] = LTRIANGLE;
  symbols["ABS_A2"] = UTRIANGLE;
  symbols["ABS_B1"] = RTRIANGLE;
  symbols["B2"] = DIAMOND;
  symbols["A1FI"] = STAR6;
  symbols["A1T"] = DTRIANGLE;
  symbols["ga"] = CIRCLE;
  symbols["gg"] = SQUARE;

  return symbols;
}

void loadData(DataManager& dm, list<ScatterVars* >& vars, list<LegendToken* >& legend_bits, Range& xrange, Range& yrange, bool stick = false, bool use_symbols = true, bool use_colours = true) {
  string last_model = "";
  string last_expt = "";
  string last_nmodel = "";
  string last_nexpt = "";
  LegendToken* tok = 0;
  LegendToken* special = 0;
  LegendToken* ensemble = 0;
  list<ScatterVars*>::iterator vars_iter;

  map<const string, int> colours = getModelColours();
  map<const string, enum SYMBOL> symbols = getExptSymbols();

  for (vars_iter = vars.begin(); vars_iter != vars.end(); vars_iter++) {
    // Set symbols to use on map and legend
    if ((*vars_iter)->model != last_model || (*vars_iter)->expt != last_expt) {
      int colour;
      bool filled = true;
      enum SYMBOL symbol;

      if (stick) {
        /* On stick plot, colour only the model passed in as --model and --expt */
        if ((*vars_iter)->model == dm.config.model && (*vars_iter)->expt == dm.config.expt) {
          colour = LTBLUE;
        } else {
          colour = TGRAY;  // TODO change this colour so it's less stupid
        }

        symbol = CIRCLE;
      } else {
        if (use_symbols) {
          string id = (*vars_iter)->expt.substr(0, 2);
          if (symbols.find(id) != symbols.end()) {
            symbol = symbols[id];
          } else {
            symbol = NONE;
          }
        } else {
          symbol = NONE;
        }
        if (use_colours) {
          assert(colours.find((*vars_iter)->model) != colours.end());
          colour = colours[(*vars_iter)->model];
        } else {
          colour = LTGRAY;
        }
      }

      if (symbols.find((*vars_iter)->expt) != symbols.end()) {
        // Special case: If this is A1FI or such...
        if (use_symbols) {
          symbol = symbols[(*vars_iter)->expt];
        } else {
          symbol = NONE;
        }
        special = new LegendToken((*vars_iter)->model + " " + (*vars_iter)->expt, colour, symbol, true);
        legend_bits.push_back(special);
      } else if ((*vars_iter)->expt.substr(2, 1) == "x") {
        // Ensemble run
        ensemble = new LegendToken((*vars_iter)->model + " " + (*vars_iter)->expt, colour, symbol, false);
        legend_bits.push_back(ensemble);
      } else {
        // Normal...

        // If the first 2 digits of the experiment match, tack the last digit onto the end of the name
        if (!stick && (*vars_iter)->model == last_nmodel && (*vars_iter)->expt.substr(0, 2) == last_nexpt.substr(0, 2)) {
          if ((*vars_iter)->expt.find("run") == string::npos) {
            tok->name += "," + (*vars_iter)->expt.substr(2, 1);
          } else {
            tok->name += "," + (*vars_iter)->expt.substr((*vars_iter)->expt.find("run") + 3, 1);
          }
        } else {
          // If they aren't the same, it's time to push a new token up
          tok = new LegendToken((*vars_iter)->model + " " + (*vars_iter)->expt, colour, symbol, filled);
          legend_bits.push_back(tok);
        }

        last_nmodel = (*vars_iter)->model;
        last_nexpt = (*vars_iter)->expt;
      }
    }

    if (symbols.find((*vars_iter)->expt) != symbols.end()) {
      (*vars_iter)->symbol = special;
    } else if ((*vars_iter)->expt.substr(2, 1) == "x") {
      (*vars_iter)->symbol = ensemble;
    } else {
      (*vars_iter)->symbol = tok;
    }

    // Set the appropriate data file bits and load data
    if ((*vars_iter)->timeslice == "zero") {
      yrange.add(0);
      (*vars_iter)->setYData(0.);
    } else {
      DataSpec s(*(*vars_iter), dm.config.timeofyear, false, dm.config.use_anomaly, dm.config.pct_change_data);
      if (dm.config.points.size() == 1) {
        // If only 1 point, treat as a grid box
        MetaData my = calcMetaData(dm, s, false, false, false);
        (*vars_iter)->setCoord(my.gridbox_long, my.gridbox_lat);
        yrange.add(my.gridbox_data);
        (*vars_iter)->setYData(my.gridbox_data);

        if ((*vars_iter)->xvariable != "") {
          DataSpec sx(*(*vars_iter), dm.config.timeofyear, true, dm.config.use_anomaly, dm.config.pct_change_data);
          MetaData mx = calcMetaData(dm, sx, false, false, false);
          xrange.add(mx.gridbox_data);
          (*vars_iter)->setXData(mx.gridbox_data);
        }
      } else {
        // If not one point, treat as a polygon
        MetaData my = calcMetaData(dm, s, false, false);
        (*vars_iter)->setCoord(-999.0, -999.0);
        yrange.add(my.wmean);
        (*vars_iter)->setYData(my.wmean);

        // Deal with the X variable, if we have one
        if ((*vars_iter)->xvariable != "") {
          DataSpec sx(*(*vars_iter), dm.config.timeofyear, true, dm.config.use_anomaly, dm.config.pct_change_data);
          MetaData mx = calcMetaData(dm, sx, false, false);
          xrange.add(mx.wmean);
          (*vars_iter)->setXData(mx.wmean);
        }
      }
    }
    if ((*vars_iter)->timeslice != "zero") {
      last_expt = (*vars_iter)->expt;
      last_model = (*vars_iter)->model;
    }
  }
}


vector<MetaData> getScenarioSetMetaData(DataManager& dm, list<ScatterVars* >& vars) {
  list<ScatterVars*>::iterator vars_iter;
  vector<MetaData> metadata_results;

  for (vars_iter = vars.begin(); vars_iter != vars.end(); ++vars_iter) {
    DataSpec s(*(*vars_iter), dm.config.timeofyear, false, dm.config.use_anomaly, dm.config.pct_change_data);
    MetaData my = calcMetaData(dm, s, true, true, true);
    metadata_results.push_back(my);
  }
  return metadata_results;
}

void setRanges(Range& xrange, Range& yrange) {
  if (xrange.min() == xrange.max()) {
    if (xrange.min() < 0)
      xrange.setmin(xrange.min() * 2);
    if (xrange.min() > 0)
      xrange.setmax(xrange.min() * 2);
    if (xrange.min() == 0) {
      xrange.setmax(0.01);
      xrange.setmin(-0.01);
    }
  }
  double xtick_spacing = tick_spacing(xrange, DESIRED_XTICKS);
  xrange.setmin(floor(xrange.min() / xtick_spacing) * xtick_spacing);
  xrange.setmax(ceil(xrange.max() / xtick_spacing) * xtick_spacing);
  if (xrange.max() < 0) xrange.setmax(0);
  if (xrange.min() > 0) xrange.setmin(0);
  if (yrange.min() == yrange.max()) {
    if (yrange.min() < 0)
      yrange.setmin(yrange.min() * 2);
    if (yrange.min() > 0)
      yrange.setmax(yrange.min() * 2);
    if (yrange.min() == 0) {
      yrange.setmax(0.01);
      yrange.setmin(-0.01);
    }
  }
  double ytick_spacing = tick_spacing(yrange, DESIRED_YTICKS);
  yrange.setmin(floor(yrange.min() / ytick_spacing) * ytick_spacing);
  yrange.setmax(ceil(yrange.max() / ytick_spacing) * ytick_spacing);
  if (yrange.max() < 0) yrange.setmax(0);
  if (yrange.min() > 0) yrange.setmin(0);
}

// Silly comparison function because C++ doesn't support anonymous functions
bool compareScatterYVar(const ScatterVars* a, const ScatterVars* b) {
  return a->daty < b->daty;
}

// Calculates weighted percentiles given a list of data
void weighted_pctile(list<ScatterVars*>& vars, double pct, ScatterVars* v) {
  list<ScatterVars*>::const_iterator i;
  assert(vars.size() > 0);
  vector<double> dat;
  vector<double> q;

  q.push_back(pct);

  for (i = vars.begin(); i != vars.end(); ++i)
    dat.push_back((*i)->daty);

  v->daty = quantile(dat, q)[0];
  v->datx = (*(vars.begin()))->datx;
  v->lon = v->lat = -999;
}

class Percentile {
public:
  Percentile(string desc, enum SYMBOL sym, double pctile): desc(desc), sym(sym), pctile(pctile) { };
  string desc;
  enum SYMBOL sym;
  double pctile;
};

void calcPercentiles(list<ScatterVars*>& vars, list<LegendToken* >& leg_tokens, const vector<Percentile > percentiles, const vector<string > years) {
  // Set up the list-map combination
  map<const string, list<ScatterVars*> > vars_by_ts;
  list<ScatterVars*> l[years.size()];
  for (unsigned int i = 0; i < years.size(); i++) {
    vars_by_ts[years[i]] = l[i];
  }

  // Populate the lists in the map
  list<ScatterVars*>::const_iterator i;
  for (i = vars.begin(); i != vars.end(); ++i) {
    vars_by_ts[(*i)->timeslice].push_back(*i);
  }

  // Sort the lists
  for (unsigned int i = 0; i < years.size(); i++) {
    vars_by_ts[years[i]].sort(compareScatterYVar);
  }

  // Create new symbols and give them data
  for (unsigned int i = 0; i < percentiles.size(); i++) {
    LegendToken* tok = new LegendToken(percentiles[i].desc, 0x00000000, percentiles[i].sym, true);
    leg_tokens.push_back(tok);
    for (unsigned int j = 0; j < years.size(); j++) {
      if (vars_by_ts[years[j]].size() > 0) {
        ScatterVars* v = new ScatterVars(percentiles[i].desc, "", years[j], "", "");
        weighted_pctile(vars_by_ts[years[j]], percentiles[i].pctile, v);
        v->symbol = tok;
        vars.push_back(v);
      }
    }
  }
}


// Nearly identical to calcPercentiles, just used for boxplots where the Max, 75th, 50th, 25th, Min are what we want
// Heinous hack. Shoehorn anyone?
void calcBoxplotPercentiles(list<ScatterVars*>& vars, list<LegendToken* >& leg_tokens) {
  const int NUM_PCTILES = 5;
  const int NUM_YEARS = 3;

  const double percentiles[] = { 1.0, 0.75, 0.50, 0.25, 0.0 };
  const string years[] = { "2020", "2050", "2080" };
  const string desc[] = { "Max", "75th percentile", "Median", "25th percentile", "Min" };
  const enum SYMBOL symbols[] = { CIRCLE, UTRIANGLE, DIAMOND, DTRIANGLE, CIRCLE };

  // Set up the list-map combination
  map<const string, list<ScatterVars*> > vars_by_ts;
  list<ScatterVars*> l[NUM_PCTILES];
  for (int i = 0; i < NUM_YEARS; i++) {
    vars_by_ts[years[i]] = l[i];
  }

  // Populate the lists in the map
  list<ScatterVars*>::const_iterator i;
  for (i = vars.begin(); i != vars.end(); ++i) {
    vars_by_ts[(*i)->timeslice].push_back(*i);
  }

  // Sort the lists
  for (int i = 0; i < NUM_YEARS; i++) {
    vars_by_ts[years[i]].sort(compareScatterYVar);
  }

  // Create new symbols and give them data
  for (int i = 0; i < NUM_PCTILES; i++) {
    LegendToken* tok = new LegendToken(desc[i], 0x00000000, symbols[i], true);
    leg_tokens.push_back(tok);
    for (int j = 0; j < NUM_YEARS; j++) {
      if (vars_by_ts[years[j]].size() > 0) {
        ScatterVars* v = new ScatterVars(desc[i], "", years[j], "", "");
        v->symbol = tok;

        if (percentiles[i] == 1.0) {
          v->datx = (vars_by_ts[years[j]].back())->datx;
          v->daty = (vars_by_ts[years[j]].back())->daty;
        } else if (percentiles[i] == 0.0) {
          v->datx = (vars_by_ts[years[j]].front())->datx;
          v->daty = (vars_by_ts[years[j]].front())->daty;
        } else {
          weighted_pctile(vars_by_ts[years[j]], percentiles[i], v);
        }
        vars.push_back(v);
      }
    }
  }
}


enum GCMINFO_OFFSETS { MODEL_OFFSET, MODEL_COUNT, JUNK1, SERIES_OFFSET, MODELNAME_OFFSET, EXPT_OFFSET,
                       S2020_OFFSET, S2050_OFFSET, S2080_OFFSET, BASELINE_OFFSET, STARTYEAR_OFFSET, ENDYEAR_OFFSET
                     };
// JMS: Horrendous hack. ScatterTimeslice now does boxplots - this is because they are almost the same plot type.
// MPN: Not called directly anymore.  Mildly sanitized.
void handleScatterTimesliceEtc(Displayer& disp, DataManager& dm, bool textOnly = false, bool doBoxPlot = false, bool includeHistorical = false, bool showBands = false) {
  int count = 0;
  int var_offset = -1;
  const DataSpec s(dm.config.model, dm.config.expt, dm.config.timeslice, dm.config.xvariable, dm.config.timeofyear);

  // Read in GCMINFO file
  list<vector<string> > gcminfo;
  readGCMInfo(dm, gcminfo);
  list<vector<string> >::iterator beg = gcminfo.begin();

  if (dm.config.yvariable == "") {
    fprintf(stderr, "genimage: error: Y variable must be set for plots\n");
    exit(1);
  }

  // Read header and identify variable
  for (count = 0; count < GCMINFO_COLS; ++count) {
    if ((*beg)[count] == dm.config.yvariable)
      var_offset = count;
  }
  ++beg;

  assert(var_offset >= 0);

  fprintf(stderr, "ANOMALY MODE: %i\n", dm.config.use_anomaly);


  // Create a list of the available expts
  list<int>::const_iterator ssb = dm.config.scenario_set.begin();
  list<ScatterVars* > vars;
  for (int i = 0; beg != gcminfo.end() && ssb != dm.config.scenario_set.end(); ++beg, ++i) {
    string model = (*beg)[MODEL_OFFSET];
    string expt = (*beg)[EXPT_OFFSET];
    bool ensemble = (expt.substr(2, 1) == "x");
    bool var_available = is_available((*beg)[var_offset]);
    bool abs_data = (expt.substr(0, 3) == "ABS");

    if (var_available && !(dm.config.use_anomaly == ANOMALY && abs_data) && *ssb == i && !((dm.config.percentiles || showBands) && ensemble)) {
      if (includeHistorical) {
        ScatterVars* v = new ScatterVars(model, expt, "zero", "", dm.config.yvariable);
        v->setXData(1975);
        vars.push_back(v);
      }
      if (is_available((*beg)[S2020_OFFSET])) {
        ScatterVars* v = new ScatterVars(model, expt, "2020", "", dm.config.yvariable);
        v->setXData(2025);
        vars.push_back(v);
      }
      if (is_available((*beg)[S2050_OFFSET])) {
        ScatterVars* v = new ScatterVars(model, expt, "2050", "", dm.config.yvariable);
        v->setXData(2055);
        vars.push_back(v);
      }
      if (is_available((*beg)[S2080_OFFSET])) {
        ScatterVars* v = new ScatterVars(model, expt, "2080", "", dm.config.yvariable);
        v->setXData(2085);
        vars.push_back(v);
      }
    }
    if (*ssb == i) {
      ++ssb;
    }
  }

  // Load in the data
  Range xrange;
  Range yrange;
  list<LegendToken* > leg_tokens;

  loadData(dm, vars, leg_tokens, xrange, yrange, false, !dm.config.percentiles && !showBands, !dm.config.percentiles || showBands);

  if (disp.range_dynamic && !showBands) {
    setRanges(xrange, yrange);
  } else {
    yrange.setmin(disp.yrange_min);
    yrange.setmax(disp.yrange_max);

  }
  if (includeHistorical) {
    xrange.setmin(1960);
  } else {
    xrange.setmin(2010);
  }
  xrange.setmax(2100);

  if (doBoxPlot) {
    calcBoxplotPercentiles(vars, leg_tokens);
  } else if (dm.config.percentiles || showBands) {
    vector<Percentile > percentiles;
    vector<string > years(3);

    percentiles.push_back(Percentile("90th percentile", UTRIANGLE, 0.9));
    if (showBands) {
      percentiles.push_back(Percentile("75th percentile", UTRIANGLE, 0.75));
      percentiles.push_back(Percentile("Median", DIAMOND, 0.5));
      percentiles.push_back(Percentile("25th percentile", UTRIANGLE, 0.25));
    } else {
      percentiles.push_back(Percentile("Median", DIAMOND, 0.5));
    }
    percentiles.push_back(Percentile("10th percentile", DTRIANGLE, 0.1));

    if (includeHistorical) {
      years.push_back(string("zero"));
    }
    years.push_back(string("2020"));
    years.push_back(string("2050"));
    years.push_back(string("2080"));
    calcPercentiles(vars, leg_tokens, percentiles, years);
  }

  if (textOnly) {
    FILE* out = fopen(dm.config.outfile.c_str(), "w");
    if (!out) {
      fprintf(stderr, "Couldn't open file %s\n", dm.config.outfile.c_str());
      exit(1);
    }

    if (out) { // else complain!
      int j = 0;
      string old = "";
      list<ScatterVars*>::iterator i = vars.begin();
      fprintf(out, "Experiment,Latitude,Longitude,2020,2050,2080,");
      for (; i != vars.end(); i++) {
        ScatterVars* s = *i;
        if (old != s->model + s->expt) {
          // Make sure that the # of columns is right
          for (; j > 0; j--) {
            fprintf(out, ",");
          }

          // Clear everything and reset
          fprintf(out, "\n");
          j = 3;

          fprintf(out, "%s %s,%0.2f,%0.2f,", s->model.c_str(), s->expt.c_str(), s->lon, s->lat);
        }
        fprintf(out, "%0.6f,", s->daty);
        old = s->model + s->expt;
        j--;
      }
      // Make sure that the # of columns is right
      for (; j > 0; j--) {
        fprintf(out, ",");
      }

      // Clear everything and reset
      fprintf(out, "\n");

      fclose(out);
    }
  } else { // !textOnly
    // Set up the plot
    if (showBands) {
      Range yrange2;
      for (list<ScatterVars*>::iterator it = vars.begin(); it != vars.end(); ++it) {
        if ( (*it)->symbol->sym != NONE )
          yrange2.add((*it)->daty);
      }
      setRanges(xrange, yrange2);
      if (includeHistorical) {
        xrange.setmin(1960);
      } else {
        xrange.setmin(2010);
      }
      xrange.setmax(2100);

      disp.setBandsOffsets();
      disp.setTicks(xrange, yrange2);
      disp.createCanvas(dm.config.fontfile);

      // Clear the plot area, draw what we want to...
      disp.clearPlot();

      disp.drawBands(vars, xrange, yrange2, s);

      // Draw tick marks, text labels, and title up X and Y axes
      disp.drawTicks(xrange, yrange2);
      disp.drawAxisTitles();

      // Draw identifying text
      disp.clearIdentifyArea();
      //disp.drawCreditText();
      disp.drawIdentifyText();

      // Fill in gaps in the map
      disp.fillMapGaps();
    } else {
      disp.setScatterOffsets();
      disp.setTicks(xrange, yrange);
      disp.createCanvas(dm.config.fontfile);

      // Clear the plot area, draw what we want to...
      disp.clearPlot();

      disp.drawScatterGrid(xrange, yrange);

      if (doBoxPlot) {
        disp.drawBoxPlot(vars, xrange, yrange);
      } else {
        disp.drawLines(vars, xrange, yrange);
        disp.drawScatter(vars, xrange, yrange);
      }

      // Draw tick marks, text labels, and title up X and Y axes
      disp.drawTicks(xrange, yrange);
      disp.drawAxisTitles();

      // Draw identifying text
      disp.clearIdentifyArea();
      disp.drawCreditText();
      disp.drawIdentifyText();

      // Fill in gaps in the map
      disp.fillMapGaps();

      // Draw the legend
      if (doBoxPlot) {
        disp.drawBoxPlotLegend();
      } else {
        disp.drawLegend(leg_tokens);
      }
    }

    // Write out the file
    disp.writePng(dm.config.outfile);
  }

  // Clean up
  list<ScatterVars* >::iterator vars_iter;
  for (vars_iter = vars.begin(); vars_iter != vars.end(); ++vars_iter) {
    delete *vars_iter;
  }

  list<LegendToken* >::iterator lt_iter;
  for (lt_iter = leg_tokens.begin(); lt_iter != leg_tokens.end(); ++lt_iter) {
    delete *lt_iter;
  }
}

void handleScenarioSetMetadata(DataManager& dm) {
  vector<double > percentiles = { 0.9, 0.75, 0.5, 0.25, 0.1 };
  vector<string > percentile_names = { "90th percentile", "75th percentile", "Median", "25th percentile", "10th percentile" };
  int count = 0, var_offset = -1;

  const DataSpec s(dm.config.model, dm.config.expt, dm.config.timeslice, dm.config.xvariable, dm.config.timeofyear);

  // Read in GCMINFO file
  list<vector<string> > gcminfo;
  readGCMInfo(dm, gcminfo);
  list<vector<string> >::iterator beg = gcminfo.begin();

  if (dm.config.yvariable == "") {
    fprintf(stderr, "genimage: error: Y variable must be set for plots\n");
    exit(1);
  }

  // Read header and identify variable
  for (count = 0; count < GCMINFO_COLS; ++count) {
    if ((*beg)[count] == dm.config.yvariable)
      var_offset = count;
  }
  ++beg;
  assert(var_offset >= 0);

  fprintf(stderr, "ANOMALY MODE: %i\n", dm.config.use_anomaly);

  // Create a list of the available expts
  list<int>::const_iterator ssb = dm.config.scenario_set.begin();
  assert(dm.config.scenario_set.begin() != dm.config.scenario_set.end());
  list<ScatterVars* > vars;
  unsigned int gcm_time_off = ((dm.config.timeslice == "2020") ? S2020_OFFSET : ((dm.config.timeslice == "2050") ? S2050_OFFSET : ((dm.config.timeslice == "2080") ? S2080_OFFSET : 0)));
  for (int i = 0; beg != gcminfo.end() && ssb != dm.config.scenario_set.end(); ++beg, ++i) {
    string model = (*beg)[MODEL_OFFSET];
    string expt = (*beg)[EXPT_OFFSET];
    bool ensemble = (expt.substr(2, 1) == "x");
    bool var_available = is_available((*beg)[var_offset]);
    bool ts_available = is_available((*beg)[gcm_time_off]);
    bool abs_data = (expt.substr(0, 3) == "ABS");

    if (var_available && ts_available && !(dm.config.use_anomaly == ANOMALY && abs_data) && *ssb == i && !ensemble)
      vars.push_back(new ScatterVars(model, expt, dm.config.timeslice, "", dm.config.yvariable));

    if (*ssb == i)
      ++ssb;
  }

  vector<MetaData > md = getScenarioSetMetaData(dm, vars);

  vector<double> spatial_min(md.size()), spatial_max(md.size()), spatial_mean(md.size()), spatial_median(md.size()), spatial_sd(md.size());
  vector<string> expt_name;

  assert(md.size() == vars.size());

  for (size_t i = 0; i < md.size(); ++i) {
    spatial_min[i] = md[i].datarange.min();
    spatial_max[i] = md[i].datarange.max();
    spatial_mean[i] = md[i].wmean;
    spatial_median[i] = md[i].wmedian;
    spatial_sd[i] = md[i].stddev;
  }

  for (list<ScatterVars*>::const_iterator vars_iter = vars.begin(); vars_iter != vars.end(); ++vars_iter) {
    expt_name.push_back((*vars_iter)->model + " " + (*vars_iter)->expt);
    delete *vars_iter;
  }

  // Generate percentiles
  vector<double> spatial_min_quantiles = quantile(spatial_min, percentiles);
  spatial_min.insert(spatial_min.end(), spatial_min_quantiles.begin(), spatial_min_quantiles.end());
  vector<double> spatial_max_quantiles = quantile(spatial_max, percentiles);
  spatial_max.insert(spatial_max.end(), spatial_max_quantiles.begin(), spatial_max_quantiles.end());
  vector<double> spatial_mean_quantiles = quantile(spatial_mean, percentiles);
  spatial_mean.insert(spatial_mean.end(), spatial_mean_quantiles.begin(), spatial_mean_quantiles.end());
  vector<double> spatial_median_quantiles = quantile(spatial_median, percentiles);
  spatial_median.insert(spatial_median.end(), spatial_median_quantiles.begin(), spatial_median_quantiles.end());
  vector<double> spatial_sd_quantiles = quantile(spatial_sd, percentiles);
  spatial_sd.insert(spatial_sd.end(), spatial_sd_quantiles.begin(), spatial_sd_quantiles.end());

  expt_name.insert(expt_name.end(), percentile_names.begin(), percentile_names.end());

  assert(spatial_min.size() == spatial_max.size() && spatial_min.size() == spatial_mean.size() && spatial_min.size() == spatial_median.size() && spatial_min.size() == spatial_sd.size() && spatial_min.size() == expt_name.size());

  FILE* out = fopen(dm.config.outfile.c_str(), "w");
  if (!out) {
    fprintf(stderr, "Couldn't open file %s\n", dm.config.outfile.c_str());
    exit(1);
  }

  fprintf(out, "Experiment,min,max,mean,median,sd\n");
  for (size_t i = 0; i < spatial_min.size(); ++i) {
    fprintf(out, "%s,%0.6f,%0.6f,%0.6f,%0.6f,%0.6f\n", expt_name[i].c_str(), spatial_min[i], spatial_max[i], spatial_mean[i], spatial_median[i], spatial_sd[i]);
  }
  fclose(out);
}

//MPN:  Less of a hack now;  proper wrappers to abstract this crap.
void handleScatterTimeslice(Displayer& disp, DataManager& dm) {
  bool textOnly = false, doBoxPlot = false, includeHistorical = false;

  return handleScatterTimesliceEtc(disp, dm, textOnly, doBoxPlot, includeHistorical);
}

void handleBoxplotTimeslice(Displayer& disp, DataManager& dm) {
  bool textOnly = false, doBoxPlot = true, includeHistorical = false;

  dm.config.percentiles = true;
  return handleScatterTimesliceEtc(disp, dm, textOnly, doBoxPlot, includeHistorical);
}

void handleScatterTimesliceText(Displayer& disp, DataManager& dm) {
  bool textOnly = true, doBoxPlot = false, includeHistorical = false;

  return handleScatterTimesliceEtc(disp, dm, textOnly, doBoxPlot, includeHistorical);
}

void handleBoxplotTimesliceText(Displayer& disp, DataManager& dm) {
  bool textOnly = true, doBoxPlot = true, includeHistorical = false;

  dm.config.percentiles = true;
  return handleScatterTimesliceEtc(disp, dm, textOnly, doBoxPlot, includeHistorical);
}

void handleScatterTimesliceHist(Displayer& disp, DataManager& dm) {  // MPN:  for planners
  bool textOnly = false, doBoxPlot = false, includeHistorical = true;

  return handleScatterTimesliceEtc(disp, dm, textOnly, doBoxPlot, includeHistorical);
}

void handleBandsTimeslice(Displayer& disp, DataManager& dm) {
  bool textOnly = false, doBoxPlot = false, includeHistorical = false, showBands = true;

  dm.config.percentiles = true;
  return handleScatterTimesliceEtc(disp, dm, textOnly, doBoxPlot, includeHistorical, showBands);
}

void handleBandsTimesliceHist(Displayer& disp, DataManager& dm) {
  bool textOnly = false, doBoxPlot = false, includeHistorical = true, showBands = true;

  dm.config.percentiles = true;
  return handleScatterTimesliceEtc(disp, dm, textOnly, doBoxPlot, includeHistorical, showBands);
}

// WIP; chainsawed copy of handleScatterTimesliceEtc
void handleStickplot(Displayer& disp, DataManager& dm) {
  //bool textOnly = false, doBoxPlot = false, includeHistorical = false; /* continuity */

  int count = 0;
  int var_offset = -1;
  bool sset_has_modelexpt = false; /* stick */

  // Read in GCMINFO file
  list<vector<string> > gcminfo;
  readGCMInfo(dm, gcminfo);
  list<vector<string> >::iterator beg = gcminfo.begin();

  // Read header and identify variable
  for (count = 0; count < GCMINFO_COLS; ++count) {
    if ((*beg)[count] == dm.config.yvariable)
      var_offset = count;
  }
  assert(var_offset != -1);// MPN  -- this wasn't getting set; no harm in it being checked...
  ++beg;

  // Create a list of the available expts
  list<int>::const_iterator ssb = dm.config.scenario_set.begin();
  list<ScatterVars* > vars;
  for (int i = 0; beg != gcminfo.end() && ssb != dm.config.scenario_set.end(); ++beg, ++i) {
    string model = (*beg)[MODEL_OFFSET];
    string expt = (*beg)[EXPT_OFFSET];

    // stickplot -- even if this is not available we want it flagged, so we don't add it later despite its unavailability.
    if (*ssb == i && model == dm.config.model && expt == dm.config.expt) {
      sset_has_modelexpt |= true;
      fprintf(stderr, "Added key model FROM SSET to stickplot: %s ; %s\n", model.c_str(), expt.c_str());
    }


    bool ensemble = (expt.substr(2, 1) == "x");
    bool var_available = is_available((*beg)[var_offset]);

    if (var_available && *ssb == i && !ensemble) {
      int offset, xdata;

      // boy oh boy is this ugly, but we want fixed horizontal positioning of plot symbols, so...
      if (dm.config.timeslice == "2020") {
        offset = S2020_OFFSET; xdata = 31337; // 2025;
      } else if (dm.config.timeslice == "2050") {
        offset = S2050_OFFSET; xdata = 31337; // 2055;
      } else if (dm.config.timeslice == "2080") {
        offset = S2080_OFFSET; xdata = 31337; // 2085;
      } else assert(0);

      if (is_available((*beg)[offset])) {
        ScatterVars* v = new ScatterVars(model, expt, dm.config.timeslice /*"2020"*/, "", dm.config.yvariable);
        v->setXData(xdata);
        vars.push_back(v);
        // DEBUG
        //  fprintf(stderr, "Added 1 var.\n");
      }
    }
    if (*ssb == i) {
      ++ssb;
    }
  }

  if (!sset_has_modelexpt) {
    fprintf(stderr, "Forcibly added key model to stickplot: %s ; %s\n", dm.config.model.c_str(), dm.config.expt.c_str());
    ScatterVars* v = new ScatterVars(dm.config.model, dm.config.expt, dm.config.timeslice, "", dm.config.yvariable);
    v->setXData(/*xdata*/ 31337);
    vars.push_back(v);
  } else {
    fprintf(stderr, "Did not forcibly add key model to stickplot.\n");
  }


  // Load in the data
  Range xrange;
  Range yrange;
  list<LegendToken* > leg_tokens;  // this is now pretty much just dummy data

  loadData(dm, vars, leg_tokens, xrange, yrange, /* stick = */ true);

  // FIXME set these to a fixed range -somewhere- ?  of course this has to be done per-var, so shouldn't be in -this- code, ever.
  //  if(disp.range_dynamic) {  FIXME -- have this set properly in modperl so it's not forced here.  for now, boring...
  setRanges(xrange, yrange);
  //} else {
  //   yrange.setmin(disp.yrange_min);
  //  yrange.setmax(disp.yrange_max);
  // }
  fprintf(stderr, "Y axis range is %lf to %lf\n", yrange.min(), yrange.max());

  xrange.setmin(31327);
  xrange.setmax(31347);
  //  xrange.setmin(2010);
  //  xrange.setmax(2100);

  // Set up the plot
  disp.setStickplotOffsets();
  disp.setTicks(xrange, yrange);
  disp.createCanvas(dm.config.fontfile);

  disp.clearCanvas(); // clear the whole thing, not just the plot area, so there aren't gaps between elements

  // Clear the plot area, draw what we want to...
  //  disp.clearPlot();

  //  disp.drawScatterGrid(xrange, yrange);
  disp.drawScatter(vars, xrange, yrange, /* stick = */ true);  // FIXME need to sort vars based on colour or something so the correct one is -on top- i.e. drawn last

  // Draw tick marks, text labels, and title up X and Y axes
  disp.drawYTicks(xrange, yrange);
  disp.drawAxisTitles();  // Should stick around anyways, just don't feed it text..

  // Draw identifying text
  //  disp.clearIdentifyArea();
  //  disp.drawCreditText();
  //  disp.drawIdentifyText();

  // Fill in gaps in the map
  disp.fillMapGaps();

  // Write out the file
  disp.writePng(dm.config.outfile);

  // Clean up
  list<ScatterVars* >::iterator vars_iter;
  for (vars_iter = vars.begin(); vars_iter != vars.end(); ++vars_iter) {
    delete *vars_iter;
  }

  list<LegendToken* >::iterator lt_iter;
  for (lt_iter = leg_tokens.begin(); lt_iter != leg_tokens.end(); ++lt_iter) {
    delete *lt_iter;
  }
}


void handleScatterVariable(Displayer& disp, DataManager& dm, bool textOnly = false) {
  int count = 0;
  int xvar_offset = -1;
  int yvar_offset = -1;

  map<const string, int> ts;

  ts["2020"] = S2020_OFFSET;
  ts["2050"] = S2050_OFFSET;
  ts["2080"] = S2080_OFFSET;

  if (dm.config.timeslice == "" || dm.config.region == "" || dm.config.xvariable == "" || dm.config.yvariable == "") {
    fprintf(stderr, "Must specify timeslice, region, xvariable, and yvariable\n");
    exit(1);
  }

  // Read in GCMINFO file
  list<vector<string> > gcminfo;
  readGCMInfo(dm, gcminfo);
  list<vector<string> >::iterator beg = gcminfo.begin();

  // Read header and identify variables
  for (count = 0; count < GCMINFO_COLS; ++count) {
    if ((*beg)[count] == dm.config.xvariable)
      xvar_offset = count;
    if ((*beg)[count] == dm.config.yvariable)
      yvar_offset = count;
  }
  ++beg;

  fprintf(stderr, "ANOMALY MODE: %i\n", dm.config.use_anomaly);

  // Create a list of the available expts
  list<int>::const_iterator ssb = dm.config.scenario_set.begin();
  list<ScatterVars* > vars;
  for (int i = 0; beg != gcminfo.end() && ssb != dm.config.scenario_set.end(); ++beg, ++i) {
    string model = (*beg)[MODEL_OFFSET];
    string expt = (*beg)[EXPT_OFFSET];
    bool xvar_available = is_available((*beg)[xvar_offset]);
    bool yvar_available = is_available((*beg)[yvar_offset]);
    bool ts_available = dm.config.timeslice == "1961_1990" || (ts.end() != ts.find(dm.config.timeslice) && is_available((*beg)[ts[dm.config.timeslice]]));
    bool abs_data = (expt.substr(0, 3) == "ABS");

    if (*ssb == i) {
      if (xvar_available && yvar_available && ts_available && !(dm.config.use_anomaly == ANOMALY && abs_data)) {
        ScatterVars* v = new ScatterVars(model, expt, dm.config.timeslice, dm.config.xvariable, dm.config.yvariable);
        vars.push_back(v);
      }
      ++ssb;
    }
  }

  // Load in the data
  Range xrange;
  Range yrange;
  list<LegendToken* > leg_tokens;

  loadData(dm, vars, leg_tokens, xrange, yrange);

  if (disp.range_dynamic) {
    setRanges(xrange, yrange);
  } else {
    xrange.setmin(disp.xrange_min);
    xrange.setmax(disp.xrange_max);
    yrange.setmin(disp.yrange_min);
    yrange.setmax(disp.yrange_max);
  }

  if (textOnly) {
    FILE* out = fopen(dm.config.outfile.c_str(), "w");
    if (!out) {
      fprintf(stderr, "Couldn't open file %s\n", dm.config.outfile.c_str());
      exit(1);
    }

    if (out) {
      list<ScatterVars*>::iterator i = vars.begin();
      fprintf(out, "Experiment,Longitude,Latitude,%s,%s,\n", dm.config.xvariable.c_str(), dm.config.yvariable.c_str());
      for (; i != vars.end(); i++) {
        ScatterVars* s = *i;
        fprintf(out, "%s %s,%0.2f,%0.2f,%0.6f,%0.6f,\n", s->model.c_str(), s->expt.c_str(), s->lon, s->lat, s->datx, s->daty);
      }
      fclose(out);
    }
  } else {
    // Set up the plot
    disp.setScatterOffsets();
    disp.setTicks(xrange, yrange);
    disp.createCanvas(dm.config.fontfile);

    // Clear the plot area, draw what we want to...
    disp.clearPlot();
    disp.drawScatterGrid(xrange, yrange);
    disp.drawScatter(vars, xrange, yrange);

    // Draw tick marks, text labels, and title up X and Y axes
    disp.drawTicks(xrange, yrange);
    disp.drawAxisTitles();

    // Draw identifying text
    disp.clearIdentifyArea();
    disp.drawCreditText();
    disp.drawIdentifyText();

    // Draw the legend
    disp.drawLegend(leg_tokens);

    // Fill in gaps in the map
    disp.fillMapGaps();

    // Write out the file
    disp.writePng(dm.config.outfile);
  }

  // Clean up
  list<ScatterVars* >::iterator vars_iter;
  for (vars_iter = vars.begin(); vars_iter != vars.end(); ++vars_iter) {
    delete *vars_iter;
  }

  list<LegendToken* >::iterator lt_iter;
  for (lt_iter = leg_tokens.begin(); lt_iter != leg_tokens.end(); ++lt_iter) {
    delete *lt_iter;
  }
}

// FIXME: This will not handle world stuff wrapped around
// FIXME: This will not handle regions that are 1 by x or x by 1
void handleGeotiff(DataManager& dm) {
  const DataSpec s(dm.config.model, dm.config.expt, dm.config.timeslice, dm.config.xvariable, dm.config.timeofyear, dm.config.use_anomaly, dm.config.pct_change_data);
  const DataGrid<int> datamask = dm.get_datamask(s, dm.config.threshold);
  DataGrid<double> data = dm.get_data(s);
  double* dv = data.values().get();
  int* dmask = datamask.values().get();
  double missing = data.missing();

  // Omit data that isn't in region
  // Can't do this very well; may be able to use a transparency mask or other methods...
  //for(int i = 0; i < datamask.grid_size(); i++)
  //  if(dmask[i] == 0)
  //    dv[i] = missing;

  GDALAllRegister();

  GDALDriver* gdal_driver = GetGDALDriverManager()->GetDriverByName("GTiff");
  assert(gdal_driver);
  GDALDataset* out_ds = gdal_driver->Create(dm.config.outfile.c_str(), data.x_size(), data.y_size(), 1 , GDT_Float32, NULL);
  assert(out_ds);

  // Get min/max and deltas for x/y
  double* x_center = get_centers_from_grid(data.x_grid().get(), data.x_size());
  double* y_center = get_centers_from_grid(data.y_grid().get(), data.y_size());
  Range xrange(x_center, data.x_size());
  Range yrange(y_center, data.y_size());
  double xdiff = xrange.range() / data.x_size();
  double ydiff = yrange.range() / data.y_size();
  delete[] x_center;
  delete[] y_center;

  // Set up projection params
  double adfGeoTransform[6] = {xrange.min(), xdiff, 0, yrange.max(), 0, -ydiff };
  assert(out_ds->SetGeoTransform(adfGeoTransform) == CE_None);
  if (data.proj4_string().length() > 0) {
    assert(out_ds->SetProjection(data.proj4_string().c_str()) == CE_None);
  }
  assert(out_ds->RasterIO(GF_Write, 0, 0, data.x_size(), data.y_size(), dv, data.x_size(), data.y_size(), GDT_Float64, 1, NULL, 0, 0, 0) == CE_None);
  out_ds->GetRasterBand(1)->SetNoDataValue(missing);

  delete out_ds;
}

int main(int argc, char ** argv) {
  Config c;
  Displayer disp;

  // Prevent NetCDF from dying when bad things happen, like variables being missing
  NcError n(NcError::silent_nonfatal);

  string configfile = "/usr/local/etc/genimage.cfg";
  ConfigFile cf(configfile);

  cf.readInto(c.gcminfofile, "gcminfofile");
  cf.readInto(c.fontfile, "fontfile");
  cf.readInto(disp.credit_text, "credit_text");
  cf.readInto(c.data_dir, "data_dir");
  cf.readInto(c.map_dir, "map_dir");

  // Read in data on regions
  map<string, RegionDescription > regions;
  string regionlist;
  assert(cf.readInto(regionlist, "regions"));
  boost::tokenizer<boost::escaped_list_separator<char> > tok(regionlist);
  boost::tokenizer<boost::escaped_list_separator<char> >::iterator beg = tok.begin();
  for (; beg != tok.end(); ++beg) {
    string temp;
    RegionDescription d;
    assert(cf.readInto(temp, "region_" + *beg + "_upperleft"));
    assert(parse_data_point(temp.c_str(), d.upperleft));
    assert(cf.readInto(temp, "region_" + *beg + "_lowerright"));
    assert(parse_data_point(temp.c_str(), d.lowerright));
    assert(cf.readInto(d.projection, "region_" + *beg + "_projection"));
    d.proj4_string = "";
    cf.readInto(d.proj4_string, "region_" + *beg + "_proj4string");
    regions[*beg] = d;
  }

  parseArgs(c, disp, argc, argv, regions);
  DerivedSubsetDataManager dm(c);

  // TODO: Check if all required args have been specified.
  switch (c.plot_type) {
  case TYPE_MAP:
    handleMap(disp, dm);
    break;

  case TYPE_TEXT:
    handleText(disp, dm);
    break;

  case TYPE_MASK:
    handleMask(disp, dm);
    break;

  case TYPE_GEOREF:
    handleGeoref(disp, dm);
    break;

  case TYPE_PLOTINFO:
    handlePlotInfo(disp, dm);
    break;

  case TYPE_REGIONONLY:
    handleRegionOnly(disp, dm);
    break;

  case TYPE_SCATTER_TIMESLICE:
    handleScatterTimeslice(disp, dm);
    break;

  case TYPE_SCATTER_TIMESLICE_HIST:
    handleScatterTimesliceHist(disp, dm);
    break;

  case TYPE_BOXPLOT_TIMESLICE:
    handleBoxplotTimeslice(disp, dm);
    break;

  case TYPE_SCATTER_VARIABLE:
    handleScatterVariable(disp, dm);
    break;

  case TYPE_SCATTER_TIMESLICE_TEXT:
    handleScatterTimesliceText(disp, dm);
    break;

  case TYPE_BOXPLOT_TIMESLICE_TEXT:
    handleBoxplotTimesliceText(disp, dm);
    break;

  case TYPE_SCATTER_VARIABLE_TEXT:
    handleScatterVariable(disp, dm, true);
    break;

  case TYPE_STICKPLOT:
    handleStickplot(disp, dm);
    break;

  case TYPE_MAP_DIFFERENCE:
    handleDiffMap(disp, dm);
    break;

  case TYPE_SCENARIO_DATA:
    handleScenariosData(dm);
    break;

  case TYPE_SLMASK_DATA:
    handleSlmaskData(dm);
    break;

  case TYPE_LATS_DATA:
    handleLatsData(dm);
    break;

  case TYPE_LONGS_DATA:
    handleLongsData(dm);
    break;

  case TYPE_GEOTIFF:
    handleGeotiff(dm);
    break;

  case TYPE_BANDS_TIMESLICE:
    handleBandsTimeslice(disp, dm);
    break;

  case TYPE_BANDS_TIMESLICE_HIST:
    handleBandsTimesliceHist(disp, dm);
    break;

  case TYPE_SCENARIO_SET_METADATA:
    handleScenarioSetMetadata(dm);
    break;

  default:
    fprintf(stderr, "Invalid plot type!\n");
    break;
  }
}
