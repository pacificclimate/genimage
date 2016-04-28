#include "genimage.h"
#include "line.h"
#include <iostream>
#include <algorithm>
#include "displayer.h"
#include "interpdatamanager.h"

using namespace std;

enum OCEAN_DRAW_TYPE {DRAW_LAND = 0, DRAW_ALL, DRAW_OCEAN};

int DataProvider::sub_y_size(const DataSpec& s, const list<Window>& l) {
  open_file(s);
  const Window& w = *(l.begin());
  return (int)(w.bottom - w.top) + 1;
}

int DataProvider::sub_x_size(const DataSpec& s, const list<Window>& l) {
  int sum = 0;
  open_file(s);
  for (list<Window>::const_iterator i = l.begin(); i != l.end(); i++)
    sum += (int)((*i).right - (*i).left) + 1;
  return sum;
}

int get_recsize_and_edges(NcVar* invar, long* edges) {
  // Construct the list of what to copy
  int recsize = 1;
  for (int i = 0; i < invar->num_dims(); i++) {
    recsize *= edges[i];
  }
  return recsize;
}

double* get_centers_from_grid(const double* grid, int size) {
  double* ret = new double[size];
  for (int i = 0; i < size; i++)
    ret[i] = (grid[i * 2] + grid[i * 2 + 1]) / 2;
  return ret;
}

NcVar* DataProvider::get_ncdf_var(const DataSpec& s, SLOT slot) {
  open_file(s);
  string identifier;
  identifier.reserve(100);
  open_file(s);

  switch (slot) {
  case DATA_SLOT:
  case BASELINE_SLOT:
    identifier.append(s.expt).append("_").append(s.timeslice).append("_").append(s.variable);
    break;
  case SLMASK_SLOT:
    identifier = "slmask";
    break;
  case LAT_BNDS_SLOT:
    identifier = "lat_bnds";
    break;
  case LON_BNDS_SLOT:
    identifier = "lon_bnds";
    break;
  case LATS_SLOT:
    identifier = "lats";
    break;
  case LONGS_SLOT:
    identifier = "longs";
    break;
  case XC_SLOT:
    identifier = "xc";
    break;
  case YC_SLOT:
    identifier = "yc";
    break;
  case PROJ_MAPPING_SLOT:
    identifier = get_projection(s);
    break;
  default:
    assert(false);
  }

  //fprintf(stderr, "get_ncdf_var: Slot %i, requesting %s\n", slot, identifier.c_str());

  NcVar* cur = 0;
  map<string, NcVar*>::iterator item = file_cache[current_filename].var_cache.find(identifier);
  if (item != file_cache[current_filename].var_cache.end()) {
    cur = item->second;
  } else {
    cur = f->get_var(identifier.c_str());
    file_cache[current_filename].var_cache[identifier] = cur;
  }
  return (cur);
}

Window DataManager::get_display_window() {
  double size_x_2 = (config.lowerright.x - config.upperleft.x) / config.zoom_factor / 2;
  double size_y_2 = (config.upperleft.y - config.lowerright.y) / config.zoom_factor / 2;
  if (config.region == "world") {
    Point ctr(config.center.x,
              clip_to_range(config.lowerright.y + size_y_2, config.upperleft.y - size_y_2, config.center.y));
    //fprintf(stderr, "get_display_window: half-sizes: %f, %f\n", size_y_2, size_x_2);
    return Window(ctr.y + size_y_2, ctr.x - size_x_2, ctr.y - size_y_2, ctr.x + size_x_2);
  } else {
    Point ctr(clip_to_range(config.upperleft.x + size_x_2, config.lowerright.x - size_x_2, config.center.x),
              clip_to_range(config.lowerright.y + size_y_2, config.upperleft.y - size_y_2, config.center.y));
    //fprintf(stderr, "get_display_window: half-sizes: %f, %f\n", size_y_2, size_x_2);
    return Window(ctr.y + size_y_2, ctr.x - size_x_2, ctr.y - size_y_2, ctr.x + size_x_2);
  }
}

DataGrid<double> DataManager::get_areagrid(const DataSpec& s) {
  fprintf(stderr, "dm::get_areagrid\n");
  DataGrid<double> g(s);
  set_up_grid(g, SLMASK_SLOT, use_alt_baseline());
  double* values = g.values().get();
  const double* x_grid = g.x_grid().get();
  const double* y_grid = g.y_grid().get();
  assert(g.grid_size() > 0);
  if (g.projection() == "equidistant_cylindrical") {
    const double er_squared = squared(EARTH_RADIUS);
    const double deg2rad = M_PI / 180;
    const int num_x = g.x_size();
    const int num_y = g.y_size();
    for (int y = 0; y < num_y; y++) {
      double* aptr = &values[y * num_x];
      const double r_weight = er_squared * deg2rad * fabs(sin(y_grid[y * 2] * deg2rad) - sin(y_grid[y * 2 + 1] * deg2rad));
      for (int x = 0; x < num_x; x++) {
        aptr[x] = r_weight * fabs(x_grid[x * 2] - x_grid[x * 2 + 1]);
      }
    }
  } else {
    // Return area in square km not square m
    double box_area = fabs((y_grid[0] - y_grid[1]) * (x_grid[0] - x_grid[1])) / 1000000;
    std::fill(values, &values[g.grid_size()], box_area);
  }
  return g;
}

DataGrid<int> DataManager::get_drawmask(const DataSpec& s) {
  fprintf(stderr, "dm::get_drawmask\n");
  DataGrid<int> g = DataManager::get_slmask(s);
  int* values = g.values().get();
  const int datasize = g.grid_size();

  switch (config.plot_over_ocean) {
  case DRAW_OCEAN:
    for (int i = 0; i < datasize; i++)
      values[i] = values[i] ^ 1;
    break;
  case DRAW_ALL:
    std::fill(values, &values[datasize], 1);
    break;
  case DRAW_LAND:
    break;
  default:
    fprintf(stderr, "Invalid plot_over_ocean value!\n");
    assert(false);
    break;
  }
  return g;
}

vector<Point> DataManager::get_projected_points(const DataSpec& s) {
  int numpoints = config.points.size();
  vector<Point> points(numpoints);
  if (numpoints > 0) {
    // FIXME: THIS IS BROKEN. NEEDS TO USE BASELINE.
    if (dp.get_projection(s) == "equidistant_cylindrical") {
      for (int i = 0; i < numpoints; i++)
        points[i] = config.points[i];
    } else {
      projPJ pj = pj_init_plus(dp.get_proj4_string(s).c_str());
      assert(pj);
      for (int i = 0; i < numpoints; i++) {
        projUV uv;
        uv.u = config.points[i].x * DEG_TO_RAD;
        uv.v = config.points[i].y * DEG_TO_RAD;
        projUV val = pj_fwd(uv, pj);
        points[i] = Point(val.u, val.v);
        fprintf(stderr, "Point %i: x -> %f, y: %f\n", i, points[i].x, points[i].y);
      }
    }
  }
  return points;
}

DataGrid<int> DataManager::get_datamask(const DataSpec& s, const double threshold) {
  fprintf(stderr, "dm::get_datamask\n");

  // Grab in all the goodies I want to use
  DataGrid<int> g = DataManager::get_drawmask(s);
  int* values = g.values().get();
  const double* x_grid = g.x_grid().get();
  const double* y_grid = g.y_grid().get();

  const int cols = g.x_size();
  const int rows = g.y_size();
  const int datasize = g.grid_size();

  // Allocate the fractional coverage grid, fill with 0
  double* r_grid = new double[datasize];
  std::fill(r_grid, &r_grid[datasize], 0.0);

  // Get points, get ready to use them
  vector<Point> points = get_projected_points(s);
  const int numpoints = points.size();

  if (numpoints == 1) {
    // Case 1.1: 1 point. Select the grid box at that point.
    std::fill(values, &(values[datasize]), 0);

    // Find the point
    int i, j;
    if (find_in_grid_var(points[0].x, cols, x_grid, i, true) &&
        find_in_grid_var(points[0].y, rows, y_grid, j, (y_grid[0] < y_grid[1]))) {

      values[(j * cols) + i] = 1;
    }
  } else if (numpoints > 2) {
    if (config.fringe_size > 0) {
      // Case 2.1: 3+ points, fringe set to > 0. Expand polygon and get coverage grid.
      int i;
      vector<Point> p = points;
      vector<Point> np(numpoints);
      Line lines[numpoints];

      // FIXME invalid assumption that these will be the same
      double box_size_x = x_grid[1] - x_grid[0];
      double box_size_y = y_grid[0] - y_grid[1];

      // Expand polygon  MPN: FIXME This is broken and blows up on high-res polygons because of line intersections disappearing into the polygon (or something like that and I don't have time to analyze it further).
      for (i = 0; i < numpoints; i++) {
        const int inext = (i + 1) % numpoints;

        // Calculate unit vectors
        const double len = sqrt(squared(p[i].x - p[inext].x) + squared(p[i].y - p[inext].y));
        const double vx = (p[i].x - p[inext].x) / len;
        const double vy = (p[i].y - p[inext].y) / len;

        // Add normal to normalized vector to points at ends of lines
        lines[i].from.x = p[i].x + vy * config.fringe_size * box_size_x;
        lines[i].from.y = p[i].y - vx * config.fringe_size * box_size_y;
        lines[i].to.x = p[inext].x + vy * config.fringe_size * box_size_x;
        lines[i].to.y = p[inext].y - vx * config.fringe_size * box_size_y;
      }

      // Define new points on basis of lines
      for (i = 0; i < numpoints; i++) {
        const int iprev = (i + numpoints - 1) % numpoints;

        // New point is intersect of line behind this one and this line
        np[i] = line_intersect(lines[iprev], lines[i]);

        // Fail if lines didn't intersect (should be impossible)
        assert(np[i].valid());

        const double len = sqrt(squared(p[i].x - np[i].x) + squared(p[i].y - np[i].y));
        const double vx = (np[i].x - p[i].x) / len;
        const double vy = (np[i].y - p[i].y) / len;
        const double fringe_len = sqrt(squared(config.fringe_size * box_size_x) + squared(config.fringe_size * box_size_y)) * 1.05;

        if (len > fringe_len) {
          np[i] = Point(p[i].x + vx * fringe_len, p[i].y + vy * fringe_len);
        }
      }
      // FIXME: Split up polygon along dividing line if necessary

      draw_polygon(rows, cols, r_grid, y_grid, x_grid, np);
    } else {
      // Case 2.2: 3+ points, fringe set to <= 0. Get coverage grid.
      draw_polygon(rows, cols, r_grid, y_grid, x_grid, points);
    }
    // Mask off all boxes which are outside the poly in the data (stipple) mask
    // Need to allow most-covered grid box to be selected if no others available
    for (int i = 0; i < rows; i++) {
      double area = fabs((y_grid[i * 2] - y_grid[i * 2 + 1]) * (x_grid[1] - x_grid[0]));
      int* val_row = &values[i * cols];
      double* r_row = &r_grid[i * cols];
      for (int j = 0; j < cols; j++) {

        // If coverage is less than a set threshold, throw this box out
        val_row[j] &= (r_row[j] / area > threshold);
      }
    }
  }

  delete[] r_grid;
  return g;
}

gdImagePtr get_image(string s) {
  FILE* infile = fopen(s.c_str(), "rb");
  if (!infile) {
    fprintf(stderr, "Couldn't read map file '%s'\n", s.c_str());
    exit(1);
  }
  gdImagePtr image = gdImageCreateFromPng(infile);
  if (!image) {
    fprintf(stderr, "Couldn't load map file '%s'\n", s.c_str());
    exit(1);
  }
  fclose(infile);
  return image;
}

void copy_partial_tile(gdImagePtr basemap, gdImagePtr image, Window& w, Window& tile_win, Window& intersect, const double size_x, const double size_y) {
  const int height = (int)ceil(((intersect.top - intersect.bottom) / size_y) * (double)image->sy);
  const int width = (int)ceil(((intersect.right - intersect.left) / size_x) * (double)image->sx);
  const int dstX = (int)floor(((intersect.left - w.left) / size_x) * (double)image->sx);
  const int dstY = (int)floor(((w.top - intersect.top) / size_y) * (double)image->sy);
  const int srcX = max(0, (int)floor(((intersect.left - tile_win.left) / size_x) * (double)image->sx));
  const int srcY = max(0, (int)floor(((tile_win.top - intersect.top) / size_y) * (double)image->sy));
  fprintf(stderr, "dstx: %i, dsty: %i, srcx: %i, srcy: %i, width: %i, height: %i, img_x: %i, img_y: %i\n", dstX, dstY, srcX, srcY, width, height, image->sx, image->sy);
  gdImageCopy(basemap, image, dstX, dstY, srcX, srcY, width, height);
}

gdImagePtr DataManager::get_basemap_image(Window& w) {
  gdImagePtr basemap = 0;
  const double size_x = (config.lowerright.x - config.upperleft.x) / config.zoom_factor;
  const double size_y = (config.upperleft.y - config.lowerright.y) / config.zoom_factor;

  const int min_x = (int)floor((w.left - config.upperleft.x) / size_x);
  const int min_y = (int)floor((config.upperleft.y - w.top) / size_y);
  const int max_x = (int)ceil((w.right - config.upperleft.x) / size_x);
  const int max_y = (int)ceil((config.upperleft.y - w.bottom) / size_y);
  assert(min_y < max_y);
  assert(min_x < max_x);
  fprintf(stderr, "Region: (%f, %f) to (%f, %f)\n", config.upperleft.x, config.upperleft.y, config.lowerright.x, config.lowerright.y);
  fprintf(stderr, "window: %f, %f, %f, %f\n", w.top, w.bottom, w.left, w.right);
  for (int y = min_y; y < max_y; y++) {
    for (int x = min_x; x < max_x; x++) {
      fprintf(stderr, "x: %i, y: %i\n", x, y);
      stringstream s;
      s << config.map_dir << config.region << "_" << config.zoom_factor << "_" << config.resolution << "_" << (x + config.zoom_factor) % config.zoom_factor << "_" << y % config.zoom_factor << ".png";

      gdImagePtr image = get_image(s.str());
      assert(image);

      if (!basemap) {
        basemap = gdImageCreate(image->sx, image->sy);

        // Ensure that colours get allocated in proper order
        gdImageColorAllocate(basemap, 255, 255, 255);
        gdImageColorAllocate(basemap, 0, 0, 0);
      }

      Window tile_win(config.upperleft.y - y * size_y, config.upperleft.x + x * size_x,
                      config.upperleft.y - (y + 1) * size_y, config.upperleft.x + (x + 1) * size_x);
      Window* intersect = tile_win.intersect(w);

      //fprintf(stderr, "window: %f, %f, %f, %f\n", w.top, w.bottom, w.left, w.right);
      fprintf(stderr, "tile window: %f, %f, %f, %f\n", tile_win.top, tile_win.bottom, tile_win.left, tile_win.right);
      fprintf(stderr, "intersect: %f, %f, %f, %f\n", intersect->top, intersect->bottom, intersect->left, intersect->right);

      copy_partial_tile(basemap, image, w, tile_win, *intersect, size_x, size_y);

      delete intersect;
      gdImageDestroy(image);
    }
  }
  return (basemap);
}

gdImagePtr DataManager::get_basemap() {
  return get_basemap_image(w);
}

DataGrid<int> DataManager::get_slmask(const DataSpec& s) {
  if (use_alt_baseline()) {
    DataSpec s2(config.baseline_model, config.baseline_expt, s.timeslice, s.variable, s.timeofyear, s.anom, s.percent_change);
    DataGrid<int> g(s);
    set_up_grid(g, SLMASK_SLOT, use_alt_baseline());
    int* values = g.values().get();
    assert(bp.get_windowed_mask_data(s2, SLMASK_SLOT, values));
    return g;
  } else {
    DataGrid<int> g(s);
    set_up_grid(g, SLMASK_SLOT);
    int* values = g.values().get();
    assert(dp.get_windowed_mask_data(s, SLMASK_SLOT, values));
    return g;
  }
}

double* DataProvider::get_x_grid(const DataSpec& s, const list<Window>& l) {
  const int total_width = sub_x_size(s, l);
  assert(total_width);

  NcVar* nc_grid_longs = get_ncdf_var(s, LON_BNDS_SLOT);
  double* values = new double[total_width * 2];
  std::string projection = get_projection(s);

  int k = 0;
  if (nc_grid_longs) {
    for (list<Window>::const_iterator i = l.begin(); i != l.end(); ++i) {
      const Window& w = *i;
      const int width = w.right - w.left + 1;
      assert(nc_grid_longs->set_cur(w.left, 0));
      assert(nc_grid_longs->get(&values[k], width, 2));

      // If we cross 180 degrees longitude, the sign will flip; correct this
      if (projection == "equidistant_cylindrical" && values[0] > values[1])
        values[0] -= 360.0;

      k += width * 2;
    }
  } else {
    NcVar* nc_x;
    double* xdat = new double[total_width];
    if (projection == "equidistant_cylindrical") {
      nc_x = get_ncdf_var(s, LONGS_SLOT);
    } else {
      nc_x = get_ncdf_var(s, XC_SLOT);
    }
    assert(nc_x);

    for (list<Window>::const_iterator i = l.begin(); i != l.end(); ++i) {
      const Window& w = *i;
      const int width = w.right - w.left + 1;

      assert(nc_x->set_cur(w.left));
      assert(nc_x->get(xdat, width));

      values[k] = xdat[0] - ((xdat[1] - xdat[0]) / 2);
      values[k + 1] = xdat[0] + ((xdat[1] - xdat[0]) / 2);
      for (int i = 1; i < width; i++) {
        values[k + i * 2] = values[i * 2 - 1];
        values[k + i * 2 + 1] = xdat[i] + ((xdat[i] - xdat[i - 1]) / 2);
      }

      k += width * 2;
    }

    delete[] xdat;
  }
  return values;
}

double* DataProvider::get_y_grid(const DataSpec& s, const list<Window>& l) {
  NcVar* nc_grid_lats = get_ncdf_var(s, LAT_BNDS_SLOT);
  const Window& w = *(l.begin());
  const int height = w.bottom - w.top + 1;
  assert(height);

  double* values = new double[height * 2];

  if (nc_grid_lats) {
    // If we have the actual grid edges, then use them
    assert(nc_grid_lats->set_cur(w.top, 0));
    assert(nc_grid_lats->get(values, height, 2));
    for (int i = 0; i < height; i++) {
      swap(values[i * 2], values[i * 2 + 1]);
    }
  } else {
    NcVar* nc_y;
    double* ydat = new double[height];
    if (get_projection(s) == "equidistant_cylindrical") {
      nc_y = get_ncdf_var(s, LATS_SLOT);
    } else {
      nc_y = get_ncdf_var(s, YC_SLOT);
    }

    assert(nc_y);
    assert(nc_y->set_cur(w.top));
    assert(nc_y->get(ydat, height));

    values[0] = ydat[0] + ((ydat[0] - ydat[1]) / 2);
    values[1] = ydat[0] - ((ydat[0] - ydat[1]) / 2);
    for (int i = 1; i < height; i++) {
      values[i * 2] = values[i * 2 - 1];
      values[i * 2 + 1] = ydat[i] - ((ydat[i - 1] - ydat[i]) / 2);
    }
    delete[] ydat;
  }
  return values;
}

DataGrid<double> DataManager::get_data(const DataSpec& s) {

  //fprintf(stderr, "dm::get_data\n");
  const bool alt_baseline = use_alt_baseline();

  // Grab in all the goodies I want to use
  DataGrid<double> g(s);
  set_up_grid(g, DATA_SLOT, alt_baseline);
  DataSpec ns(s.model, s.expt, g.base_period(), s.variable, s.timeofyear, ABSOLUTE);
  double* values = g.values().get();
  if (alt_baseline) {
    std::fill(values, &values[g.x_size() * g.y_size()], g.missing());
    // Grab in original data.
    DataGrid<double> gprime(s);
    set_up_grid(gprime, DATA_SLOT);
    double* input_values = gprime.values().get();
    assert(dp.get_windowed_data(s, DATA_SLOT, input_values));

    // If it's a percentage anomaly, we have to convert it to a non-percentage anomaly by multiplying by the baseline, interpolate the non-percentage anomaly and the baseline, then divide interpolated non-percent anomaly by the baseline. Yay.
    if (gprime.r.anom == ANOMALY && gprime.r.percent_change) {
      fprintf(stderr, "Alt baseline pct change code activated\n");

      // Get baseline data at input data resolution.
      DataGrid<double> gprime_base(ns);
      set_up_grid(gprime_base, DATA_SLOT);
      double* input_values_base = gprime_base.values().get();
      assert(dp.get_windowed_data(ns, DATA_SLOT, input_values_base));

      fprintf(stderr, "grid size: %i\n", gprime.grid_size());

      // Convert percentage anomaly data to non-percent anomaly data.
      for (int i = 0; i < gprime.grid_size(); ++i) {
        fprintf(stderr, "(%f,%f) ", input_values[i], input_values_base[i]);
        input_values[i] = (input_values[i] / 100) * input_values_base[i];
      }
      fprintf(stderr, "\n");

      // Set up baseline data to interpolate grid to...
      DataGrid<double> g_base(ns);
      set_up_grid(g_base, BASELINE_SLOT, alt_baseline);
      double* values_base = g_base.values().get();
      std::fill(values_base, &values_base[g_base.x_size() * g_base.y_size()], g_base.missing());

      // Interpolate both baseline data and (non-percent) anomaly data to alt baseline res.
      interpolate_data(gprime, g);
      interpolate_data(gprime_base, g_base);

      // Convert non-percent anomaly data to percent at new new res.
      for (int i = 0; i < g.grid_size(); ++i)
        values[i] = (values[i] / values_base[i]) * 100;

    } else {
      interpolate_data(gprime, g);
    }
  } else {
    assert(dp.get_windowed_data(s, DATA_SLOT, values));
  }

  // If we need a baseline, fetch one
  DataGrid<double> bg(ns);
  double* baseline = 0;
  if ((g.anomaly() && s.anom == ABSOLUTE) || (!g.anomaly() && s.anom == ANOMALY)) {
    if (alt_baseline) {
      set_up_grid(bg, DATA_SLOT, alt_baseline);
      baseline = bg.values().get();
      if (g.anomaly() && s.anom == ABSOLUTE) {
        // If we want absolute numbers, add the alternate baseline
        DataSpec base_spec(config.baseline_model, config.baseline_expt, g.base_period(), s.variable, s.timeofyear, ABSOLUTE);
        double missing = bp.get_missing<double>(base_spec, BASELINE_SLOT);
        bg.set_data(bg.x_size(), bg.y_size(), missing, bg.values(), bg.x_grid(), bg.y_grid(), bg.proj4_string(), bg.projection(), bg.anomaly(), bg.base_period());
        assert(bp.get_windowed_data(base_spec, BASELINE_SLOT, baseline));
      } else {
        // If we want an anomaly, remove the model's baseline, but at the alternate baseline's resolution
        std::fill(baseline, &baseline[bg.x_size() * bg.y_size()], bg.missing());
        DataGrid<double> bgprime(ns);
        set_up_grid(bgprime, DATA_SLOT);
        double* input_values = bgprime.values().get();
        assert(dp.get_windowed_data(ns, BASELINE_SLOT, input_values));
        interpolate_data(bgprime, bg);
      }
    } else {
      set_up_grid(bg, DATA_SLOT);
      baseline = bg.values().get();
      assert(dp.get_windowed_data(ns, BASELINE_SLOT, baseline));
    }
  }
  double data_missing = g.missing();
  double baseline_missing = bg.missing();
  //  fprintf(stderr, "data missing: %f, baseline missing: %f\n", data_missing, baseline_missing);

  // If desired anomaly expressed, modify data to desired anomaly type
  if (g.anomaly() && s.anom == ABSOLUTE) {
    if (s.percent_change)
      for (int i = 0; i < g.grid_size(); i++) {
        if (baseline[i] == baseline_missing || values[i] == data_missing) {
          values[i] = data_missing;
        } else {
          if (baseline[i] == 0)
            values[i] = 0;
          else
            values[i] = (values[i] / 100) * baseline[i] + baseline[i];
        }
      }
    else
      for (int i = 0; i < g.grid_size(); i++)
        if (baseline[i] == baseline_missing || values[i] == data_missing) {
          values[i] = data_missing;
        } else
          values[i] += baseline[i];
  } else if (!g.anomaly() && s.anom == ANOMALY) {
    if (s.percent_change) {
      for (int i = 0; i < g.grid_size(); i++)
        if (baseline[i] == baseline_missing || values[i] == data_missing) {
          values[i] = data_missing;
        } else
          values[i] = ((values[i] - baseline[i]) / baseline[i]) * 100;
    } else {
      for (int i = 0; i < g.grid_size(); i++)
        if (baseline[i] == baseline_missing || values[i] == data_missing) {
          values[i] = data_missing;
        } else
          values[i] -= baseline[i];
    }
  }

  // Set anomaly to appropriate value
  if ((g.anomaly() && s.anom == ABSOLUTE) || (!g.anomaly() && s.anom == ANOMALY))
    g.set_data(g.x_size(), g.y_size(), g.missing(), g.values(), g.x_grid(), g.y_grid(), g.proj4_string(), g.projection(), s.anom, g.base_period());

  return (g);
}

bool DataProvider::get_windowed_data(const DataSpec& s, SLOT slot, double* values) {
  int k = 0;
  const int out_width = sub_x_size(s);
  NcVar* data = get_ncdf_var(s, slot);

  // Special case size=1 to avoid copy penalty
  if (l.size() == 1) {
    assert(get_slice(values, data, s.timeofyear, *(l.begin())));
  } else {
    for (list<Window>::const_iterator i = l.begin(); i != l.end(); ++i) {
      const Window& w = *i;
      const int height = w.bottom - w.top + 1;
      const int width = w.right - w.left + 1;
      double* temp = new double[height * width];

      assert(get_slice(temp, data, s.timeofyear, w));

      for (int j = 0; j < height; j++) {
        std::copy(&temp[j * width], &temp[(j + 1) * width], &values[j * out_width + k]);
      }
      delete[] temp;
      k += width;
    }
  }
  return true;
}

bool DataProvider::get_windowed_mask_data(const DataSpec& s, SLOT slot, int* values) {
  int k = 0;
  const int out_width = sub_x_size(s);
  NcVar* data = get_ncdf_var(s, slot);

  // Special case size=1 to avoid copy penalty
  if (l.size() == 1) {
    const Window& w = *(l.begin());
    const int height = w.bottom - w.top + 1;
    const int width = w.right - w.left + 1;
    assert(data->set_cur(w.top, w.left));
    assert(data->get(values, height, width));
  } else {
    for (list<Window>::const_iterator i = l.begin(); i != l.end(); ++i) {
      const Window& w = *i;
      const int height = w.bottom - w.top + 1;
      const int width = w.right - w.left + 1;
      int* temp = new int[height * width];
      assert(data->set_cur(w.top, w.left));
      assert(data->get(temp, height, width));

      for (int j = 0; j < height; j++) {
        std::copy(&temp[j * width], &temp[(j + 1) * width], &values[j * out_width + k]);
      }
      delete[] temp;
      k += width;
    }
  }
  return true;
}

bool DataProvider::get_slice(double* values, NcVar* data, const int timeofyear, const Window& w) {
  assert(data->set_cur(timeofyear, w.top, w.left));
  float* fvals;
  const int xsize = w.right - w.left + 1;
  const int ysize = w.bottom - w.top + 1;
  const int datasize = xsize * ysize;
  switch (data->type()) {
  case ncDouble:
    assert(data->get(values, 1, ysize, xsize));
    break;
  case ncFloat:
    fvals = new float[datasize];
    assert(data->get(fvals, 1, ysize, xsize));
    copy(fvals, &fvals[datasize], values);
    delete[] fvals;
    break;
  default:
    fprintf(stderr, "Unhandled input data type\n");
    assert(false);
    break;
  }
  return true;
}

#define NUM_WIN 3
list<Window> DataProvider::get_subsets(const DataSpec& s) {
  list<Window> lw;
  list<int> l[NUM_WIN];
  list<int> h;
  int i;
  const int xsize = x_size(s);
  const int ysize = y_size(s);
  const int xgrid_size = xsize * 2;
  const int ygrid_size = ysize * 2;

  list<Window> tmp;
  tmp.push_back(Window(0, 0, ysize - 1, xsize - 1));
  fprintf(stderr, "x,y sizes: %i, %i\n", xsize, ysize);
  const double* x_grid = get_x_grid(s, tmp);
  const double* y_grid = get_y_grid(s, tmp);
  //std::copy(g.x_grid(), &g.x_grid()[x_size], x_grid);

  if (get_projection(s) == "equidistant_cylindrical") {
    fprintf(stderr, "sdm::get_subsets: About to shift longs\n");

    // Accumulate lists of offsets of gridboxes that are within the window
    // FIXME: Should be replaced by some form of binary search; this is inefficient
    // and this will be done more than once per data file with new architecture
    int j = 0;
    double k = 0;
    for (i = 0; i < xgrid_size && j < 3; i += 2) {
      if ((w.left + k <= x_grid[i] && x_grid[i] < w.right + k) ||
          (w.left + k <= x_grid[i + 1] && x_grid[i + 1] < w.right + k) ||
          (x_grid[i] <= w.left + k && w.left + k < x_grid[i + 1])) {
        l[j].push_back(i / 2);
      }
      //if((x_grid[i] <= w.left && w.left < x_grid[i + 1]) ||
      // (x_grid[i] <= w.right && w.right < x_grid[i + 1])) {
      if (x_grid[i] <= w.right && w.right < x_grid[i + 1]) {
        j++;
        k += 360;
      }
    }

    // Find the latitudes within this window
    for (i = 0; i < ygrid_size; i += 2) {
      //if((w.bottom <= y_grid[i] && w.top > y_grid[i]) ||
      // (w.bottom <= y_grid[i + 1] && w.top > y_grid[i + 1])) {
      if ((w.bottom <= y_grid[i] && w.top > y_grid[i]) ||
          (w.bottom <= y_grid[i + 1] && w.top > y_grid[i + 1]) ||
          (y_grid[i + 1] <= w.bottom && w.bottom < y_grid[i])) {
        h.push_back(i / 2);
      }
    }

    assert(h.size() > 0);

    // Add windows into the data to the list
    for (i = NUM_WIN - 1; i >= 0; i--) {
      if (l[i].size() > 0) {
        lw.push_back(Window(*(min_element(h.begin(), h.end())), *(min_element(l[i].begin(), l[i].end())),
                            *(max_element(h.begin(), h.end())), *(max_element(l[i].begin(), l[i].end()))));
        fprintf(stderr, "Added window\n");
      }
    }
  } else {
    int top, bottom, left, right;
    if (!find_in_grid_var(w.top, ysize, y_grid, top, (y_grid[0] < y_grid[1])))
      top = 0;
    if (!find_in_grid_var(w.bottom, ysize, y_grid, bottom, (y_grid[0] < y_grid[1])))
      bottom = ysize - 1;
    if (!find_in_grid_var(w.left, xsize, x_grid, left, (x_grid[0] < x_grid[1])))
      left = 0;
    if (!find_in_grid_var(w.right, xsize, x_grid, right, (x_grid[0] < x_grid[1])))
      right = xsize - 1;

    if (top > bottom) swap(top, bottom);
    if (left > right) swap(left, right);

    lw.push_back(Window(top, left, bottom, right));
  }

  if (used_for_bilinear) {
    Window& first = lw.front();
    Window& last = lw.back();
    int new_top = (first.top == 0) ? 0 : first.top - 1;
    int new_bottom = (first.bottom == (xsize - 1)) ? xsize - 1 : first.bottom + 1;

    for (list<Window>::iterator i = lw.begin(); i != lw.end(); ++i) {
      Window& w = *i;
      w.top = new_top;
      w.bottom = new_bottom;
    }

    if (first.left == 0) {
      lw.push_front(Window(new_top, xsize - 1, new_bottom, xsize - 1));
    } else {
      --first.left;
    }

    if (last.right == xsize - 1) {
      lw.push_back(Window(new_top, 0, new_bottom, 0));
    } else {
      ++last.right;
    }
  }

  delete[] x_grid;
  delete[] y_grid;

  return lw;
}

