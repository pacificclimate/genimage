#include "support.h"
#include "line.h"

enum CONN_STATES { UNDEF, USED, LINE, TEMP };

vector<double> quantile(const vector<double> indata, const vector<double> quantiles) {
  vector<double> result;
  if (quantiles.size() == 0)
    return result;

  const int n = indata.size();

  // Make a copy because nth_element modifies the data.
  vector<double> data = indata;

  // Constants for quantiles. Can be modified if needed.
  const double a = 1.0 / 3.0;
  const double b = 1.0 / 3.0;
  const double fuzz = 4 * std::numeric_limits<double>::epsilon();

  for (vector<double>::const_iterator quantile_iter = quantiles.begin(); quantile_iter != quantiles.end(); ++quantile_iter) {
    if (n < 2) {
      result.push_back(std::numeric_limits<double>::quiet_NaN());
      continue;
    }
    const double quantile = *quantile_iter;

    const double nppm = a + quantile * (n + 1 - a - b) - 1;
    const size_t j = (size_t)std::max(0.0, floor(nppm + fuzz));

    // Variance from R: Should probably be <= not < here.
    const double h = (fabs(nppm - (double)j) <= fuzz) ? 0 : nppm - (double)j;
    size_t right_elem = std::max(0ul, std::min(j + 1, data.size() - 1));
    size_t left_elem = std::max(0ul, std::min(j, data.size() - 1));

    if (h == 1) {
      std::nth_element(data.begin(), data.begin() + right_elem, data.end());
      result.push_back((data[right_elem]));
    } else {
      // No guarantee that 2nd nth_element call will preserve order such that the pointer used by the 1st call still points to the same thing; so store the result before calling nth_element again.
      std::nth_element(data.begin(), data.begin() + left_elem, data.end());
      const double left = data[left_elem];
      if (h == 0) {
        result.push_back(left);
      } else {
        std::nth_element(data.begin(), data.begin() + right_elem, data.end());
        const double right = data[right_elem];
        result.push_back((1 - h) * left + h * right);
      }
    }
  }
  return result;
}

// Affects the input data
void shift_longs(double* longs, const int size, const double x_center) {
  for (int i = 0; i < size; i++)
    longs[i] = fmod(longs[i] - x_center + 540, 360) + x_center - 180;
  if (longs[0] == longs[size - 1]) {
    if ((fabs(longs[size - 1] - longs[size - 2]) + 1) < fabs(longs[0] - longs[1])) {
      longs[0] -= 360;
    } else {
      longs[size - 1] += 360;
    }
  }
}

bool find_in_range(double num, int max, const double* vals, int& ret, bool ascending) {
  if (ascending) {
    for (ret = 0; ret < max; ret++)
      if (vals[ret] <= num && vals[ret + 1] > num)
        break;
  } else {
    for (ret = 0; ret < max; ret++)
      if (vals[ret + 1] < num && vals[ret] >= num)
        break;
  }
  return (ret != max);
}

// Finds a lat/lon in a grid variable
bool find_in_grid_var(double num, int max, const double* vals, int& ret, bool ascending) {
  if (ascending) {
    for (ret = 0; ret < max; ret++)
      if (vals[ret * 2] <= num && num < vals[ret * 2 + 1])
        break;
  } else {
    for (ret = 0; ret < max; ret++)
      if (vals[ret * 2 + 1] < num && num <= vals[ret * 2])
        break;
  }
  return (ret != max);
}

double poly_area(const vector<Point>& points) {
  double area = 0;
  if (points.size() == 0) return (0.0);

  Point s = points[points.size() - 1];
  for (vector<Point>::const_iterator p = points.begin(); p != points.end(); ++p) {
    const Point &e = *p;
    area += (s.x + e.x) * (s.y - e.y);
    s = e;
  }

  return area * .5;
}

void split_on_x(const double x, const vector<Point>& in, vector<Point>& left, vector<Point>& right) {
  Point s = in[in.size() - 1];
  for (unsigned int j = 0; j < in.size(); j++) {
    const Point& e = in[j];
    if (e.x > x) { // ends to right of line
      const double intersect = s.y + ( ((e.y - s.y) / (e.x - s.x)) * (x - s.x) );
      if (s.x < x) { // started to the left of the line
        left.push_back(Point( x, intersect ));
        right.push_back(Point( x, intersect ));
      }
      right.push_back(e);
    } else if (e.x == x) { // ends ON boundary, therefore NO intersection.
      left.push_back(e);
      right.push_back(e);
    } else if (s.x > x) { // started to right of line, ended to left of line
      const double intersect = s.y + ( ((e.y - s.y) / (e.x - s.x)) * (x - s.x) );
      right.push_back(Point( x, intersect ));
      left.push_back(Point( x, intersect ));
      left.push_back(e);
    } else { // started and ended to left of line
      left.push_back(e);
    }
    s = e;
  }
}

void split_on_y(const double y, const vector<Point>& in, vector<Point>& top, vector<Point>& bottom) {
  Point s = in[in.size() - 1];
  for (unsigned int j = 0; j < in.size(); j++) {
    const Point& e = in[j];
    if (e.y > y) { // ends above the line
      const double intersect = s.x + ( ((e.x - s.x) / (e.y - s.y)) * (y - s.y) );
      if (s.y < y) { // started below the line
        top.push_back(Point( intersect, y ));
        bottom.push_back(Point( intersect, y ));
      }
      top.push_back(e);
    } else if (e.y == y) { // ends ON boundary, therefore NO intersection.
      top.push_back(e);
      bottom.push_back(e);
    } else if (s.y > y) { // started above the line, ended below the line
      const double intersect = s.x + ( ((e.x - s.x) / (e.y - s.y)) * (y - s.y) );
      top.push_back(Point( intersect, y ));
      bottom.push_back(Point( intersect, y ));
      bottom.push_back(e);
    } else { // started and ended below the line
      bottom.push_back(e);
    }
    s = e;
  }
}

void partition_recurse_y(const vector<Point>& in, const int x_off, const int y_off, const int x_size, const int y_size, const int row_length, double* grid, const double* x_grid, const double* y_grid, int lev);
void partition_recurse_x(const vector<Point>& in, const int x_off, const int y_off, const int x_size, const int y_size, const int row_length, double* grid, const double* x_grid, const double* y_grid, int lev) {
  // Case for no coverage
  if (in.size() == 0)
    return;

  // Base case for recursion
  if (x_size == 1 && y_size == 1) {
    grid[(y_off * row_length) + x_off] = fabs(poly_area(in));
    return;
  }

  // Keep recursing on other dimension if size of dimension we are recursing on is 1
  if (x_size == 1) {
    partition_recurse_y(in, x_off, y_off, x_size, y_size, row_length, grid, x_grid, y_grid, lev + 1);
    return;
  }

  if (in.size() == 4) {
    // If the diagonals of the polygon and of the grid box match, then they must be the same
    const double hyplen1 = pow(in[0].x - in[2].x, 2) + pow(in[0].y - in[2].y, 2);
    const double hyplen2 = pow(in[1].x - in[3].x, 2) + pow(in[1].y - in[3].y, 2);
    const double boxhyp = pow(x_grid[x_off * 2] - x_grid[(x_off + x_size) * 2 - 1], 2) + pow(y_grid[y_off * 2] - y_grid[(y_off + y_size) * 2 - 1], 2);
    if (fabs(boxhyp / hyplen2) - 1 < 1e-6 && fabs(boxhyp / hyplen1) - 1 < 1e-6) {
      //if(boxhyp == hyplen2 && boxhyp == hyplen1) {
      // If diagonals match, they're equal; so fill it in.
      for (int y = y_off; y < y_off + y_size; y++)
        std::fill(&grid[(y * row_length) + x_off], &grid[(y * row_length) + x_off + x_size], fabs((y_grid[y * 2] - y_grid[y * 2 + 1]) * (x_grid[x_off * 2] - x_grid[x_off * 2 + 1])));

      return;
    }
  }

  const int center_x = x_off + x_size / 2;
  const int left_size = x_size / 2;
  const int right_size = (x_size + 1) / 2;
  const double center_x_value = x_grid[center_x * 2];
  vector<Point> left, right;
  split_on_x(center_x_value, in, left, right);
  //fprintf(stderr, "level %i: in points: %li, left points: %li, right points: %li, center: %f, center index: %i, x_size: %i, y_size: %i, x_off: %i, y_off: %i\n", lev, in.size(), left.size(), right.size(), center_x_value, center_x, x_size, y_size, x_off, y_off);
  partition_recurse_y(left, x_off, y_off, left_size, y_size, row_length, grid, x_grid, y_grid, lev + 1);
  partition_recurse_y(right, x_off + left_size, y_off, right_size, y_size, row_length, grid, x_grid, y_grid, lev + 1);
}

void partition_recurse_y(const vector<Point>& in, const int x_off, const int y_off, const int x_size, const int y_size, const int row_length, double* grid, const double* x_grid, const double* y_grid, int lev) {
  // Case for no coverage
  if (in.size() == 0)
    return;

  // Base case for recursion
  if (x_size == 1 && y_size == 1) {
    grid[(y_off * row_length) + x_off] = fabs(poly_area(in));
    return;
  }

  // Keep recursing on other dimension if size of dimension we are recursing on is 1
  if (y_size == 1) {
    partition_recurse_x(in, x_off, y_off, x_size, y_size, row_length, grid, x_grid, y_grid, lev + 1);
    return;
  }

  if (in.size() == 4) {
    // If the diagonals of the polygon and of the grid box match, then they must be the same
    const double hyplen1 = pow(in[0].x - in[2].x, 2) + pow(in[0].y - in[2].y, 2);
    const double hyplen2 = pow(in[1].x - in[3].x, 2) + pow(in[1].y - in[3].y, 2);
    const double boxhyp = pow(x_grid[x_off * 2] - x_grid[(x_off + x_size) * 2 - 1], 2) + pow(y_grid[y_off * 2] - y_grid[(y_off + y_size) * 2 - 1], 2);
    if (fabs(boxhyp / hyplen2) - 1 < 1e-6 && fabs(boxhyp / hyplen1) - 1 < 1e-6) {
      //if(boxhyp == hyplen2 && boxhyp == hyplen1) {
      // If diagonals match, they're equal; so fill it in.
      for (int y = y_off; y < y_off + y_size; y++)
        std::fill(&grid[(y * row_length) + x_off], &grid[(y * row_length) + x_off + x_size], fabs((y_grid[y * 2] - y_grid[y * 2 + 1]) * (x_grid[x_off * 2] - x_grid[x_off * 2 + 1])));

      return;
    }
  }

  const int center_y = y_off + y_size / 2;
  const int top_size = y_size / 2;
  const int bottom_size = (y_size + 1) / 2;
  const double center_y_value = y_grid[center_y * 2];
  vector<Point> top, bottom;
  split_on_y(center_y_value, in, top, bottom);
  //fprintf(stderr, "level %i: in points: %li, top points: %li, bottom points: %li, center: %f, center index: %i, x_size: %i, y_size: %i, x_off: %i, y_off: %i\n", lev, in.size(), top.size(), bottom.size(), center_y_value, center_y, x_size, y_size, x_off, y_off);
  partition_recurse_x(top, x_off, y_off, x_size, top_size, row_length, grid, x_grid, y_grid, lev + 1);
  partition_recurse_x(bottom, x_off, y_off + top_size, x_size, bottom_size, row_length, grid, x_grid, y_grid, lev + 1);
}

// Cohen-Sutherland algorithm
// See http://www.cc.gatech.edu/grads/h/Hao-wei.Hsieh/Haowei.Hsieh/mm.html
// Also see http://en.wikipedia.org/wiki/Sutherland-Hodgman_clipping_algorithm
vector<Point> sutherland_clip(const double top, const double left, const double bottom, const double right, const vector<Point>& points) {
  assert(points.size() >= 3); // GRRR
  assert(top != bottom);
  assert(right != left);
  assert(top > bottom);
  assert(right > left);

  vector<Point> a, b;

  { // Left
    const vector<Point>& in = points;
    vector<Point>& out      = a;

    Point s = in[in.size() - 1];
    for (unsigned int j = 0; j < in.size(); j++) {
      const Point& e = in[j];
      if (e.x > left) { // ends fully inside
        if (s.x < left) { // started outside, therefore intersection
          double intersect = s.y + ( ((e.y - s.y) / (e.x - s.x)) * (left - s.x) );
          out.push_back(Point( left, intersect ));
        }
        out.push_back(e);
      } else if (e.x == left) { // ends ON boundary, therefore NO intersection.
        out.push_back(e);
      } else if (s.x > left) { // started fully inside (and did not end inside)
        double intersect = s.y + ( ((e.y - s.y) / (e.x - s.x)) * (left - s.x) );
        out.push_back(Point(   left, intersect ));
      }
      s = e;
    }

    if (out.size() == 0) return (out);
  }

  { // Right
    const vector<Point>& in = a;
    vector<Point>& out      = b;

    Point s = in[in.size() - 1];
    for (unsigned int j = 0; j < in.size(); j++) {
      const Point& e = in[j];
      if (e.x < right) { // ends fully inside
        if (s.x > right) { // started outside, therefore intersection
          double intersect = s.y + ( ((e.y - s.y) / (e.x - s.x)) * (right - s.x) );
          out.push_back(Point( right, intersect));
        }
        out.push_back(e);
      } else if (e.x == right) { // ends ON boundary, therefore NO intersection.
        out.push_back(e);
      } else if (s.x < right) { // started fully inside (and did not end inside)
        double intersect = s.y + ( ((e.y - s.y) / (e.x - s.x)) * (right - s.x) );
        out.push_back(Point( right, intersect));
      }
      s = e;
    }

    if (out.size() == 0) return (out);
  }

  { // Bottom
    const vector<Point>& in = b;
    vector<Point>& out      = a;
    out.clear();

    Point s = in[in.size() - 1];
    for (unsigned int j = 0; j < in.size(); j++) {
      const Point& e = in[j];
      if (e.y > bottom) { // ends fully inside
        if (s.y < bottom) { // started outside, therefore intersection
          double intersect = s.x + ( ((e.x - s.x) / (e.y - s.y)) * (bottom - s.y) );
          out.push_back(Point( intersect, bottom));
        }
        out.push_back(e);
      } else if (e.y == bottom) { // ends ON boundary, therefore NO intersection.
        out.push_back(e);
      } else if (s.y > bottom) { // started fully inside (and did not end inside)
        double intersect = s.x + ( ((e.x - s.x) / (e.y - s.y)) * (bottom - s.y) );
        out.push_back(Point(intersect, bottom));
      }
      s = e;
    }

    if (out.size() == 0) return (out);
  }

  { // Top
    const vector<Point>& in = a;
    vector<Point>& out      = b;
    out.clear();

    Point s = in[in.size() - 1];
    for (int unsigned j = 0; j < in.size(); j++) {
      const Point& e = in[j];
      if (e.y < top) { // ends fully inside
        if (s.y > top) { // started outside, therefore intersection
          double intersect = s.x + ( ((e.x - s.x) / (e.y - s.y)) * (top - s.y) );
          out.push_back(Point(intersect, top));
        }
        out.push_back(e);
      } else if (e.y == top) { // ends ON boundary, therefore NO intersection.
        out.push_back(e);
      } else if (s.y < top) { // started fully inside (and did not end inside)
        double intersect = s.x + ( ((e.x - s.x) / (e.y - s.y)) * (top - s.y) );
        out.push_back(Point(intersect, top));
      }
      s = e;
    }

    return (out);
  }
}

void draw_polygon(int y_size, int x_size, double* grid, const double* y_grid, const double* x_grid, const vector<Point>& points) {
  // ASSUMPTION: Non-overlapping polygons
  int min_y = 0, max_y = y_size - 1;
  int min_x = 0, max_x = x_size - 1;
  Range r;
  Range xr(x_grid, x_size * 2);
  Range yr(y_grid, y_size * 2);
  Range prx;
  Range pry;

  for (unsigned int i = 0; i < points.size(); i++) {
    prx.add(points[i].x);
    pry.add(points[i].y);
  }

  fprintf(stderr, "x grid min: %f, max: %f\n", xr.min(), xr.max());
  fprintf(stderr, "y grid min: %f, max: %f\n", yr.min(), yr.max());
  fprintf(stderr, "points x min: %f, max: %f\n", prx.min(), prx.max());
  fprintf(stderr, "points y min: %f, max: %f\n", pry.min(), pry.max());

  if (!find_in_grid_var(pry.min(), y_size, y_grid, min_y, (y_grid[0] < y_grid[1])))
    min_y = (y_grid[0] < y_grid[1]) ? 0 : y_size - 1;

  if (!find_in_grid_var(pry.max(), y_size, y_grid, max_y, (y_grid[0] < y_grid[1])))
    max_y = (y_grid[0] < y_grid[1]) ? y_size - 1 : 0;

  if (!find_in_grid_var(prx.min(), x_size, x_grid, min_x, (x_grid[0] < x_grid[1])))
    min_x = 0;

  if (!find_in_grid_var(prx.max(), x_size, x_grid, max_x, (x_grid[0] < x_grid[1])))
    max_x = x_size - 1;

  if (y_grid[min_y * 2] > y_grid[max_y * 2])
    assert(false);

  // Note: min_y and max_y are the indices of the minimum and maximum VALUES of Y, not (necessarily) the minimum and maximum indices covered by the polygon.
  fprintf(stderr, "xmin: %i, xmax: %i, ymin: %i, ymax: %i\n", min_x, max_x, min_y, max_y);
  const double top = (y_grid[max_y * 2] > y_grid[max_y * 2 + 1]) ? y_grid[max_y * 2] : y_grid[max_y * 2 + 1];
  const double bottom = (y_grid[min_y * 2] > y_grid[min_y * 2 + 1]) ? y_grid[min_y * 2 + 1] : y_grid[min_y * 2];
  fprintf(stderr, "xmin_lon: %f, xmax_lon: %f, ymin_lat: %f, ymax_lat: %f\n", x_grid[min_x * 2], x_grid[max_x * 2 + 1], bottom, top);
  const vector<Point> cp = sutherland_clip(top, x_grid[min_x * 2], bottom, x_grid[max_x * 2 + 1], points);
  if (min_y < max_y) {
    partition_recurse_x(cp, min_x, min_y, max_x - min_x + 1, max_y - min_y + 1, x_size, grid, x_grid, y_grid, 0);
  } else {
    partition_recurse_x(cp, min_x, max_y, max_x - min_x + 1, min_y - max_y + 1, x_size, grid, x_grid, y_grid, 0);
  }

  /*
  for(int j = min_y; j <= max_y; j++) {
    for(int i = min_x; i <= max_x; i++) {
      // cohen_sutherland_clip takes top, left, bottom, right, points.
      grid[(j * x_size) + i] = fabs(poly_area(sutherland_clip(top, x_grid[i * 2], bottom, x_grid[i * 2 + 1], points)));
      r.add(grid[(j * x_size) + i]);
    }
  }
  fprintf(stderr, "grid area min: %f, max: %f\n", r.min(), r.max());
  */
}
