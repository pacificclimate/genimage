#ifndef __SUPPORT_H
#define __SUPPORT_H
#include <iostream>
#include <list>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <assert.h>
#include "point.h"
#include "range.h"

using namespace std;

#define QUAD_SIDES 4

enum DIRECTIONS { LEFT, RIGHT };

// Compute quantiles.
vector<double> quantile(const vector<double> indata, const vector<double> quantiles);

// Alters longitudes to not wrap around
void shift_longs(double* longs, const int size, const double x_center);

// Finds a value between 2 numbers in a sorted range
bool find_in_range(double num, int max, const double* vals, int& ret, bool ascending);

bool find_in_grid_var(double num, int max, const double* vals, int& ret, bool ascending);

void draw_triangle(int rows, int cols, double* grid, const double* lats, const double* lons, Point p1, Point p2, Point p3);

void draw_polygon(int rows, int cols, double* grid, const double* lats, const double* lons, const vector<Point>& points);

double clip_precision(double number, int num_mantissa_bits);
double clip_precision(double number, int num_mantissa_bits);

inline double squared(const double number) {
  return number * number;
}

inline int sign(double number) {
  if(number == 0) {
    return 0;
  } else {
    return (int)(number / fabs(number));
  }
}

//bool compute_outline_change(list<Line >& outline, const Line& polyline, const Line& boxline, Point topright, Point bottomleft, bool& inside);

bool within_bbox(Point p, 
		 Point topright, 
		 Point bottomleft);

// Returns area of (admittedly non-spherical) triangle
double triangle_area(Point p1, Point p2, Point p3);

class WPoint {
 public:
  WPoint() { p = 0; w = 0; }
  WPoint(double p, double w): p(p), w(w) { }
  bool operator < (const WPoint& pt) const __attribute__ ((pure)) { return p < pt.p; }
  bool operator > (const WPoint& pt) const { return p > pt.p; }
  bool operator <= (const WPoint& pt) const { return p <= pt.p; }
  bool operator >= (const WPoint& pt) const { return p >= pt.p; }
  bool operator == (const WPoint& pt) const { return p == pt.p; }
  double p, w;
};

inline double tick_spacing(const Range& range, const int desired_ticks) {
  const double desired_min = log10(desired_ticks);
  double r = range.range();
  if(r == 0) {
    r = 0.01;
  }
  const double range_log = log10(r);
  const double base = floor(range_log - desired_min);
  const double rb_diff = range_log - base;
  const double log2 = log10(2);
  const double log5 = log10(5);
  
  // Figure out how to get the right # of tickmarks on the graph
  // Comparison against -1E12 is to compensate for error
  if(rb_diff - (log5 + desired_min) > -1E-12) {
    return pow(10, base) * 5;
  } else if(rb_diff - (log2 + desired_min) > -1E-12) {
    return pow(10, base) * 2;
  } else {
    return pow(10, base);
  }
}

#endif
