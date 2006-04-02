#ifndef __SUPPORT_H
#define __SUPPORT_H
#include <iostream>
#include <list>
#include "point.h"

using namespace std;

#define QUAD_SIDES 4

enum DIRECTIONS { LEFT, RIGHT };

void draw_triangle(int rows, int cols, double* grid, const double* lats, const double* lons, Point p1, Point p2, Point p3);

void draw_polygon(int rows, int cols, double* grid, const double* lats, const double* lons, Point **points, int numlines);

double clip_precision(double number, int num_mantissa_bits);
double clip_precision(double number, int num_mantissa_bits);

inline double squared(double number) {
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
// Must write function for spherical area
double triangle_area(Point p1, Point p2, Point p3);

class WPoint {
 public:
  WPoint() { p = 0; w = 0; }
  WPoint(double p, double w) { this->p = p; this->w = w; }
  WPoint(const WPoint& pt) { p = pt.p; w = pt.w; }
  bool operator < (const WPoint& pt) { return p < pt.p; }
  bool operator > (const WPoint& pt) { return p > pt.p; }
  bool operator <= (const WPoint& pt) { return p <= pt.p; }
  bool operator >= (const WPoint& pt) { return p >= pt.p; }
  bool operator == (const WPoint& pt) { return p == pt.p; }
  double p, w;
};

#endif
