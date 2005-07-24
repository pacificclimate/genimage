#ifndef __SUPPORT_H
#define __SUPPORT_H
#include <iostream>
#include <list>
#include "line.h"

using namespace std;

#define QUAD_SIDES 4

enum DIRECTIONS { LEFT, RIGHT };

void draw_triangle(int rows, int cols, float* grid, float* lats, float* lons, Point<float> p1, Point<float> p2, Point<float> p3);

void draw_polygon(int rows, int cols, float* grid, float* lats, float* lons, Point<float> **points, int numlines);

Point<float>* 
ls_intersect(
const Line<Point<float> >& l1, 
const Line<Point<float> >& l2, 
bool ignore_first_range = false);

Point<float>* ls_intersect(const Point<float>& p00, const Point<float>& p01, const Point<float>& p10, const Point<float>& p11, bool ignore_first_range = false);

float clip_precision(float number, int num_mantissa_bits);
double clip_precision(double number, int num_mantissa_bits);

inline double squared(double number) {
  return number * number;
}

inline float squared(float number) {
  return number * number;
}

//bool compute_outline_change(list<Line<Point<float> > >& outline, const Line<Point<float> >& polyline, const Line<Point<float> >& boxline, Point<float> topright, Point<float> bottomleft, bool& inside);

bool within_bbox(Point<float> p, 
		 Point<float> topright, 
		 Point<float> bottomleft);

//list<Line<Point<float> > > get_box_outline(Point<float> topright, Point<float> bottomleft);
Line<Point<float> >* get_box_outline(Point<float> topright, Point<float> bottomleft);

// Returns area of (admittedly non-spherical) triangle
// Must write function for spherical area
float triangle_area(Point<float> p1, Point<float> p2, Point<float> p3);

// Returns list of points composing bounding box in clockwise order
Line<Point<float> >* get_box_outline(Point<float> topright, Point<float> bottomleft);

class WPoint {
 public:
  WPoint() { p = 0; w = 0; }
  WPoint(float p, float w) { this->p = p; this->w = w; }
  WPoint(const WPoint& pt) { p = pt.p; w = pt.w; }
  bool operator < (const WPoint& pt) { return p < pt.p; }
  bool operator > (const WPoint& pt) { return p > pt.p; }
  bool operator <= (const WPoint& pt) { return p <= pt.p; }
  bool operator >= (const WPoint& pt) { return p >= pt.p; }
  bool operator == (const WPoint& pt) { return p == pt.p; }
  float p, w;
};

#endif
