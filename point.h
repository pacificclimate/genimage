#ifndef __POINT_H
#define __POINT_H
#include <iostream>
#include <math.h>

#ifndef INFINITY
#define INFINITY infinityf()
#endif

using namespace std;

class Point {
 public:
  // Constructors
  Point() { x = 0; y = 0; selected = false; l_used = r_used = false; }
  Point(const double x, const double y, bool sel = false) { this->x = x; this->y = y; selected = sel; }
  Point(const Point& in) { x = in.x; y = in.y; selected = in.selected; }

  // Miscellaneous
  double dist_from(const Point& p) { return hypot(x - p.x, y - p.y); }
  double cross(const Point& p) { double temp; temp = (x * p.y) - (p.x * y); return temp; }

  // Operators
  Point operator - (const Point& p) const { Point temp; temp.x = x - p.x; temp.y = y - p.y; return temp; }
  Point operator + (const Point& p) const { Point temp; temp.x = x + p.x; temp.y = y + p.y; return temp; }
  double operator * (const Point& p) const { return x*p.x +  y*p.y; }
  Point operator * (const double p) const { return Point(x*p, y*p); }
  Point operator / (const Point& p) const { return Point(x/p.x, y/p.y); }
  Point operator / (const double p) const { return Point(x/p, y/p); }
  Point& operator -= (const Point& p) { x -= p.x; y -= p.y; return *this; }
  Point& operator += (const Point& p) { x += p.x; y += p.y; return *this; }
  bool operator == (const Point& p) const { if(p.x == x && p.y == y) { return true; } else { return false; } }
  bool operator != (const Point& p) const { if(p.x != x || p.y != y) { return true; } else { return false; } }
  bool operator <= (const Point& p) const { if((p.x*p.x + p.y*p.y) >= (x*x + y*y)) { return true; } else { return false; } }
  bool operator < (const Point& p) const { if((p.x*p.x + p.y*p.y) > (x*x + y*y)) { return true; } else { return false; } }
  bool operator >= (const Point& p) const { if((p.x*p.x + p.y*p.y) <= (x*x + y*y)) { return true; } else { return false; } }
  bool operator > (const Point& p) const { if((p.x*p.x + p.y*p.y) < (x*x + y*y)) { return true; } else { return false; } }

  friend ostream& operator<<(ostream& os, const Point& f);
  double x, y;
  bool l_used; // Left edge used
  bool r_used; // Right edge used
  bool selected;
};

#endif
