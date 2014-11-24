#ifndef __LINE_H
#define __LINE_H
#include <iostream>
#include "point.h"
#include "support.h"

using namespace std;

class Line;

Point ls_intersect(const Line& l1, const Line& l2, bool ignore_first_range = false);

Point ls_intersect(const Point& p00, const Point& p01, const Point& p10, const Point& p11, bool ignore_first_range = false);

Point line_intersect(const Line& l1, const Line& l2);

Point line_intersect(const Point& p00, const Point& p01, const Point& p10, const Point& p11);

// Returns list of points composing bounding box in clockwise order
Line* get_box_outline(Point topright, Point bottomleft);

class Line {
public:
  Point from, to;
  int inside;
  float slope;
  bool used;
  friend ostream& operator<< (ostream& os, const Line& f);

  Line() { used = false; };

  Line(const Point& ifrom, const Point& ito, bool calc_slope = false, bool used = false) { 
    this->used = used;
    from.x = ifrom.x; 
    from.y = ifrom.y; 
    to.x = ito.x; 
    to.y = ito.y; 
    inside = -1;
    if(calc_slope) {
      if(to.x == from.x) {
	slope = INFINITY * sign(to.y - from.y);
      } else {
	slope = (to.y - from.y) / (to.x - from.x);
      }
    }
  }
};

#endif
