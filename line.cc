#include "support.h"
#include "line.h"

ostream& operator<<(ostream& os, const Line& l) {
  os << "[" << l.from << " - " << l.to << "]";
  return os;
}

Point* ls_intersect(const Line& l1, const Line& l2, bool ignore_first_range) {
  return ls_intersect(l1.from, l1.to, l2.from, l2.to, ignore_first_range);
}

Point* ls_intersect(const Point& p00, const Point& p01, const Point& p10, const Point& p11, bool ignore_first_range) {
  float mA, mB, x, y, s, t;
  int vlA = 0, vlB = 0;
  Point* point = 0;

  // Check for vertical lines

  if((p00.x - p01.x) == 0) {
    // A is vertical line
    vlA = 1;
    mA = INFINITY;
  } else {
    mA = (p00.y - p01.y) / (p00.x - p01.x);
  }

  if((p10.x - p11.x) == 0) {
    // B is vertical line
    vlB = 1;
    mB = INFINITY;
  } else {
    mB = (p10.y - p11.y) / (p10.x - p11.x);
  }

  //mA = clip(mA);
  //mB = clip(mB);

  // If parallel or antiparallel...
  if(mA == mB) {
    return 0;
  }

  if(vlA) {
    //cout << "(" << p00.x << " - " << p11.x << ") = " << (p00.x - p11.x) << endl;

    // Line A (p00 - p01) is a vertical line
    x = p00.x;
    y = ((p00.x - p11.x) * mB + p11.y);
    s = ((y - p00.y) / (p01.y - p00.y));
    t = ((x - p10.x) / (p11.x - p10.x));
  } else if(vlB) {
    // Line B (p10 - p11) is a vertical line
    x = p10.x;
    y = ((p10.x - p01.x) * mA + p01.y);
    s = ((x - p00.x) / (p01.x - p00.x));
    t = ((y - p10.y) / (p11.y - p10.y));
  } else {
    // Neither line is a vertical line
    x = ( - mB * p10.x + p10.y + mA * p00.x - p00.y ) / ( mA - mB );
    y = (mA * ( x - p00.x ) + p00.y);
    s = ((x - p00.x) / (p01.x - p00.x));
    t = ((x - p10.x) / (p11.x - p10.x));
  }

  // Check range
  if((ignore_first_range || (s >= 0 && s <= 1)) && t >= 0 && t <= 1) {
    point = new Point(x, y);
  }

  // Return whatever we got, or 0 (NULL) if no intersect
  return point;
}

// Returns list of points composing bounding box in clockwise order
Line* get_box_outline(Point topright, Point bottomleft) {
  Line* lines = new Line[QUAD_SIDES];
  
  lines[0] = Line(Point(bottomleft.x, topright.y), Point(topright.x, topright.y), true);
  lines[1] = Line(Point(topright.x, topright.y), Point(topright.x, bottomleft.y), true);
  lines[2] = Line(Point(topright.x, bottomleft.y), Point(bottomleft.x, bottomleft.y), true);
  lines[3] = Line(Point(bottomleft.x, bottomleft.y), Point(bottomleft.x, topright.y), true);
  
  return lines;
}

