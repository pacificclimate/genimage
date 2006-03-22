#include "point.h"

ostream& operator<<(ostream& os, const Point& f) {
  os << "(" << f.x << "," << f.y << ")";
  return os;
}
