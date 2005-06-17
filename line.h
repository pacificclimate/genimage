#ifndef __LINE_H
#define __LINE_H
#include <iostream>
#include "point.h"
#include "support.h"

using namespace std;

template <class Q>
class Line {
public:
  Q from, to;
  int inside;
  float slope;
  bool used;
  friend ostream& operator<< <> (ostream& os, const Line<Q>& f);

  Line() { used = false; };

  Line(const Q& ifrom, const Q& ito, bool calc_slope = false, bool used = false) { 
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

template <class Q>
ostream& 
operator<<(ostream& os, const Line<Q>& l) {
  os << "[" << l.from << " - " << l.to << "]";
  return os;
}

#endif
