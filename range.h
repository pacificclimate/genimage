#include <math.h>

#ifndef __GENIMAGE_RANGE_H
#define __GENIMAGE_RANGE_H

class range {
 public:
  range() { _modified = 0; }
  double max() { return _max; }
  double min() { return _min; }
  void round(double factor) {
    _min = factor * floorf(_min / factor);
    _max = factor * ceilf(_max / factor);
  }
  void add(double in) { 
    if(!_modified) {
      _max = _min = in;
      _modified = 1;
    } else {
      if(in > _max) { _max = in; }; 
      if(in < _min) { _min = in; };
    }
  }
 private:
  int _modified;
  double _max;
  double _min;
};

#endif
