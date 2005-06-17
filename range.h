#ifndef __GENIMAGE_RANGE_H
#define __GENIMAGE_RANGE_H

class range {
 public:
  range() { _modified = 0; }
  float max() { return _max; }
  float min() { return _min; }
  void round(float factor) {
    _min = factor * floorf(_min / factor);
    _max = factor * ceilf(_max / factor);
  }
  void add(float in) { 
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
  float _max;
  float _min;
};

#endif
