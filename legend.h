#include <math.h>
#ifndef __GENIMAGE_LEGEND_H
#define __GENIMAGE_LEGEND_H

class legend {
 public:
  legend(double min, double max, int legend_no = 0, int reversed = 0);
  inline int numcolours() { return num_colours; }
  inline int lookup(int idx) {
    if(idx < 0) {
      return colours[0];
    }
    if(idx >= num_colours) {
      return colours[num_colours - 1];
    }
    return colours[idx]; 
  }
  inline int lookup(double in) {
    if(in < min) {
      return colours[0];
    }
    if(in > max) {
      return colours[num_colours - 1];
    }
    return colours[(int)roundf((in - min) * scale_factor)]; 
  }
  inline void setmax(double newmax) { max = newmax; adjust_range(); adjust_scale_factor(); }
  inline void setmin(double newmin) { min = newmin; adjust_range(); adjust_scale_factor(); }
 private:
  void adjust_scale_factor() { scale_factor = ((num_colours - 1) / range); }
  void adjust_range() { range = max - min; }
  const int* colours;
  int num_colours;
  double scale_factor;
  double min;
  double max;
  double range;
};

#endif
