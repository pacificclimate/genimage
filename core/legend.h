#include "legends.h"
#include "range.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef __GENIMAGE_LEGEND_H
#define __GENIMAGE_LEGEND_H

#include <iostream>
#include <algorithm>

class Legend {
public:
  Legend(const Range& r, int legend_no = CONTINUOUS, int reversed = NORMAL): range(r), r_min(r.min()), r_range(r.range()) {
    switch (legend_no) {
    case CONTINUOUS:
      colours = continuous[reversed];
      num_colours = sizeof(continuous[reversed]) / sizeof(int);
      break;
    case STEPWISE:
      colours = stepwise[reversed];
      num_colours = sizeof(stepwise[reversed]) / sizeof(int);
      break;
    default:
      fprintf(stderr, "Invalid legend number!\n");
      exit(2);
      break;
    }
  }
  inline int numcolours() const { return num_colours; }
  inline int lookup(int idx) const {
    if (idx < 0) {
      return colours[0];
    }
    if (idx >= num_colours) {
      return colours[num_colours - 1];
    }
    return colours[idx];
  }
  inline int lookup(double in) const {
    int snap = (int)(((in - r_min) / r_range) * (num_colours - 1) + 0.5);
    int colour = MAX(0, MIN(num_colours - 1, snap));
    assert(colour >= 0 && colour < num_colours);
    return colours[colour];
  }

  const Range range;

private:
  const int* colours;
  int num_colours;
  const double r_min;
  const double r_range;
};

#endif
