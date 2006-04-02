#include "legends.h"
#include "range.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef __GENIMAGE_LEGEND_H
#define __GENIMAGE_LEGEND_H

class Legend {
 public:
  Legend(Range r, int legend_no = CONTINUOUS, int reversed = NORMAL) {
    range = r;
    switch(legend_no) {
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
    if(idx < 0) {
      return colours[0];
    }
    if(idx >= num_colours) {
      return colours[num_colours - 1];
    }
    return colours[idx]; 
  }
  inline int lookup(double in) const {
    int colour;
    if(in < range.min()) {
      colour = 0;
    } else if(in > range.max()) {
      colour = num_colours - 1;
    } else {
      colour = (int)roundf(((in - range.min()) / range.range()) * (num_colours - 1));
    }
    return colours[colour];
  }

  Range range;

 private:
  const int* colours;
  int num_colours;
};

#endif
