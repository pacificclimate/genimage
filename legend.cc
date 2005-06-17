#include "legend.h"
#include "legends.h"

legend::legend(float min, float max, int legend_no, int reversed) {
  this->min = min;
  this->max = max;
  if(legend_no == 0) {
    // "Continuous" colourmap
    if(reversed) {
      colours = continuous_rev;
    } else {
      colours = continuous;
    }
    num_colours = 50;
  } else if(legend_no == 1) {
    if(reversed) {
      // "Stepwise" colourmap
      colours = stepwise_rev;
    } else {
      colours = stepwise;
    }
    num_colours = 12;
  }
  adjust_range();
  adjust_scale_factor();
}
