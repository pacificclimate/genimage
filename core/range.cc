#include "range.h"
#include <math.h>
#include <algorithm>
#include <list>
#include <utility>

Range::Range(const double* values, const int num_values, RangeType type, double missing) {
  this->type = type;
  this->missing = missing;
  setRange(values, num_values);
}
Range::Range(const double* values, const int rows, const int cols, RangeType type, double missing) {
  this->type = type;
  this->missing = missing;
  setRange(values, rows, cols);
}
Range::Range(const double min, const double max, RangeType type, double missing) {
  _max = max;
  _min = min;
  this->type = type;
  this->missing = missing;
}
void Range::round(const double factor) {
  _min = factor * floor(_min / factor);
  _max = factor * ceil(_max / factor);
}
void Range::add(const double in) { 
  if(in < missing) {
    if(in > _max) { _max = in; }; 
    if(in < _min) { _min = in; };
  }
}

void Range::setRange(const double* values, const int rows, const int cols) {
  // FIXME: What if values[0] is missing?
  std::list<std::pair<double, int> > mins;
  std::list<std::pair<double, int> > maxes;
  int offsets[5][2] = { {0, 0}, {0, 99}, {99, 0}, {99, 99}, {50, 50} };
  for(int i = 0; i < rows; i += 100) {
    for(int j = 0; j < cols; j += 100) {
      for(int k = 0; k < 5; k++) {
	int off_k = MIN(j + offsets[k][0], cols) + (MIN(i + offsets[k][1], rows) * cols);
	const double val = values[off_k];
	if(val != missing) {
	  if(mins.empty()) {
	    mins.push_front(std::pair<double, int>(val, (i * cols) + j));
	    maxes.push_front(std::pair<double, int>(val, (i * cols) + j));
	    break;
	  } else if(mins.front().first > val) {
	    mins.push_front(std::pair<double, int>(val, (i * cols) + j));
	    if(mins.size() > 5)
	      mins.pop_back();
	    break;
	  } else if(maxes.front().first < val) {
	    maxes.push_front(std::pair<double, int>(val, (i * cols) + j));
	    if(maxes.size() > 5)
	      maxes.pop_back();
	    break;
	  }
	}
      }
    }
  }


  const int num_values = rows * cols;
  double minv, maxv;
  maxv = -INFINITY;
  minv = INFINITY;
  if(!num_values) {
    return;
  }
  const double* vals = values + num_values;
  do {
    vals--;
    const double val = *vals;
    if(val != missing) {
     if(val > maxv) maxv = val;
      else if(val < minv) minv = val;
    }
  } while(vals != values);
  _max = maxv;
  _min = minv;
}

// This is not fast but I see no way to optimize it further without making assumptions about the data
void Range::setRange(const double* values, const int num_values) {
  double minv, maxv;
  maxv = -INFINITY;
  minv = INFINITY;
  if(!num_values) {
    return;
  }

  const double* vals = values + num_values;
  do {
    vals--;
    const double val = *vals;
    if(val != missing) {
      if(val > maxv) maxv = val;
      if(val < minv) minv = val;
    }
  } while(vals != values);
  _max = maxv;
  _min = minv;
}
