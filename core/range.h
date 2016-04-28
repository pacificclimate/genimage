#ifndef __GENIMAGE_RANGE_H
#define __GENIMAGE_RANGE_H

#ifndef _ISOC99_SOURCE
#define _ISOC99_SOURCE
#endif

#include <math.h>


#ifndef MAX
#define MAX(a,b)   (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b)   (((a) < (b)) ? (a) : (b))
#endif


class Range {
public:
  enum RangeType { RANGE_LAT, RANGE_LON, RANGE_DATA };

  Range(const Range& r) { _min = r.min(); _max = r.max(); missing = r.missing; type = r.type; }
  Range(const double* values, const int num_values, RangeType type = RANGE_DATA, double missing = 1e20);
  Range(const double* values, const int rows, const int cols, RangeType type = RANGE_DATA, double missing = 1e20);
  Range(const double min, const double max, RangeType type = RANGE_DATA, double missing = 1e20);
  Range(RangeType type = RANGE_DATA, double missing = 1e20) { this->type = type; _min = INFINITY; _max = -INFINITY; this->missing = missing; }
  void setRange(const double* values, const int num_values);
  void setRange(const double* values, const int rows, const int cols);
  double max() const { return _max; }
  double min() const { return _min; }
  void setmax(double max) { _max = max; }
  void setmin(double min) { _min = min; }
  double range() const { return _max - _min; }
  void round(double factor);
  void add(double in);
private:
  double _max;
  double _min;
  RangeType type;
  double missing;
};

class LatLonRange {
public:
  LatLonRange() {
    min = INFINITY;
    max = -INFINITY;
    minlat = 0;
    minlong = 0;
    maxlat = 0;
    maxlong = 0;
  }
  void addPoint(const double data, const double newlat, const double newlong) {
    if (data > max) {
      max = data;
      maxlong = newlong;
      maxlat = newlat;
    }
    if (data < min) {
      min = data;
      minlong = newlong;
      minlat = newlat;
    }
  }
  double getMax() const {
    return max;
  }
  double getMaxLat() const {
    return maxlat;
  }
  double getMaxLong() const {
    return maxlong;
  }
  double getMin() const {
    return min;
  }
  double getMinLat() const {
    return minlat;
  }
  double getMinLong() const {
    return minlong;
  }

private:
  double min, minlat, minlong;
  double max, maxlat, maxlong;
};


#endif
