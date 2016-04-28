#include "interpdatamanager.h"
#include "support.h"  // for find_in_range
#include <assert.h>
#include <math.h>
#include <stdlib.h> // for rand()

using namespace std;

bool InterpRect::isvalid() const {
  if (src_minlong >= src_maxlong || src_minlat >= src_maxlat)
    return false;
  if (dst_boxlong <= 0 || dst_boxlat <= 0)
    return false;
  if (src_maxlong - src_minlong < dst_boxlong)
    return false;
  if (src_maxlat - src_minlat < dst_boxlat)
    return false;

  return true;
}

int InterpRect::y_size() const {
  return (int)floor((src_maxlat - src_minlat) / dst_boxlat);
}

int InterpRect::x_size() const {
  return (int)floor((src_maxlong - src_minlong) / dst_boxlong);
}

double* InterpRect::get_xgrid() const {
  assert(isvalid());
  const int count = x_size() * 2;
  double* vals = new double[count];
  for (int i = 0; i < count; i += 2) {
    vals[i] = src_minlong + ((float)i - 0.5) * dst_boxlong;
    vals[i + 1] = src_minlong + ((float)i + 0.5) * dst_boxlong;
  }
  return vals;
}

double* InterpRect::get_ygrid() const {
  assert(isvalid());
  const int count = y_size() * 2;
  double* vals = new double[count];
  for (int i = 0; i < count; i += 2) {
    vals[i] = src_maxlat - ((float)i - 0.5) * dst_boxlat;
    vals[i + 1] = src_maxlat - ((float)i + 0.5) * dst_boxlat;
  }
  return vals;
}
