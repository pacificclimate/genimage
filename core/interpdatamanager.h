#ifndef __GENIMAGE_INTERPDATAMANAGER_H
#define __GENIMAGE_INTERPDATAMANAGER_H

#include "datamanager.h"

class InterpRect {
 public:
  // Does rectangle include some area, and make sense?
  bool isvalid() const;
  double* get_xgrid() const;
  double* get_ygrid() const;
  int x_size() const;
  int y_size() const;

  // Source window which we wish to interpolate
  double src_minlong, src_minlat, src_maxlong, src_maxlat;
  // Width/height of a single gridbox on grid data will be interpolated to
  double dst_boxlong, dst_boxlat;
};

template <typename T> void interpolate_data(const DataGrid<T>& src, DataGrid<T>& dest) {
  // Get source data
  const int src_xsize = src.x_size(), src_ysize = src.y_size();
  const T* src_data = src.values().get();
  const int dest_xsize = dest.x_size(), dest_ysize = dest.y_size();
  T* dest_data = dest.values().get();
  const double missing = src.missing();
  const double* src_x_center = get_centers_from_grid(src.x_grid().get(), src_xsize);
  const double* src_y_center = get_centers_from_grid(src.y_grid().get(), src_ysize);
  const double* dest_x = get_centers_from_grid(dest.x_grid().get(), dest_xsize);
  const double* dest_y = get_centers_from_grid(dest.y_grid().get(), dest_ysize);
  
  for(int j = 0; j < dest_ysize; ++j) {
    bool found;
    double y, y1, y2;
    int iy1, iy2;
    // Find y1
    y = dest_y[j];
    found = find_in_range(y, src_ysize - 1, src_y_center, iy1, false);
    if(!found) continue;
    y1 = src_y_center[iy1];
    
    // y2
    iy2 = iy1+1;
    if(iy2 >= src_ysize) continue;
    y2 = src_y_center[iy2];

    for(int i = 0; i < dest_xsize; ++i) {
      double x, x1, x2;
      int ix1, ix2;
      
      // Find x1
      x = dest_x[i];
      found = find_in_range(x, src_xsize - 1, src_x_center, ix1, true);
      if(!found) continue;
      x1 = src_x_center[ix1];

      // x2
      ix2 = ix1+1;
      if(ix2 >= src_xsize) continue;
      x2 = src_x_center[ix2];
      
      // Do the linear interpolation
      // http://en.wikipedia.org/wiki/Bilinear_interpolation
      
      #define SRC_VALS(x, y)  (src_data[(y) * src_xsize + (x)])
      const double r11 = SRC_VALS(ix1, iy1);
      const double r21 = SRC_VALS(ix2, iy1);
      const double r12 = SRC_VALS(ix1, iy2);
      const double r22 = SRC_VALS(ix2, iy2);
      if(r11 == missing || r21 == missing || r12 == missing || r22 == missing)
	continue;

      const double r1 = ((x2 - x) / (x2 - x1)) * r11 + ((x - x1) / (x2 - x1)) * r21;
      const double r2 = ((x2 - x) / (x2 - x1)) * r12 + ((x - x1) / (x2 - x1)) * r22;
      const double val = ((y2 - y) / (y2 - y1)) * r1 + ((y - y1) / (y2 - y1)) * r2;
      dest_data[(j * dest_xsize) + i] = val;
    }
  }
  delete[] dest_x;
  delete[] dest_y;
  delete[] src_x_center;
  delete[] src_y_center;
}

template <typename T> DataGrid<T> interpolate_grid(const DataGrid<T>& g, const InterpRect& rect) {
  const int dest_xsize = rect.x_size(), dest_ysize = rect.y_size();
  double* dest_xgrid = rect.get_xgrid();
  double* dest_ygrid = rect.get_ygrid();

  T* dest_data = new T[dest_xsize * dest_ysize];
  std::fill(dest_data, &dest_data[dest_xsize * dest_ysize], g.missing());

  DataGrid<T> dg(g.r);
  dg.set_data(dest_xsize, dest_ysize, g.missing(), dest_data, dest_xgrid, dest_ygrid, g.proj4_string(), g.projection(), g.anomaly(), g.base_period());
  interpolate_data(g, dg);

  return dg;
}

#endif //#define __GENIMAGE_INTERPDATAMANAGER_H
