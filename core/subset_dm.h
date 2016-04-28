#ifndef __GENIMAGE_SUBSETDATAMANAGER_H
#define __GENIMAGE_SUBSETDATAMANAGER_H

#include "datamanager.h"

class SubsetDataManager : public DataManager {
public:
  SubsetDataManager(Config& c) : DataManager(c) { }
  virtual DataGrid<int> get_slmask(const DataSpec& s) { return (get_data_subset(s, DataManager::get_slmask(s))); }
  virtual DataGrid<int> get_drawmask(const DataSpec& s) { return (get_data_subset(s, DataManager::get_drawmask(s))); }
  virtual DataGrid<int> get_datamask(const DataSpec& s, const double threshold) { return (get_data_subset(s, DataManager::get_datamask(s, threshold))); }
  virtual DataGrid<double> get_areagrid(const DataSpec& s) { return (get_data_subset(s, DataManager::get_areagrid(s))); }
  virtual DataGrid<double> get_data(const DataSpec& s) { return (get_data_subset(s, DataManager::get_data(s))); }

protected:
  int sub_x_size(const list<Window>& l);
  int sub_y_size(const list<Window>& l);
  Window get_data_window();
  void shift_longs(double* longs, const int size);

  //template <typename T> double* get_x_subset(const DataGrid<T>& indat, const list<Window>& l);
  //template <typename T> double* get_y_subset(const DataGrid<T>& indat, list<Window>& l);
  //template <typename T> list<Window> get_subsets(const DataGrid<T>& g);
  //template <typename T> DataGrid<T> get_data_subset(const DataSpec& s, const DataGrid<T>& indat);


  template <typename T> double* get_x_subset(const DataGrid<T>& indat, const list<Window>& l) {
    double* outdat = new double[sub_x_size(l) * 2];
    fprintf(stderr, "old x_grid size: %i, new x_grid size: %i\n", indat.x_size() * 2, sub_x_size(l) * 2);
    const double* invals = indat.x_grid();
    int j = 0;
    for (list<Window>::const_iterator i = l.begin(); i != l.end(); ++i) {
      const Window& w = *i;
      std::copy(&invals[(int)w.left * 2], &invals[(int)w.right * 2 + 1], &outdat[j * 2]);
      j += (int)w.right - (int)w.left + 1;
    }
    return outdat;
  }

  template <typename T> double* get_y_subset(const DataGrid<T>& indat, list<Window>& l) {
    double* outdat = new double[sub_y_size(l) * 2];
    fprintf(stderr, "old y_grid size: %i, new y_grid size: %i\n", indat.y_size() * 2, sub_y_size(l) * 2);
    const double* invals = indat.y_grid();
    const Window& w = *(l.begin());
    std::copy(&invals[(int)w.top * 2], &invals[(int)w.bottom * 2 + 1], outdat);
    return outdat;
  }

  template <typename T> DataGrid<T> get_data_subset(const DataSpec& s, const DataGrid<T>& indat) {
    open_file(get_filename(s));
    list<Window> l = get_subsets(indat);
    list<Window>::const_iterator i;
    DataGrid<T> g(s);
    int incols = indat.x_size();
    int outcols = sub_x_size(l);
    T* invals = indat.values();
    fprintf(stderr, "sub_x_size: %i, sub_y_size: %i\n", sub_x_size(l), sub_y_size(l));
    T* outvals = new T[sub_x_size(l) * sub_y_size(l)];

    int long_off = 0;
    for (i = l.begin(); i != l.end(); i++) {
      const Window& w = *i;

      for (int i = (int)w.top; i <= (int)w.bottom; i++) {
        const int in_offset = i * incols;
        const int out_offset = (i - (int)w.top) * outcols + long_off;
        std::copy(&invals[in_offset + (int)w.left], &invals[in_offset + (int)w.right + 1], &outvals[out_offset]);
      }
      long_off += (int)w.right - (int)w.left + 1;
    }
    g.set_data(sub_x_size(l), sub_y_size(l), indat.missing(), outvals, get_x_subset(indat, l), get_y_subset(indat, l), indat.proj4_string(), indat.projection());
    return (g);
  }

// Returns a list of windows into the data, shifted to eliminate the need for stitching/shifting later
// Effectively gives tools needed to perform subsetting and shifting
  template <typename T> list<Window> get_subsets(const DataGrid<T>& g) {
    list<Window> lw;
    list<int> l[2];
    list<int> h;
    int i;
    const int x_size = g.xgrid_size();
    const int y_size = g.ygrid_size();
    double* x_grid = new double[x_size];
    const double* y_grid = g.y_grid();
    Window w = get_data_window();
    std::copy(g.x_grid(), &g.x_grid()[x_size], x_grid);

    fprintf(stderr, "ll sizes: %i, %i\n", x_size, y_size);
    if (g.projection() == "equidistant_cylindrical") {
      fprintf(stderr, "sdm::get_subsets: About to shift longs\n");
      shift_longs(x_grid, x_size);

      // Accumulate lists of offsets of gridboxes that are within the window
      // FIXME: Should be replaced by some form of binary search; this is inefficient
      // and this will be done more than once per data file with new architecture
      int j = 0;
      for (i = 0; i < x_size && j < 2; i += 2) {
        if (x_grid[i] > x_grid[i + 1]) {
          j++;
        }
        // FIXME: Does using fmod here cause problems in practice?
        if ((fmod(w.left, 180) <= x_grid[i] && fmod(w.right, 180) > x_grid[i]) ||
            (fmod(w.left, 180) <= x_grid[i + 1] && fmod(w.right, 180) > x_grid[i + 1])) {
          l[j].push_back(i / 2);
        }
      }

      // Find the latitudes within this window
      for (i = 0; i < y_size; i += 2) {
        if ((w.bottom <= y_grid[i] && w.top > y_grid[i]) ||
            (w.bottom <= y_grid[i + 1] && w.top > y_grid[i + 1])) {
          h.push_back(i / 2);
        }
      }

      // Add windows into the data to the list
      for (i = 1; i >= 0; i--) {
        if (l[i].size() > 0) {
          lw.push_back(Window(*(min_element(h.begin(), h.end())), *(min_element(l[i].begin(), l[i].end())),
                              *(max_element(h.begin(), h.end())), *(max_element(l[i].begin(), l[i].end()))));
          fprintf(stderr, "Added window\n");
        }
      }
    } else {
      lw.push_back(Window(0, 0, (y_size / 2) - 1, (x_size / 2) - 1));
    }

    return lw;
  }
};

#endif
