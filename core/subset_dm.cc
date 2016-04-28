#include "subset_dm.h"
#include "genimage.h"
#include "line.h"
#include <iostream>
#include <algorithm>

using namespace std;

enum OCEAN_DRAW_TYPE {DRAW_LAND = 0, DRAW_ALL, DRAW_OCEAN};

Window SubsetDataManager::get_data_window() {
  double size_x_2 = (config.lowerright.x - config.upperleft.x) / 2;
  double size_y_2 = (config.upperleft.y - config.lowerright.y) / 2;
  Point ctr(clip_to_range(config.upperleft.x + size_x_2, config.lowerright.x - size_x_2, config.center.x),
            clip_to_range(config.lowerright.y + size_y_2, config.upperleft.y - size_y_2, config.center.y));
  return Window(ctr.y + size_y_2, ctr.x - size_x_2, ctr.y - size_y_2, ctr.x + size_x_2);
}

// Affects the input data
void SubsetDataManager::shift_longs(double* longs, const int size) {
  for (int i = 0; i < size; i++)
    longs[i] = fmod(longs[i] - config.center.x + 540, 360) + config.center.x - 180;
  if (longs[0] == longs[size - 1])
    longs[0] -= 360;
}

int SubsetDataManager::sub_y_size(const list<Window>& l) {
  const Window& w = *(l.begin());
  return (int)(w.bottom - w.top) + 1;
}

int SubsetDataManager::sub_x_size(const list<Window>& l) {
  int sum = 0;
  for (list<Window>::const_iterator i = l.begin(); i != l.end(); i++)
    sum += (int)((*i).right - (*i).left) + 1;
  return sum;
}
