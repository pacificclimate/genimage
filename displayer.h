#ifndef __GENIMAGE_DISPLAYER_H
#define __GENIMAGE_DISPLAYER_H

#include <string>
#include "canvas.h"
#include "legend.h"
#include "point.h"
#include "canvas.h"

using namespace std;

#define LEG_BOTTOM_TEXT_HEIGHT 55
#define LEG_LABELS_TEXT_HEIGHT 20
#define LEG_TOP_TEXT_HEIGHT 20
#define LEG_TOP_SPACE_HEIGHT 35
#define LEG_HEIGHT 120
#define LAT_WIDTH 40
#define LAT_EXTRA_HEIGHT 10
#define LEG_EXTRA_RWIDTH 20
#define BORDER_WIDTH 1

#define POINT_SIZE 4
#define LINE_WIDTH 3

class Displayer {
public:
  Displayer() {
    colour_map = 0;
    colour_map_rev = 0;
    lon_text_spacing = 20;
    lat_text_spacing = 10;
    leg_dec_places = 1;
    range_dynamic = false;
    grid = false;
    c = 0;
  }
  ~Displayer() {
    if(c)
      delete c;
  }
  int grid, range_dynamic;
  int colour_map, colour_map_rev;
  int img_width, img_height;
  int lon_text_spacing, lat_text_spacing;
  int leg_dec_places;
  double range_min, range_max;
  string leg_text;
  string credit_text;
  string identify_text;

  int map_width, map_height;
  int leg_width, leg_height;
  int lat_width, lat_height;

  int map_offset_x, map_offset_y;
  int leg_offset_x, leg_offset_y;
  int lat_offset_x, lat_offset_y;

  // Set offsets and such
  void setOffsets(const gdImagePtr basemap);

  // Set ranges
  void setRanges(const Range& datarange);

  // Fill in gaps in the map
  void fillGaps();

  // Fill in gaps in the region map
  void fillGapsRegionMap();

  Legend* getLegend(const Range& datarange);

  // Get a canvas to draw on
  void createCanvas(string fontfile);

  // Draw a scale on the map indicating the range of values displayed (for maps)
  void drawScale(Legend& leg_colours);

  // Draw a legend on the map mapping symbols to names (for plots)
  void drawLegend();
  
  // Draw (and possibly label) tickmarks on the map/plot
  void drawTicks(const Range& xrange, const Range& yrange);

  // Draw a map
  void drawMap(const gdImagePtr basemap, const Legend& leg_colours, const int* draw_mask, const int* data_mask, const double* data, const double* grid_lats, const double* grid_longs, const int rows, const int cols);

  // Copy basemap to map
  void copyMap(const gdImagePtr basemap);

  // Draw the polygon on the map
  void drawPolygon(int numpoints, Point** points, Range xrange, Range yrange);

  // Draw the decal on the map
  void plotDecal(string decalfile, Point dpoint, Range xrange, Range yrange);
    
  // Draw identity text on map
  void drawIdentifyText();

  // Write out PNG
  bool writePng(string filename);

private:
  Canvas* c;
};

#endif
