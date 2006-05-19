#ifndef __GENIMAGE_DISPLAYER_H
#define __GENIMAGE_DISPLAYER_H

#include <string>
#include <list>
#include "canvas.h"
#include "legend.h"
#include "point.h"
#include "canvas.h"
#include "scattervars.h"
#include "datamanager.h"

// Maps
#define MAP_LEG_BOTTOM_TEXT_HEIGHT 30
#define MAP_LEG_LABEL_HEIGHT 15
#define MAP_LEG_HEIGHT 70

// Scatter plots
#define SCATTER_WIDTH 525
#define SCATTER_HEIGHT 525
#define SCATTER_LEG_WIDTH 140

// General vars
#define XAXIS_EXTRA_WIDTH 20
#define XAXIS_HEIGHT 30
#define XAXIS_TITLE_HEIGHT 15

#define YAXIS_EXTRA_HEIGHT 10
#define YAXIS_TITLE_WIDTH 15
#define YAXIS_WIDTH 40

#define BORDER_WIDTH 1
#define IDENTIFY_HEIGHT 20

#define POINT_SIZE 4
#define LINE_WIDTH 3

#define DESIRED_XTICKS 5
#define DESIRED_YTICKS 5

class Displayer {
public:
  Displayer() {
    colour_map = 0;
    colour_map_rev = 0;
    x_text_spacing = 20;
    y_text_spacing = 10;
    range_dynamic = false;
    grid = false;
    c = 0;
    latlon_plot = true;
  }
  ~Displayer() {
    if(c)
      delete c;
  }

  bool latlon_plot;

  int grid, range_dynamic;
  int colour_map, colour_map_rev;
  int img_width, img_height;
  double x_text_spacing, y_text_spacing;
  double xrange_min, xrange_max;
  double yrange_min, yrange_max;
  std::string xaxis_text;
  std::string yaxis_text;
  std::string credit_text;
  std::string identify_text;

  int plot_width, plot_height;
  int leg_width, leg_height;
  int identify_width, identify_height;
  int xaxis_width, xaxis_height;
  int yaxis_width, yaxis_height;

  int plot_offset_x, plot_offset_y;
  int leg_offset_x, leg_offset_y;
  int identify_offset_x, identify_offset_y;
  int xaxis_offset_x, xaxis_offset_y;
  int yaxis_offset_x, yaxis_offset_y;

  // Sets offsets for map plots
  void setOffsets(const gdImagePtr basemap);

  // Sets offsets for scatter plots
  void setScatterOffsets();

  // Fill in gaps in the map
  void fillMapGaps();

  // Fill in gaps in the region map
  void fillGapsRegionMap();

  Legend* getLegend(const Range& datarange);

  // Get a canvas to draw on
  void createCanvas(std::string fontfile);

  // Draw a scale on the map indicating the range of values displayed (for maps)
  void drawScale(Legend& leg_colours);

  // Draw a legend on the map mapping symbols to names (for plots)
  void drawLegend(list<LegendToken* >& vars);

  // Draws the titles on the axes
  void drawAxisTitles();
  
  // Draw (and possibly label) tickmarks on the map/plot
  void drawTicks(const Range& xrange, const Range& yrange);

  // Clears the area to be plotted
  void clearPlot();

  // Draw a map
  void drawMap(const gdImagePtr basemap, const Legend& leg_colours, const int* draw_mask, const int* data_mask, const double* data, const double* grid_lats, const double* grid_longs, const int rows, const int cols);

  // Draw scatter plot
  void drawScatter(list<ScatterVars* >& vars, const Range& xrange, const Range& yrange);

  // Draws grid on scatter plot
  void drawScatterGrid(const Range& xrange, const Range& yrange);

  // Draw lines for scatter plot
  void drawLines(std::list<ScatterVars*>& vars, const Range& xrange, const Range& yrange);
  
  // Set intervals for tickmarks
  void setTicks(const Range& xrange, const Range& yrange);

  // Copy basemap to map
  void copyMap(const gdImagePtr basemap);

  // Draw the polygon on the map
  void drawPolygon(int numpoints, Point** points, Range xrange, Range yrange);

  // Draw the decal on the map
  void plotDecal(std::string decalfile, Point dpoint, Range xrange, Range yrange);

  // Clear (to white) the identify text area
  void clearIdentifyArea();
    
  // Draw identity text on map
  void drawIdentifyText();

  // Draw credit text on map
  void drawCreditText();

  // Write out PNG
  bool writePng(std::string filename);

private:
  Canvas* c;
};

#endif
