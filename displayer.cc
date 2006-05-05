#include <iostream>
#include <stdio.h>
#include "genimage.h"
#include "displayer.h"
#include "legend.h"
#include "range.h"
#include "legends.h"

// Set offsets and such
void Displayer::setOffsets(const gdImagePtr basemap) {
  plot_width = basemap->sx + 2 * BORDER_WIDTH;
  plot_height = basemap->sy + 2 * BORDER_WIDTH;

  plot_offset_x = YAXIS_WIDTH;
  plot_offset_y = YAXIS_EXTRA_HEIGHT;

  yaxis_width = YAXIS_WIDTH;
  yaxis_height = plot_height + YAXIS_EXTRA_HEIGHT;

  yaxis_offset_x = 0;
  yaxis_offset_y = 0;

  xaxis_width = plot_width + XAXIS_EXTRA_WIDTH + YAXIS_WIDTH;
  xaxis_height = XAXIS_HEIGHT;

  xaxis_offset_x = 0;
  xaxis_offset_y = yaxis_height + yaxis_offset_y;

  leg_width = xaxis_width;
  leg_height = MAP_LEG_HEIGHT;

  leg_offset_x = 0;
  leg_offset_y = xaxis_offset_y + xaxis_height;

  identify_width = leg_width;
  identify_height = IDENTIFY_HEIGHT;

  identify_offset_x = 0;
  identify_offset_y = leg_offset_y + leg_height;

  img_width = leg_width;
  img_height = leg_height + yaxis_height + xaxis_height + identify_height;
}

// Set offsets and such
void Displayer::setScatterOffsets() {
  latlon_plot = false;

  plot_width = SCATTER_WIDTH + 2 * BORDER_WIDTH;
  plot_height = SCATTER_HEIGHT + 2 * BORDER_WIDTH;

  plot_offset_x = YAXIS_WIDTH + YAXIS_TITLE_WIDTH;
  plot_offset_y = YAXIS_EXTRA_HEIGHT;

  yaxis_width = YAXIS_WIDTH + YAXIS_TITLE_WIDTH;
  yaxis_height = plot_height + YAXIS_EXTRA_HEIGHT;

  yaxis_offset_x = 0;
  yaxis_offset_y = 0;

  leg_width = SCATTER_LEG_WIDTH + 2 * BORDER_WIDTH;
  leg_height = plot_height;

  leg_offset_x = plot_offset_x + plot_width + XAXIS_EXTRA_WIDTH;
  leg_offset_y = plot_offset_y;

  xaxis_width = yaxis_width + plot_width + leg_width + XAXIS_EXTRA_WIDTH;
  xaxis_height = XAXIS_HEIGHT + XAXIS_TITLE_HEIGHT;

  xaxis_offset_x = 0;
  xaxis_offset_y = yaxis_height + yaxis_offset_y;

  identify_width = xaxis_width;
  identify_height = IDENTIFY_HEIGHT;

  img_width = xaxis_width;
  img_height = yaxis_height + xaxis_height + identify_height;

  identify_offset_x = 0;
  identify_offset_y = xaxis_offset_y + xaxis_height;
}

// Copy basemap to map
void Displayer::copyMap(const gdImagePtr basemap) {
  c->copy(basemap, plot_offset_x, plot_offset_y, 0, 0, basemap->sx, basemap->sy);
}

#define LINE_HEIGHT 11

// Draw a legend on the map mapping symbols to names (for plots)
void Displayer::drawLegend(list<LegendToken* >& vars) {
  const int true_leg_height = leg_height - 2 * BORDER_WIDTH;
  const int true_leg_width = leg_width - 2 * BORDER_WIDTH;

  const int true_leg_offset_x = leg_offset_x + BORDER_WIDTH;
  const int true_leg_offset_y = leg_offset_y + BORDER_WIDTH;

  c->colour = 0x00FFFFFF;
  c->setAlpha(0);
  c->fillRect(leg_offset_x, yaxis_offset_y, leg_width, leg_offset_y - yaxis_offset_y);
  c->fillRect(true_leg_offset_x, true_leg_offset_y, true_leg_width, true_leg_height);
  c->setAlpha(1);

  list<LegendToken* >::iterator liter = vars.begin();

  c->fontsize = 8;
  int x = true_leg_offset_x + 20;
  int y = true_leg_offset_y + 1;
  for(; liter != vars.end(); liter++) {
    c->colour = (*liter)->colour;
    if((*liter)->filled) {
      c->fillSymbol((*liter)->sym, x - 15, y);
    } else {
      c->drawSymbol((*liter)->sym, x - 15, y);
    }
    c->colour = 0x00000000;
    c->drawText((*liter)->name, x, y, Canvas::TOP, Canvas::LEFT);
    y += LINE_HEIGHT;
  }
}

void Displayer::fillMapGaps() {
  c->setAlpha(0);
  c->fillRectAbs(plot_offset_x, yaxis_offset_y, 
		 plot_offset_x + plot_width, plot_offset_y - 1, 0x00FFFFFF);
  c->fillRect(plot_offset_x + plot_width, yaxis_offset_y, 
	      XAXIS_EXTRA_WIDTH, yaxis_height, 0x00FFFFFF);
  c->setAlpha(1);
}

void Displayer::fillGapsRegionMap() {
  c->setAlpha(0);
  c->fillRect(leg_offset_x, leg_offset_y, leg_width, leg_height, 0x00FFFFFF);
  c->setAlpha(1);
}

// Create canvas
void Displayer::createCanvas(string fontfile) {
  c = new Canvas(fontfile, img_width, img_height);
}

int dec_places_needed(double value, int dec_places) {
  while(value < 1) {
    dec_places++;
    value *= 10;
  }
  return dec_places;
}

// Draw grid for scatter plot
void Displayer::drawScatterGrid(const Range& xrange, const Range& yrange) {
  // Lon grid
  double factor = (plot_width - 1) / xrange.range();
  double min_x = xrange.min();
  double max_x = xrange.max();
  double min_y = yrange.min();
  double max_y = yrange.max();

  c->colour = gdStyled;
  
  c->setStyle(DASHED);
  for(double i = floor(max_x / x_text_spacing) * x_text_spacing; i >= min_x; i -= x_text_spacing) {
    const int x = plot_offset_x + (int)round((i - min_x) * factor);
    const int y = plot_offset_y;

    // Draw vertical lines
    if(i > -1E-12 && i < 1E-12) {
      c->drawLine(x, y, x, y + plot_height - 1, 0x00808080);
    } else {
      c->drawLine(x, y, x, y + plot_height - 1);
    }
  }
    
  // Lat points
  factor = (plot_height - 1) / yrange.range();

  for(double i = floor(max_y / y_text_spacing) * y_text_spacing; i >= min_y; i -= y_text_spacing) {
    const int x = plot_offset_x;
    const int y = plot_offset_y + (int)round((max_y - i) * factor);

    // Draw horizontal lines
    if(i > -1E-12 && i < 1E-12) {
      c->drawLine(x, y, x + plot_width - 1, y, 0x00808080);
    } else {
      c->drawLine(x, y, x + plot_width - 1, y);
    }
  }
}

// Draw the titles on the axes
void Displayer::drawAxisTitles() {
  c->fontsize = 12;
  c->drawText(xaxis_text, plot_offset_x + ((plot_width - 1) / 2), xaxis_offset_y + xaxis_height - 1, Canvas::BOTTOM, Canvas::CENTER);
  c->drawText(yaxis_text, yaxis_offset_x + 5, plot_offset_y + ((plot_width - 1) / 2), Canvas::MIDDLE, Canvas::LEFT, 0.5 * M_PI);
}
  
// Draw (and possibly label) tickmarks on the map/plot
void Displayer::drawTicks(const Range& xrange, const Range& yrange) {
  // Lon grid
  double factor = (plot_width - 1) / xrange.range();
  double min_x = xrange.min();
  double max_x = xrange.max();
  double min_y = yrange.min();
  double max_y = yrange.max();
  char format[32];

  c->setAlpha(0);
  c->fillRect(xaxis_offset_x, xaxis_offset_y, xaxis_width, xaxis_height, 0x00FFFFFF);
  c->setAlpha(1);
  c->colour = 0x00000000;
  c->fontsize = 10;

  // Set the legend's decimal places
  sprintf(format, "%%0.%if%%s", dec_places_needed(x_text_spacing, 0));

  for(double i = floor(max_x / x_text_spacing) * x_text_spacing; i >= min_x; i -= x_text_spacing) {
    const int x = plot_offset_x + (int)round((i - min_x) * factor);
    const int y = xaxis_offset_y;

    // Put bottom tick marks on
    c->drawLine(x, y, x, y + 3);
    
    // Format the lon string nicely
    char output[32];

    if(!latlon_plot || i == 0) {
      sprintf(output, format, i, "");
    } else if(i < 0) {
      sprintf(output, format, i, "W");
    } else if(i > 0) {
      sprintf(output, format, i, "E");
    }
    
    c->drawText(output, x, y + 4, Canvas::TOP, Canvas::CENTER);
  }

  c->setAlpha(0);
  c->fillRect(yaxis_offset_x, yaxis_offset_y, yaxis_width, yaxis_height, 0x00FFFFFF);
  c->setAlpha(1);
    
  // Lat points
  // Set the legend's decimal places
  sprintf(format, "%%0.%if%%s", dec_places_needed(y_text_spacing, 0));

  factor = (plot_height - 1) / yrange.range();

  for(double i = floor(max_y / y_text_spacing) * y_text_spacing; i >= min_y; i -= y_text_spacing) {
    const int x = yaxis_offset_x;
    const int y = plot_offset_y + (int)round((max_y - i) * factor);

    c->drawLine(yaxis_width - 4, y, yaxis_width - 1, y);

    // Format the lon string nicely
    char output[128];
    if(!latlon_plot || i == 0) {
      sprintf(output, format, i, "");
    } else if(i < 0) {
      sprintf(output, format, fabs(i), "S");
    } else if(i > 0) {
      sprintf(output, format, i, "N");
    }

    c->drawText(output, x + yaxis_width - 4, y, Canvas::MIDDLE, Canvas::RIGHT);
  }
}

// Clear plot
void Displayer::clearPlot() {
  c->colour = 0x00FFFFFF;
  c->setAlpha(0);
  c->fillRect(plot_offset_x, plot_offset_y, plot_width, plot_height);
  c->setAlpha(1);
}

// Draw scatter plot
void Displayer::drawScatter(std::list<ScatterVars*>& vars, const Range& xrange, const Range& yrange) {
  list<ScatterVars*>::iterator i = vars.begin();
  double xfactor = (double)(plot_width) / xrange.range();
  double yfactor = (double)(plot_height) / yrange.range();

  for(; i != vars.end(); i++) {
    int x = (int)round(((*i)->datx - xrange.min()) * xfactor) + plot_offset_x;
    int y = (int)round((yrange.max() - (*i)->daty) * yfactor) + plot_offset_y;

    c->colour = (*i)->symbol->colour;
    if((*i)->symbol->filled) {
      c->fillSymbol((*i)->symbol->sym, x - 5, y - 5);
    } else {
      c->drawSymbol((*i)->symbol->sym, x - 5, y - 5);
    }
  }

  // Clean up the rect around the map
  c->drawRect(plot_offset_x, plot_offset_y, plot_width, plot_height, 0x00000000);
}

// Draw lines for timeslice scatter plot
void Displayer::drawLines(std::list<ScatterVars*>& vars, const Range& xrange, const Range& yrange) {
  list<ScatterVars*>::iterator i = vars.begin();
  ScatterVars* old = 0;
  double xfactor = (double)(plot_width) / xrange.range();
  double yfactor = (double)(plot_height) / yrange.range();
  int old_x = 0, old_y = 0;
  for(; i != vars.end(); i++) {
    int x = (int)round(((*i)->datx - xrange.min()) * xfactor) + plot_offset_x;
    int y = (int)round((yrange.max() - (*i)->daty) * yfactor) + plot_offset_y;

    if(old && (old->model == (*i)->model && old->expt == (*i)->expt)) {
      c->drawLine(old_x, old_y, x, y, (*i)->symbol->colour);
    }

    old_x = x;
    old_y = y;
    old = *i;
  }
}

// Sets appropriate tick mark distances
void Displayer::setTicks(const Range& xrange, const Range& yrange) {
  x_text_spacing = tick_spacing(xrange, DESIRED_XTICKS);
  y_text_spacing = tick_spacing(yrange, DESIRED_YTICKS);
}

// Draw a map
void Displayer::drawMap(const gdImagePtr basemap, const Legend& leg_colours, const int* draw_mask, const int* data_mask, const double* data, const double* grid_lats, const double* grid_longs, const int rows, const int cols) {
  int* xpoint = new int[cols + 1];
  int* ypoint = new int[rows + 1];
  const int datasize = rows * cols;

  const Range xrange(grid_longs, cols + 1);
  const Range yrange(grid_lats, rows + 1);
  const double minlong = xrange.min();
  const double difflong = xrange.range();
  const double maxlat = yrange.max();
  const double difflat = yrange.range();

  unsigned int* data_colours = new unsigned int[datasize];
  unsigned char** map = basemap->pixels;
  int** img = c->img->tpixels;

  const int true_plot_height = plot_height - 2 * BORDER_WIDTH;
  const int true_plot_width = plot_width - 2 * BORDER_WIDTH;

  const int true_plot_offset_x = plot_offset_x + BORDER_WIDTH;
  const int true_plot_offset_y = plot_offset_y + BORDER_WIDTH;

  // Generate the array of Y axis grid points
  for(int y = 0; y <= rows; y++) {
    ypoint[y] = (int)round(((maxlat - grid_lats[y]) / difflat) * true_plot_height);
  }

  // Generate the array of X axis grid points
  for(int x = 0; x <= cols; x++) {
    xpoint[x] = (int)round(((grid_longs[x] - minlong) / difflong) * true_plot_width);
  }

  // Preprocess the data array with mask file
  for(int i = 0; i < datasize; i++) {
    if(draw_mask[i]) {
      data_colours[i] = leg_colours.lookup(data[i]);
    } else {
      data_colours[i] = 0x00FFFFFF;
    }
  }

  const int* datamask_ptr = data_mask;
  int y = ypoint[0];
  for(int j = 0; j < rows; j++) {
    const int yend = ypoint[j + 1];
    const unsigned int* data_colours_ptr = &data_colours[j * cols];
    for(; y < yend; y++) {
      const unsigned char* mlineptr = map[y];
      int* lineptr = img[y + true_plot_offset_y] + true_plot_offset_x;
      int x = xpoint[0];
      datamask_ptr = &data_mask[j * cols];
      for(int i = 0; i < cols; i++) {
 	// Look up colour
 	const unsigned int colour = data_colours_ptr[i];
 	const int xend = xpoint[i + 1];
 	// Map is either 1 (white) or 0 (black)
 	// Negation of 1 is 0xffffffff; negation of 0 is 0x00000000
 	// AND'ing that with the colour gives either black or the colour
 	if(*datamask_ptr) {
 	  // If the data mask is 1 at the location, then the location is not masked
 	  for(; x < xend; x++) {
 	    // Apply colour to map
 	    lineptr[x] = colour & -((int)mlineptr[x]);
 	  }
 	} else {
 	  // If the data mask is 0 at the location, then the location is masked
 	  for(; x < xend; x++) {
 	    // Apply stipple and colour to map
 	    //lineptr[x] = colour & -((int)mlineptr[x]);
 	    lineptr[x] = colour & -((const int)mlineptr[x] & x);
 	  }
	}
 	datamask_ptr++;
      }
    }
  }

  // After we're done with the map, we put the grid on it
  c->colour = 0x00A0A0A0;
  if(grid) {
    for(int j = 1; j < rows; j++) {
      c->drawLine(true_plot_offset_x, true_plot_offset_y + ypoint[j], true_plot_offset_x + plot_width, true_plot_offset_y + ypoint[j]);
    }
    for(int j = 0; j < cols; j++) {
      c->drawLine(true_plot_offset_x + xpoint[j], true_plot_offset_y, true_plot_offset_x + xpoint[j], true_plot_offset_y + plot_height);
    }
  }

  delete[] data_colours;
  delete[] xpoint;
  delete[] ypoint;
}

void Displayer::drawPolygon(int numpoints, Point** points, Range xrange, Range yrange) {
  c->colour = 0x00000000;

  const int true_plot_height = plot_height - 2 * BORDER_WIDTH;
  const int true_plot_width = plot_width - 2 * BORDER_WIDTH;

  const int true_plot_offset_x = plot_offset_x + BORDER_WIDTH;
  const int true_plot_offset_y = plot_offset_y + BORDER_WIDTH;

  // Lines to go on graph
  // This code is acceptable now
  c->setAntiAliased();
  c->setLineThickness(LINE_WIDTH);
  int x1, x2, y1, y2;
  for(int i = 0; i < numpoints; i++) {
    // Convert to XY from latlon
    x1 = lontox(true_plot_width, 0, xrange.min(), xrange.max(), points[i]->x);
    y1 = lattoy(true_plot_height, 0, yrange.min(), yrange.max(), points[i]->y);
    
    x2 = lontox(true_plot_width, 0, xrange.min(), xrange.max(), points[(i + 1) % numpoints]->x);
    y2 = lattoy(true_plot_height, 0, yrange.min(), yrange.max(), points[(i + 1) % numpoints]->y);
    
    // Lines
    c->drawLine(x1 + true_plot_offset_x, y1 + true_plot_offset_y, x2 + true_plot_offset_x, y2 + true_plot_offset_y);
  }
  c->setLineThickness(1);
  for(int i = 0; i < numpoints; i++) {
    // Convert to XY from latlon
    x1 = lontox(plot_width, 0, xrange.min(), xrange.max(), points[i]->x);
    y1 = lattoy(plot_height, 0, yrange.min(), yrange.max(), points[i]->y);
    
    // Point
    if(points[i]->selected) {
      c->fillRectAbs(x1 - (POINT_SIZE - 2) + true_plot_offset_x, y1 - (POINT_SIZE - 2) + true_plot_offset_y, x1 + (POINT_SIZE - 2) + true_plot_offset_x, y1 + (POINT_SIZE - 2) + true_plot_offset_y, 0x00FF0000);
      
    } else {
      c->fillRectAbs(x1 - 3 + true_plot_offset_x, y1 - 3 + true_plot_offset_y, x1 + 3 + true_plot_offset_x, y1 + 3 + true_plot_offset_y, 0x00000000);
    }
    
    // Box around point
    c->drawRect(x1 - (POINT_SIZE - 1) + true_plot_offset_x, y1 - (POINT_SIZE - 1) + true_plot_offset_y, 2 * (POINT_SIZE - 1), 2 * (POINT_SIZE - 1), 0x00FFFFFF);
    c->drawRect(x1 - POINT_SIZE + true_plot_offset_x, y1 - POINT_SIZE + true_plot_offset_y, 2 * POINT_SIZE, 2 * POINT_SIZE);
  }      
}

void Displayer::plotDecal(string decalfile, Point dpoint, Range xrange, Range yrange) {
  // Display the data point
  // First, load the decal image (erroring if DNE)

  const int true_plot_offset_x = plot_offset_x + BORDER_WIDTH;
  const int true_plot_offset_y = plot_offset_y + BORDER_WIDTH;

  FILE* infile = fopen(decalfile.c_str(), "r");
  if(!infile) {
    exit(3);
  }
  gdImagePtr decalimg = gdImageCreateFromPng(infile);
  if(!decalimg) {
    exit(3);
  }
  fclose(infile);
  
  // Apply the decal image to the right spot
  c->copy(decalimg, 
	  lontox(plot_width, 0, xrange.min(), xrange.max(), dpoint.x) - 
	  ((decalimg->sx - 1) >> 1) + true_plot_offset_x, 
	  lattoy(plot_height, 0, yrange.min(), yrange.max(), dpoint.y) - 
	  ((decalimg->sy - 1) >> 1) + true_plot_offset_y, 
	  0, 0, decalimg->sx, decalimg->sy);
}

Legend* Displayer::getLegend(const Range& datarange) {
  Legend* leg_colours = new Legend(datarange, colour_map, colour_map_rev);
  double leg_range_min = xrange_min;
  double leg_range_max = xrange_max;

  if(range_dynamic) {
    leg_range_min = datarange.min();
    leg_range_max = datarange.max();
  }
  leg_colours->range.setmin(leg_range_min);
  leg_colours->range.setmax(leg_range_max);

  return leg_colours;
}

void Displayer::drawScale(Legend& leg_colours) {
  // Generate the colour map
  int num_leg_segments = 10;
  int leg_dec_places = this->leg_dec_places;
  int x, y;

  if(range_dynamic && colour_map == STEPWISE) {
    num_leg_segments = 12;
  }

  // Set the legend's decimal places
  leg_dec_places = dec_places_needed((leg_colours.range.max() - leg_colours.range.min()) / (num_leg_segments), range_dynamic);

  int leg_left = plot_offset_x + BORDER_WIDTH;
  int leg_right = plot_offset_x + plot_width - BORDER_WIDTH - 1;

  int leg_top = leg_offset_y + BORDER_WIDTH;
  int leg_bottom = leg_offset_y + leg_height - MAP_LEG_BOTTOM_TEXT_HEIGHT - BORDER_WIDTH;

  int width = plot_width - 2 * BORDER_WIDTH + 1;
  int height = leg_height - MAP_LEG_BOTTOM_TEXT_HEIGHT + 1;
  int numcolours = leg_colours.numcolours();

  // White out the white areas
  c->colour = 0x00FFFFFF;
  c->setAlpha(0);
  c->fillRect(leg_offset_x, leg_offset_y, yaxis_width, leg_width);
  c->fillRect(leg_right + BORDER_WIDTH + 1, leg_offset_y, XAXIS_EXTRA_WIDTH + 1, height + 2 * BORDER_WIDTH);
  c->fillRect(leg_offset_x, leg_bottom + 2 * BORDER_WIDTH, leg_width, leg_height - height);
  c->setAlpha(1);
  
  // Draw the title for the legend
  c->colour = 0x00000000;
  c->drawText(xaxis_text, leg_left + width / 2, leg_bottom + MAP_LEG_LABEL_HEIGHT, Canvas::TOP, Canvas::CENTER);
      
  // Draw the legend
  double factor = (double)(width - 1) / (double)numcolours;
  for(int i = numcolours - 1; (i + 1); i--) {
    int colour = leg_colours.lookup(i);
    c->fillRectAbs((int)round(i * factor) + leg_left, leg_top, (int)round((i + 1) * factor) + leg_left, leg_bottom, colour);
  }
      
  double lbl_offset;
  if(range_dynamic) {
    lbl_offset = 0;
    factor = (double)(width - 1) / (num_leg_segments);
  } else {
    lbl_offset = factor;
    factor = (double)(width - (2 * lbl_offset) - 1) / num_leg_segments;
  }
  
  // Create the tick marks in the legend and label them
  double scale_factor = (leg_colours.range.max() - leg_colours.range.min()) / (num_leg_segments);
  char formatbuf[128];
  char output[128];
  sprintf(formatbuf, "%%.%if", leg_dec_places);
  c->colour = 0x00000000;
  for(int i = 0; i <= num_leg_segments; i++) {
    x = leg_left + (int)round(i * factor + lbl_offset);
    y = leg_bottom;
    c->drawLine(x, leg_top, x, leg_top + 4);
    c->drawLine(x, leg_bottom - 4, x, leg_bottom);

    sprintf(output, formatbuf, scale_factor * i + leg_colours.range.min());
    c->drawText(output, x, y + 3, Canvas::TOP, Canvas::CENTER);
  }
}

void Displayer::clearIdentifyArea() {
  c->setAlpha(0);
  c->fillRect(identify_offset_x, identify_offset_y, identify_width, identify_height, 0x00FFFFFF);
  c->setAlpha(1);
}

void Displayer::drawIdentifyText() {
  int x, y;

  c->colour = 0x00000000;
  x = identify_offset_x;
  y = identify_offset_y + identify_height;
  c->fontsize = 8;
  c->drawText(identify_text, x, y, Canvas::BOTTOM, Canvas::LEFT);
}

void Displayer::drawCreditText() {
  int x, y;

  c->colour = 0x00000000;
  x = identify_offset_x + identify_width;
  y = identify_offset_y + identify_height;
  c->fontsize = 10;
  c->drawText(credit_text, x, y, Canvas::BOTTOM, Canvas::RIGHT);
}

bool Displayer::writePng(string filename) {
  FILE* f;
  f = fopen(filename.c_str(), "wb");
  if(!f) {
    cerr << "Could not open output file!\n";
    return false;
  }
  gdImagePng(c->img, f);
  fclose(f);
  return true;
}
