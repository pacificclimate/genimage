#include <iostream>
#include <stdio.h>
#include "genimage.h"
#include "displayer.h"
#include "legend.h"
#include "range.h"
#include "legends.h"

// Set offsets and such
void Displayer::setOffsets(const gdImagePtr basemap) {
  map_width = basemap->sx;
  map_height = basemap->sy;
  leg_width = map_width + LAT_WIDTH + LEG_EXTRA_RWIDTH;
  leg_height = LEG_HEIGHT;
  lat_width = LAT_WIDTH;
  lat_height = map_height + LAT_EXTRA_HEIGHT;
  img_width = leg_width;
  img_height = leg_height + map_height + LAT_EXTRA_HEIGHT;
  map_offset_x = lat_width + BORDER_WIDTH;
  map_offset_y = LAT_EXTRA_HEIGHT - BORDER_WIDTH;
  leg_offset_x = 0;
  leg_offset_y = map_height + LAT_EXTRA_HEIGHT;
  lat_offset_x = 0;
  lat_offset_y = 0;
}

// Copy basemap to map
void Displayer::copyMap(const gdImagePtr basemap) {
  c->copy(basemap, map_offset_x, map_offset_y, 0, 0, basemap->sx, basemap->sy);
}


// Draw a legend on the map mapping symbols to names (for plots)
void Displayer::drawLegend() {
}

void Displayer::fillGaps() {
  c->fillRect(map_offset_x - BORDER_WIDTH, lat_offset_y, map_width + LEG_EXTRA_RWIDTH, LAT_EXTRA_HEIGHT - 2 * BORDER_WIDTH, 0x00FFFFFF);
  c->fillRect(map_offset_x + map_width + BORDER_WIDTH, map_offset_y - BORDER_WIDTH, LEG_EXTRA_RWIDTH, map_height + 2 * BORDER_WIDTH, 0x00FFFFFF);
}

void Displayer::fillGapsRegionMap() {
  c->fillRect(leg_offset_x, leg_offset_y + LEG_TOP_TEXT_HEIGHT, leg_width, leg_height, 0x00FFFFFF);
}

// Create canvas
void Displayer::createCanvas(string fontfile) {
  c = new Canvas(fontfile, img_width, img_height);
}
  
// Draw (and possibly label) tickmarks on the map/plot
void Displayer::drawTicks(const Range& xrange, const Range& yrange) {
  // Lon grid
  double factor = map_width / xrange.range();
  double minlon = xrange.min();
  double maxlon = xrange.max();
  double minlat = yrange.min();
  double maxlat = yrange.max();
  int leg_left = LAT_WIDTH + BORDER_WIDTH + leg_offset_x;

  c->fillRect(leg_offset_x, leg_offset_y, leg_width, LEG_TOP_TEXT_HEIGHT, 0x00FFFFFF);
  c->colour = 0x00000000;
  c->fontsize = 10;

  for(int i = (int)floorf(maxlon / lon_text_spacing) * lon_text_spacing; i > minlon; i -= lon_text_spacing) {
    int x = leg_left + (int)round((i - minlon) * factor);
    int y = leg_offset_y;

    // Put bottom tick marks on
    c->drawLine(x, y, x, y + 3);
    
    // Format the lon string nicely
    char output[128];
    if(i < 0) {
      sprintf(output, "%iW", abs(i));
    } else if(i == 0) {
      sprintf(output, "%i", i);
    } else if(i > 0) {
      sprintf(output, "%iE", i);
    }
    
    c->drawText(output, x, y + 4, Canvas::TOP, Canvas::CENTER);
  }

  c->fillRect(0 + lat_offset_x, 0 + lat_offset_y, lat_width, lat_height, 0x00FFFFFF);
    
  // Lat points
  factor = lat_height / yrange.range();;

  for(int i = (int)floorf(maxlat / lat_text_spacing) * lat_text_spacing; i > minlat; i -= lat_text_spacing) {
    int x = lat_offset_x;
    int y = lat_offset_y + (lat_height - (int)round((i - minlat) * factor));

    c->drawLine(lat_width - 4, y, lat_width - 1, y);

    // Format the lon string nicely
    char output[128];
    if(i < 0) {
      sprintf(output, "%iS", abs(i));
    } else if(i == 0) {
      sprintf(output, "%i", i);
    } else if(i > 0) {
      sprintf(output, "%iN", i);
    }

    c->drawText(output, x + lat_width - 4, y, Canvas::MIDDLE, Canvas::RIGHT);
  }
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

  // Generate the array of Y axis grid points
  for(int y = 0; y <= rows; y++) {
    ypoint[y] = (int)round(((maxlat - grid_lats[y]) / difflat) * map_height);
  }

  // Generate the array of X axis grid points
  for(int x = 0; x <= cols; x++) {
    xpoint[x] = (int)round(((grid_longs[x] - minlong) / difflong) * map_width);
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
      int* lineptr = img[y + map_offset_y] + map_offset_x;
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
      c->drawLine(map_offset_x, map_offset_y + ypoint[j], map_offset_x + map_width, map_offset_y + ypoint[j]);
    }
    for(int j = 0; j < cols; j++) {
      c->drawLine(map_offset_x + xpoint[j], map_offset_y, map_offset_x + xpoint[j], map_offset_y + map_height);
    }
  }

  delete[] data_colours;
  delete[] xpoint;
  delete[] ypoint;
}

void Displayer::drawPolygon(int numpoints, Point** points, Range xrange, Range yrange) {
  c->colour = 0x00000000;
  
  // Lines to go on graph
  // This code is acceptable now
  c->setAntiAliased();
  c->setLineThickness(LINE_WIDTH);
  int x1, x2, y1, y2;
  for(int i = 0; i < numpoints; i++) {
    // Convert to XY from latlon
    x1 = lontox(map_width, 0, xrange.min(), xrange.max(), points[i]->x);
    y1 = lattoy(map_height, 0, yrange.min(), yrange.max(), points[i]->y);
    
    x2 = lontox(map_width, 0, xrange.min(), xrange.max(), points[(i + 1) % numpoints]->x);
    y2 = lattoy(map_height, 0, yrange.min(), yrange.max(), points[(i + 1) % numpoints]->y);
    
    // Lines
    c->drawLine(x1 + map_offset_x, y1 + map_offset_y, x2 + map_offset_x, y2 + map_offset_y);
  }
  c->setLineThickness(1);
  for(int i = 0; i < numpoints; i++) {
    // Convert to XY from latlon
    x1 = lontox(map_width, 0, xrange.min(), xrange.max(), points[i]->x);
    y1 = lattoy(map_height, 0, yrange.min(), yrange.max(), points[i]->y);
    
    // Point
    if(points[i]->selected) {
      c->fillRectAbs(x1 - (POINT_SIZE - 2) + map_offset_x, y1 - (POINT_SIZE - 2) + map_offset_y, x1 + (POINT_SIZE - 2) + map_offset_x, y1 + (POINT_SIZE - 2) + map_offset_y, 0x00FF0000);
      
    } else {
      c->fillRectAbs(x1 - 3 + map_offset_x, y1 - 3 + map_offset_y, x1 + 3 + map_offset_x, y1 + 3 + map_offset_y, 0x00000000);
    }
    
    // Box around point
    c->drawRect(x1 - (POINT_SIZE - 1) + map_offset_x, y1 - (POINT_SIZE - 1) + map_offset_y, 2 * (POINT_SIZE - 1), 2 * (POINT_SIZE - 1), 0x00FFFFFF);
    c->drawRect(x1 - POINT_SIZE + map_offset_x, y1 - POINT_SIZE + map_offset_y, 2 * POINT_SIZE, 2 * POINT_SIZE);
  }      
}

void Displayer::plotDecal(string decalfile, Point dpoint, Range xrange, Range yrange) {
  // Display the data point
  // First, load the decal image (erroring if DNE)
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
	  lontox(map_width, 0, xrange.min(), xrange.max(), dpoint.x) - 
	  ((decalimg->sx - 1) >> 1) + map_offset_x, 
	  lattoy(map_height, 0, yrange.min(), yrange.max(), dpoint.y) - 
	  ((decalimg->sy - 1) >> 1) + map_offset_y, 
	  0, 0, decalimg->sx, decalimg->sy);
}

void Displayer::setRanges(const Range& datarange) {
}

Legend* Displayer::getLegend(const Range& datarange) {
  Legend* leg_colours = new Legend(datarange, colour_map, colour_map_rev);
  double leg_range_min = range_min;
  double leg_range_max = range_max;

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
  leg_dec_places = range_dynamic;
  double tmp = (leg_colours.range.max() - leg_colours.range.min()) / (num_leg_segments);
  while(tmp < 1) {
    leg_dec_places++;
    tmp *= 10;
  }

  int leg_left = LAT_WIDTH + BORDER_WIDTH + leg_offset_x;
  int leg_right = leg_width - LEG_EXTRA_RWIDTH - 2 * BORDER_WIDTH + leg_offset_x;

  int leg_top = LEG_TOP_TEXT_HEIGHT + leg_offset_y + 1;
  int leg_bottom = leg_height - LEG_BOTTOM_TEXT_HEIGHT + leg_offset_y;

  int width = leg_right - leg_left;
  int numcolours = leg_colours.numcolours();
  
  
  // White out the white areas
  c->colour = 0x00FFFFFF;
  c->fillRectAbs(0 + leg_offset_x, leg_top - 2, leg_left - 2, leg_bottom + 2);
  c->fillRectAbs(leg_right + 2, leg_top - 2, leg_width, leg_bottom + 2);
  c->fillRectAbs(0 + leg_offset_x, leg_bottom + 2, leg_width + leg_offset_x, leg_height + leg_offset_y);
  
  // Draw the title for the legend
  c->colour = 0x00000000;
  c->drawText(leg_text, leg_left + width / 2, leg_bottom + LEG_TOP_TEXT_HEIGHT, Canvas::TOP, Canvas::CENTER);
      
  // Draw the legend
  double factor = (double)width / (double)numcolours;
  for(int i = numcolours - 1; (i + 1); i--) {
    int colour = leg_colours.lookup(i);
    c->fillRectAbs((int)round(i * factor) + leg_left, leg_top, (int)round((i + 1) * factor) + leg_left, leg_bottom, colour);
  }
      
  double lbl_offset;
  if(range_dynamic) {
    lbl_offset = 0;
    factor = (double)(width) / (num_leg_segments);
  } else {
    lbl_offset = factor;
    factor = (double)(width - (2 * lbl_offset)) / num_leg_segments;
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

void Displayer::drawIdentifyText() {
  int x, y;

  // Identify text
  // Get ready to render text
  c->colour = 0x00000000;
  x = leg_offset_x;
  y = leg_offset_y + leg_height;
  c->fontsize = 8;
  c->drawText(identify_text, x, y, Canvas::BOTTOM, Canvas::LEFT);

  x = leg_offset_x + leg_width;
  y = leg_offset_y + leg_height;
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
