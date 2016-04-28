#include <assert.h>
#include <iostream>
#include <stdio.h>
#include "genimage.h"
#include "displayer.h"
#include "legend.h"
#include "range.h"
#include "datamanager.h"
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

  xaxis_width = yaxis_width + plot_width + leg_width + XAXIS_EXTRA_WIDTH;
  xaxis_height = XAXIS_HEIGHT + XAXIS_TITLE_HEIGHT + 20;

  xaxis_offset_x = 0;
  xaxis_offset_y = yaxis_height + yaxis_offset_y;

  identify_width = xaxis_width;
  identify_height = IDENTIFY_HEIGHT;

  img_width = xaxis_width;
  img_height = yaxis_height + xaxis_height + identify_height;

  leg_offset_x = plot_offset_x + plot_width + XAXIS_EXTRA_WIDTH;
  leg_offset_y = plot_offset_y;

  identify_offset_x = 0;
  identify_offset_y = xaxis_offset_y + xaxis_height;
}

// Set offsets and such
void Displayer::setBandsOffsets() {
  latlon_plot = false;

  plot_width = BAND_WIDTH + 2 * BORDER_WIDTH;
  plot_height = BAND_HEIGHT + 2 * BORDER_WIDTH;

  yaxis_width = YAXIS_WIDTH + YAXIS_TITLE_WIDTH;
  yaxis_height = plot_height + YAXIS_EXTRA_HEIGHT;

  plot_offset_x = yaxis_width;
  plot_offset_y = YAXIS_EXTRA_HEIGHT;

  yaxis_offset_x = 0;
  yaxis_offset_y = 0;

  leg_width = 0;
  leg_height = 0;

  xaxis_width = yaxis_width + plot_width + leg_width + XAXIS_EXTRA_WIDTH;
  xaxis_height = XAXIS_HEIGHT + XAXIS_TITLE_HEIGHT + 0;

  xaxis_offset_x = 0;
  xaxis_offset_y = yaxis_height + yaxis_offset_y;

  identify_width = xaxis_width;
  identify_height = IDENTIFY_HEIGHT;

  img_width = xaxis_width;
  img_height = yaxis_height + xaxis_height + identify_height;

  leg_offset_x = plot_offset_x + plot_width + XAXIS_EXTRA_WIDTH;
  leg_offset_y = plot_offset_y;

  identify_offset_x = 0;
  identify_offset_y = xaxis_offset_y + xaxis_height;
}

void Displayer::setStickplotOffsets() {
  latlon_plot = false;

  plot_width = STICK_WIDTH + 2 * BORDER_WIDTH;
  plot_height = STICK_HEIGHT + 2 * BORDER_WIDTH;

  plot_offset_x = YAXIS_WIDTH + YAXIS_TITLE_WIDTH;
  plot_offset_y = YAXIS_EXTRA_HEIGHT;

  yaxis_width = YAXIS_WIDTH + YAXIS_TITLE_WIDTH;
  yaxis_height = plot_height + YAXIS_EXTRA_HEIGHT;

  yaxis_offset_x = 0;
  yaxis_offset_y = 0;

  // FIXME:  leg_width etc probably need to be set anyways
  /*
  leg_width = SCATTER_LEG_WIDTH + 2 * BORDER_WIDTH;
  leg_height = plot_height;
  */
  xaxis_width = yaxis_width + plot_width + /* leg_width + */ XAXIS_EXTRA_WIDTH;
  xaxis_height = XAXIS_HEIGHT + XAXIS_TITLE_HEIGHT + 20;

  xaxis_offset_x = 0;
  xaxis_offset_y = yaxis_height + yaxis_offset_y;

  identify_width = xaxis_width;
  identify_height = IDENTIFY_HEIGHT;

  img_width = xaxis_width;
  img_height = yaxis_height + xaxis_height + identify_height;

  leg_offset_x = plot_offset_x + plot_width + XAXIS_EXTRA_WIDTH;
  leg_offset_y = plot_offset_y;

  identify_offset_x = 0;
  identify_offset_y = xaxis_offset_y + xaxis_height;
}

// Copy basemap to map
void Displayer::copyMap(const gdImagePtr basemap) {
  c->copy(basemap, plot_offset_x + BORDER_WIDTH, plot_offset_y + BORDER_WIDTH, 0, 0, basemap->sx, basemap->sy);
}

#define LINE_HEIGHT 11
#define BP_LINE_HEIGHT 12

// Draw a legend on the map mapping symbols to names (for plots)
void Displayer::drawBoxPlotLegend() {
  c->colour = 0x00FFFFFF;
  c->setAlpha(0);
  c->fillRect(leg_offset_x, yaxis_offset_y, leg_width, leg_offset_y - yaxis_offset_y);
  c->fillRect(leg_offset_x, leg_offset_y, leg_width, leg_height);
  c->setAlpha(1);
  c->colour = 0x00000000;

  c->fontsize = 9;
  int x = leg_offset_x;
  int y = leg_offset_y + 1;
  c->drawText("Whiskers:", x, y, Canvas::TOP, Canvas::LEFT);
  c->drawText("min/max unless + at top", x + 10, y + BP_LINE_HEIGHT * 1, Canvas::TOP, Canvas::LEFT);
  c->drawText("+/- 1.5 times IQR if + at top", x + 10, y + BP_LINE_HEIGHT * 2, Canvas::TOP, Canvas::LEFT);

  c->drawText("Box:", x, y + BP_LINE_HEIGHT * 4, Canvas::TOP, Canvas::LEFT);
  c->drawText("25th and 75th percentile", x + 10, y + BP_LINE_HEIGHT * 5, Canvas::TOP, Canvas::LEFT);
  c->drawText("line: median", x + 10, y + BP_LINE_HEIGHT * 6, Canvas::TOP, Canvas::LEFT);

  c->drawText("IQR: interquartile range", x, y + BP_LINE_HEIGHT * 9, Canvas::TOP, Canvas::LEFT);
}

// Draw a legend on the map mapping symbols to names (for plots)
void Displayer::drawLegend(list<LegendToken* >& vars) {
  const int true_leg_height = vars.size() * LINE_HEIGHT;
  const int true_leg_width = leg_width - 2 * BORDER_WIDTH;

  const int true_leg_offset_x = leg_offset_x + BORDER_WIDTH;
  const int true_leg_offset_y = leg_offset_y + BORDER_WIDTH;

  c->colour = 0x00FFFFFF;
  c->setAlpha(0);
  c->fillRect(leg_offset_x, yaxis_offset_y, leg_width, leg_offset_y - yaxis_offset_y);
  c->fillRect(true_leg_offset_x, true_leg_offset_y, true_leg_width, true_leg_height);
  c->fillRect(leg_offset_x, true_leg_offset_y + true_leg_height + BORDER_WIDTH,
              leg_width, leg_height - true_leg_height - BORDER_WIDTH * 2);

  c->colour = 0x00000000;
  c->drawRect(true_leg_offset_x - BORDER_WIDTH, true_leg_offset_y - BORDER_WIDTH,
              true_leg_width + 2 * BORDER_WIDTH, true_leg_height + 2 * BORDER_WIDTH);
  c->setAlpha(1);


  list<LegendToken* >::iterator liter = vars.begin();

  c->fontsize = 8;
  int x = true_leg_offset_x + 20;
  int y = true_leg_offset_y + 1;
  for (; liter != vars.end(); liter++) {
    c->colour = (*liter)->colour;
    if ((*liter)->filled) {
      c->fillSymbol((*liter)->sym, x - 15, y);
    } else {
      c->drawSymbol((*liter)->sym, x - 15, y);
    }
    c->colour = 0x00000000;
    if ((*liter)->sym != NONE) {
      c->drawText((*liter)->name, x, y, Canvas::TOP, Canvas::LEFT);
      y += LINE_HEIGHT;
    }
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


// JMS: FIXME:  Fixed temporarily, needs a looking over!
int dec_places_needed(double value, int dec_places) {
  value = fabs(value);
  if (value == 0.0)
    return 1;

  while (value < 1) {
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
  for (double i = floor(max_x / x_text_spacing) * x_text_spacing; i >= min_x; i -= x_text_spacing) {
    const int x = plot_offset_x + (int)round((i - min_x) * factor);
    const int y = plot_offset_y;

    // Draw vertical lines
    if (i > -1E-12 && i < 1E-12) {
      c->drawLine(x, y, x, y + plot_height - 1, 0x00808080);
    } else {
      c->drawLine(x, y, x, y + plot_height - 1);
    }
  }

  // Lat points
  factor = (plot_height - 1) / yrange.range();

  for (double i = floor(max_y / y_text_spacing) * y_text_spacing; i >= min_y; i -= y_text_spacing) {
    const int x = plot_offset_x;
    const int y = plot_offset_y + (int)round((max_y - i) * factor);

    // Draw horizontal lines
    if (i > -1E-12 && i < 1E-12) {
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
  c->drawText(yaxis_text, yaxis_offset_x + 5, plot_offset_y + ((plot_height - 1) / 2), Canvas::MIDDLE, Canvas::LEFT, 0.5 * M_PI);
}

void Displayer::fillTickAreas() {
  c->setAlpha(0);
  c->fillRect(xaxis_offset_x, xaxis_offset_y, xaxis_width, xaxis_height, 0x00FFFFFF);
  c->fillRect(yaxis_offset_x, yaxis_offset_y, yaxis_width, yaxis_height, 0x00FFFFFF);
  c->setAlpha(1);
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

  fillTickAreas();
  c->colour = 0x00000000;
  c->fontsize = 9;

  // Set the legend's decimal places
  sprintf(format, "%%0.%if%%s", dec_places_needed(x_text_spacing, 0));

  for (double i = floor(max_x / x_text_spacing) * x_text_spacing; i >= min_x; i -= x_text_spacing) {
    const int x = plot_offset_x + (int)round((i - min_x) * factor);
    const int y = xaxis_offset_y;

    // Put bottom tick marks on
    c->drawLine(x, y, x, y + 3);

    // Format the lon string nicely
    char output[32];

    if (!latlon_plot || i == 0) {
      sprintf(output, format, i, "");
    } else if (i < 0) {
      sprintf(output, format, -i, "W");
    } else if (i > 0) {
      sprintf(output, format, i, "E");
    }

    c->drawText(output, x, y + 4, Canvas::TOP, Canvas::CENTER);
  }

  // Lat points
  // Set the legend's decimal places
  sprintf(format, "%%0.%if%%s", dec_places_needed(y_text_spacing, 0));

  factor = (plot_height - 1) / yrange.range();

  for (double i = floor(max_y / y_text_spacing) * y_text_spacing; i >= min_y; i -= y_text_spacing) {
    const int x = yaxis_offset_x;
    const int y = plot_offset_y + (int)round((max_y - i) * factor);

    c->drawLine(yaxis_width - 4, y, yaxis_width - 1, y);

    // Format the lon string nicely
    char output[128];
    if (!latlon_plot || i == 0) {
      sprintf(output, format, i, "");
    } else if (i < 0) {
      sprintf(output, format, -i, "S");
    } else if (i > 0) {
      sprintf(output, format, i, "N");
    }

    c->drawText(output, x + yaxis_width - 4, y, Canvas::MIDDLE, Canvas::RIGHT);
  }
}

// Draw (and possibly label) tickmarks on the map/plot -- Y AXIS ONLY
void Displayer::drawYTicks(const Range& xrange, const Range& yrange) {
  // Lon grid
  double factor = (plot_width - 1) / xrange.range();
  //  double min_x = xrange.min();
  //  double max_x = xrange.max();
  double min_y = yrange.min();
  double max_y = yrange.max();
  char format[32];

  fillTickAreas();
  c->colour = 0x00000000;
  c->fontsize = 10;

  // Set the legend's decimal places
  sprintf(format, "%%0.%if%%s", dec_places_needed(x_text_spacing, 0));
  /*
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
        sprintf(output, format, -i, "W");
      } else if(i > 0) {
        sprintf(output, format, i, "E");
      }

      c->drawText(output, x, y + 4, Canvas::TOP, Canvas::CENTER);
    }
  */
  // Lat points
  // Set the legend's decimal places
  sprintf(format, "%%0.%if%%s", dec_places_needed(y_text_spacing, 0));

  factor = (plot_height - 1) / yrange.range();

  for (double i = floor(max_y / y_text_spacing) * y_text_spacing; i >= min_y; i -= y_text_spacing) {
    const int x = yaxis_offset_x;
    const int y = plot_offset_y + (int)round((max_y - i) * factor);

    c->drawLine(yaxis_width - 4, y, yaxis_width - 1, y);

    // Format the lon string nicely
    char output[128];
    if (!latlon_plot || i == 0) {
      sprintf(output, format, i, "");
    } else if (i < 0) {
      sprintf(output, format, -i, "S");
    } else if (i > 0) {
      sprintf(output, format, i, "N");
    }

    c->drawText(output, x + yaxis_width - 4, y, Canvas::MIDDLE, Canvas::RIGHT);
  }
}

// Clear entire canvas.  Only use this if necessary, as it's inefficient, but it's the elegant thing to do for some plot types.
void Displayer::clearCanvas() {
  c->colour = 0x00FFFFFF;
  c->setAlpha(0);
  c->fillRect(0, 0, img_width, img_height);
  c->setAlpha(1);
}

// Clear plot
void Displayer::clearPlot() {
  c->colour = 0x00FFFFFF;
  c->setAlpha(0);
  c->fillRect(plot_offset_x, plot_offset_y, plot_width, plot_height);
  c->setAlpha(1);
}

// Draw scatter plot
void Displayer::drawScatter(std::list<ScatterVars*>& vars, const Range& xrange, const Range& yrange, bool stick) {
  list<ScatterVars*>::iterator i = vars.begin();
  double xfactor = (double)(plot_width) / xrange.range();
  double yfactor = (double)(plot_height) / yrange.range();

  if (!stick) { // normal
    for (; i != vars.end(); i++) {
      int x = (int)round(((*i)->datx - xrange.min()) * xfactor) + plot_offset_x;
      int y = (int)round((yrange.max() - (*i)->daty) * yfactor) + plot_offset_y;

      c->colour = (*i)->symbol->colour;
      if ((*i)->symbol->filled) {
        c->fillSymbol((*i)->symbol->sym, x - 5, y - 5);
      } else {
        c->drawSymbol((*i)->symbol->sym, x - 5, y - 5);
      }
    }
  } else { // draw things with colours last! -- FIXME this is a terrible kludge but sorting a set also makes no sense -- should store the one point separately.
    for (; i != vars.end(); i++) {
      if ((*i)->symbol->colour == TGRAY /* ew, ew, ew, and more ew.  I should be checking this properly but there's no other flag. */) {
        int x = (int)round(((*i)->datx - xrange.min()) * xfactor) + plot_offset_x;
        int y = (int)round((yrange.max() - (*i)->daty) * yfactor) + plot_offset_y;

        c->colour = (*i)->symbol->colour; // TODO still need to do alpha blending here
        if ((*i)->symbol->filled) {
          c->fillSymbol((*i)->symbol->sym, x - 5, y - 5);
        } else {
          c->drawSymbol((*i)->symbol->sym, x - 5, y - 5);
        }
      }
    }
    i = vars.begin();
    for (; i != vars.end(); i++) { // now draw things with colours
      if ((*i)->symbol->colour != TGRAY) {
        int x = (int)round(((*i)->datx - xrange.min()) * xfactor) + plot_offset_x;
        int y = (int)round((yrange.max() - (*i)->daty) * yfactor) + plot_offset_y;

        c->colour = (*i)->symbol->colour;
        if ((*i)->symbol->filled) {
          c->fillSymbol((*i)->symbol->sym, x - 5, y - 5);
        } else {
          c->drawSymbol((*i)->symbol->sym, x - 5, y - 5);
        }
      }
    }
  }

  // Clean up the rect around the map
  if (!stick)
    c->drawRect(plot_offset_x, plot_offset_y, plot_width, plot_height, 0x00000000);
}

static bool compareByTimesliceAndY(const ScatterVars* a, const ScatterVars* b) {
  if (a->datx == b->datx)
    return a->daty < b->daty;
  return a->datx < b->datx;
}

// Draw bands plot
void Displayer::drawBands(std::list<ScatterVars*>& vars, const Range& xrange, const Range& yrange, const DataSpec& s) {
  list<ScatterVars*> median_line;

  // Screen coordinate max and min
  enum { MAX = 0, TOP, MID, BOT, MIN };
  // Real coordinate max and min
  enum { CMIN = 0, CBOT, CMID, CTOP, CMAX };

  std::list< ScatterVars*> pruned;
  std::list<ScatterVars*> chosen_model;
  for (list<ScatterVars*>::iterator it = vars.begin(); it != vars.end(); ++it) {
    ScatterVars* sv = *it;
    if ( sv->symbol->sym != NONE )
      pruned.push_back(sv);
    if ( sv->model == s.model && sv->expt == s.expt )
      chosen_model.push_back(sv);

  }
  fprintf(stderr, "chosen model list size: %zi\n", chosen_model.size());
  assert(chosen_model.size() == 4);
  chosen_model.sort(compareByTimesliceAndY);
  assert(pruned.size() % 5 == 0);
  pruned.sort(compareByTimesliceAndY);

  double xfactor = (double)(plot_width) / xrange.range();
  double yfactor = (double)(plot_height) / yrange.range();

  list<ScatterVars*>::iterator it = pruned.begin();
  ScatterVars* prev_coords[5];
  ScatterVars* coords[5];
  int prev_y[5];
  int y[5];
  for (int j = MAX; j <= MIN; ++j) {
    coords[j] = *it++;
    y[j] = (int)round((yrange.max() - coords[j]->daty) * yfactor) + plot_offset_y;
  }
  median_line.push_back(coords[CMID]);
  int x = (int)round((coords[0]->datx - xrange.min()) * xfactor) + plot_offset_x;
  int prev_x;
  for (int i = 1; i < floor(pruned.size() / 5); ++i) {
    for (int j = MAX; j <= MIN; ++j) {
      prev_coords[j] = coords[j];
      coords[j] = *it++;
      prev_y[j] = y[j];
      y[j] = (int)round((yrange.max() - coords[j]->daty) * yfactor) + plot_offset_y;
    }
    assert(y[MAX] >= y[TOP] && y[TOP] >= y[MID] && y[MID] >= y[BOT] && y[BOT] >= y[MIN]);

    prev_x = x;
    x = (int)round((coords[0]->datx - xrange.min()) * xfactor) + plot_offset_x;

    median_line.push_back(coords[CMID]);
    gdPoint quad1[] = { {prev_x, prev_y[CMIN]}, {x, y[CMIN]}, {x, y[CBOT]}, {prev_x, prev_y[CBOT]} };
    gdPoint quad2[] = { {prev_x, prev_y[CBOT]}, {x, y[CBOT]}, {x, y[CTOP]}, {prev_x, prev_y[CTOP]} };
    gdPoint quad3[] = { {prev_x, prev_y[CTOP]}, {x, y[CTOP]}, {x, y[CMAX]}, {prev_x, prev_y[CMAX]} };
    c->fillPoly(quad1, 4, LTGRAY);
    c->fillPoly(quad2, 4, DKGRAY);
    c->fillPoly(quad3, 4, LTGRAY);
  }

  drawScatterGrid(xrange, yrange);

  // Draw median line
  drawLines(median_line, xrange, yrange);

  // Draw chosen model line
  drawLines(chosen_model, xrange, yrange);

  // Clean up the rect around the map
  c->drawRect(plot_offset_x, plot_offset_y, plot_width, plot_height, 0x00000000);
}

// Draw box plot
#define WIDTH 10
#define BOXWIDTH 80
void Displayer::drawBoxPlot(std::list<ScatterVars*>& vars, const Range& xrange, const Range& yrange) {
  double xfactor = (double)(plot_width) / xrange.range();
  double yfactor = (double)(plot_height) / yrange.range();

  std::list<ScatterVars*> pruned;
  for (list<ScatterVars*>::iterator it = vars.begin(); it != vars.end(); ++it) {
    if ( (*it)->symbol->sym != NONE )
      pruned.push_back(*it);
  }
  fprintf(stderr, "pruned size: %zi\n", pruned.size());
  assert(pruned.size() % 5 == 0);
  pruned.sort(compareByTimesliceAndY);

  list<ScatterVars*>::iterator it = pruned.begin();
  for (int i = 0; i < floor(pruned.size() / 5); ++i) {
    // Screen coordinate max and min
    enum { MAX = 0, TOP, MID, BOT, MIN };

    // Real coordinate max and min
    enum { CMIN = 0, CBOT, CMID, CTOP, CMAX };
    ScatterVars* coords[5];
    int y[5];

    for (int j = MAX; j <= MIN; ++j) {
      coords[j] = *it++;
      y[j] = (int)round((yrange.max() - coords[j]->daty) * yfactor) + plot_offset_y;
    }
    assert(y[MAX] >= y[TOP] && y[TOP] >= y[MID] && y[MID] >= y[BOT] && y[BOT] >= y[MIN]);

    int x = (int)round((coords[0]->datx - xrange.min()) * xfactor) + plot_offset_x;

    double iqr = coords[CTOP]->daty - coords[CBOT]->daty;

    // If the max is more than 1.5IQR from the top of the box...
    if (coords[CMAX]->daty - coords[CTOP]->daty > 1.5 * iqr) {
      // Set the max to top + 1.5 IQR and plot a plus
      y[MIN] = (int)round((yrange.max() - (coords[CTOP]->daty + 1.5 * iqr)) * yfactor) + plot_offset_y;
      c->drawLine(x - 5, y[MIN] - 10, x + 5, y[MIN] - 10, 0x00000000);
      c->drawLine(x, y[MIN] - 15, x, y[MIN] - 5, 0x00000000);
    }

    // If the max is more than 1.5IQR from the top of the box...
    if (coords[CBOT]->daty - coords[CMIN]->daty > 1.5 * iqr) {
      // Set the max to top + 1.5 IQR and plot a plus
      y[MAX] = (int)round((yrange.max() - (coords[CBOT]->daty - 1.5 * iqr)) * yfactor) + plot_offset_y;
      c->drawLine(x - 5, y[MAX] + 10, x + 5, y[MAX] + 10, 0x00000000);
      c->drawLine(x, y[MAX] + 15, x, y[MAX] + 5, 0x00000000);
    }

    // Plot box, lines, etc
    c->drawLine(x - WIDTH, y[MAX], x + WIDTH, y[MAX], 0x00000000);
    c->drawLine(x, y[MAX], x, y[TOP], 0x00000000);
    c->drawRectAbs(x - BOXWIDTH, y[TOP], x + BOXWIDTH, y[MID], 0x00000000);
    c->drawRectAbs(x - BOXWIDTH, y[MID], x + BOXWIDTH, y[BOT], 0x00000000);
    c->drawLine(x, y[BOT], x, y[MIN], 0x00000000);
    c->drawLine(x - 10, y[MIN], x + 10, y[MIN], 0x00000000);
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
  for (; i != vars.end(); i++) {
    int x = (int)round(((*i)->datx - xrange.min()) * xfactor) + plot_offset_x;
    int y = (int)round((yrange.max() - (*i)->daty) * yfactor) + plot_offset_y;

    if (old && (old->model == (*i)->model && old->expt == (*i)->expt)) {
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
void Displayer::drawMap(const gdImagePtr basemap, const Legend& leg_colours, const DataGrid<int>& draw_mask, const DataGrid<int>& data_mask, const DataGrid<double>& data, const Window& w) {
  const int x_size = data.x_size();
  const int y_size = data.y_size();
  const double* x_grid_orig = data.x_grid().get();
  const double* y_grid = data.y_grid().get();
  const double missing = data.missing();
  const int xpsize = x_size * 2;
  const int ypsize = y_size * 2;
  int* xpoint = new int[xpsize + 2];
  int* ypoint = new int[ypsize];
  const int datasize = data.grid_size();

  const double minlong = w.left;
  const double difflong = w.right - w.left;
  const double maxlat = w.top;
  const double difflat = w.top - w.bottom;

  fprintf(stderr, "left: %f, right: %f, top: %f, bottom: %f\n", w.left, w.right, w.top, w.bottom);
  fprintf(stderr, "diff_y: %f, diff_x: %f\n", difflat, difflong);

  unsigned int* data_colours = new unsigned int[datasize];
  unsigned char** map = basemap->pixels;
  int** img = c->img->tpixels;

  const int true_plot_height = plot_height - 2 * BORDER_WIDTH;
  const int true_plot_width = plot_width - 2 * BORDER_WIDTH;

  const int true_plot_offset_x = plot_offset_x + BORDER_WIDTH;
  const int true_plot_offset_y = plot_offset_y + BORDER_WIDTH;

  // IDEA: If rows > plot width, invert xpoint stuff to have each pixel refer to a data offset instead in the X dimension

  // Generate the array of Y axis grid points
  for (int y = 0; y < ypsize; y++) {
    ypoint[y] = (int)round((clip_to_range(0.0, difflat, (maxlat - y_grid[y])) / difflat) * true_plot_height);
  }

  // Generate the array of X axis grid points
  const double* x_grid;
  if (data.projection() == "equidistant_cylindrical") {
    double* x_grid_new = new double[data.xgrid_size()];
    std::copy(x_grid_orig, &x_grid_orig[data.xgrid_size()], x_grid_new);
    shift_longs(x_grid_new, data.xgrid_size(), (w.left + w.right) / 2);
    x_grid = x_grid_new;
  } else {
    x_grid = x_grid_orig;
  }
  for (int x = 0; x < xpsize; x++) {
    xpoint[x] = (int)round((clip_to_range(0.0, difflong, (x_grid[x] - minlong)) / difflong) * true_plot_width);
  }

  // Wrap around if start and end gridlongs match
  xpoint[xpsize] = xpoint[xpsize - 1];
  if (x_grid[0] + 360 == x_grid[xpsize - 1]) {
    xpoint[xpsize + 1] = true_plot_width;
  } else {
    xpoint[xpsize + 1] = xpoint[xpsize];
  }

  if (data.projection() == "equidistant_cylindrical") {
    delete[] x_grid;
  }

  assert(ypoint[0] < ypoint[ypsize - 1]); // ypoints should increase

  //DEBUG

  fprintf(stderr, "xpoint values:");  for (int i = 0; i < xpsize + 2; i++) fprintf(stderr, " %d", xpoint[i]);  fprintf(stderr, "\n");

  //DEBUG
  if (!(xpoint[0] < xpoint[xpsize - 1]))
    fprintf(stderr, "Error: xpoint[0] (%d) >= xpoint[xpsize - 1 (%d)] (%d)\n", xpoint[0], xpsize - 1, xpoint[xpsize - 1]);


  assert(xpoint[0] < xpoint[xpsize - 1]); // ditto

  // Get the data colours
  for (int j = 0; j < y_size; j++) {
    if (ypoint[j * 2] != ypoint[j * 2 + 1]) {
      unsigned int* data_colours_line = &data_colours[j * x_size];
      const double* data_line = &data.values().get()[j * x_size];
      const int* datamask_line = &data_mask.values().get()[j * x_size];
      const int* drawmask_line = &draw_mask.values().get()[j * x_size];
      for (int i = 0; i < x_size; i++) {
        if (xpoint[i * 2] != xpoint[i * 2 + 1]) {
          if (drawmask_line[i] && data_line[i] != missing) {
            data_colours_line[i] = leg_colours.lookup(data_line[i]);
          } else {
            data_colours_line[i] = 0x00FFFFFF;
          }
          if (!datamask_line[i]) {
            // Half-brightness colour
            data_colours_line[i] = (data_colours_line[i] >> 1) & 0x7F7F7F7F;
          }
        }
      }
    }
  }

  // Plot the data on the map
  for (int j = 0; j < y_size; j++) {
    int y = ypoint[j * 2];
    const int yend = ypoint[j * 2 + 1];
    const unsigned int* data_colours_ptr = &data_colours[j * x_size];
    for (; y < yend; y++) {
      const unsigned char* mlineptr = map[y];
      int* lineptr = img[y + true_plot_offset_y] + true_plot_offset_x;
      for (int i = 0; i <= x_size; i++) {
        // Look up colour
        int x = xpoint[i * 2];
        const int xend = xpoint[i * 2 + 1];
        const unsigned int colour = data_colours_ptr[i % x_size];
        // Map is either 1 (white) or 0 (black)
        // Negation of 1 is 0xffffffff; negation of 0 is 0x00000000
        // AND'ing that with the colour gives either black or the colour
        for (; x < xend; x++) {
          // Apply colour to map
          lineptr[x] = colour & -((int)mlineptr[x] ^ 1);
        }
      }
    }
  }

  // Copy any remaining basemap to the map
  if (ypoint[0] > 0) {
    c->copy(basemap, true_plot_offset_x, true_plot_offset_y, 0, 0, true_plot_width, ypoint[0]);
  }
  if (xpoint[0] > 0) {
    c->copy(basemap, true_plot_offset_x, true_plot_offset_y, 0, 0, xpoint[0], true_plot_height);
  }
  if (ypoint[ypsize - 1] < true_plot_height) {
    c->copy(basemap, true_plot_offset_x, true_plot_offset_y + ypoint[ypsize - 1], 0, ypoint[ypsize - 1], true_plot_width, true_plot_height - ypoint[ypsize - 1]);
  }
  if (xpoint[xpsize + 1] < true_plot_width) {
    c->copy(basemap, true_plot_offset_x + xpoint[xpsize + 1], true_plot_offset_y, xpoint[xpsize + 1], 0, true_plot_width - xpoint[xpsize + 1], true_plot_height);
  }

  // Fill in the center block if there's a gap
  //for(int i = 1; i <= y_size; i++) {
  // if(xpoint[i * 2 - 1] < xpoint[i * 2]) {
  //   c->copy(basemap, true_plot_offset_x + xpoint[i * 2 - 1], true_plot_offset_y + ypoint[0], xpoint[i * 2 - 1], ypoint[0], xpoint[i * 2] - xpoint[i * 2 - 1], ypoint[ypsize - 1] - ypoint[0]);
  //  }
  //}

  // After we're done with the map, we put the grid on it
  c->colour = 0x00A0A0A0;
  if (grid) {
    for (int j = 1; j < y_size; j++) {
      c->drawLine(true_plot_offset_x, true_plot_offset_y + ypoint[j], true_plot_offset_x + plot_width, true_plot_offset_y + ypoint[j]);
    }
    for (int j = 0; j < x_size; j++) {
      c->drawLine(true_plot_offset_x + xpoint[j], true_plot_offset_y, true_plot_offset_x + xpoint[j], true_plot_offset_y + plot_height);
    }
  }

  delete[] data_colours;
  delete[] xpoint;
  delete[] ypoint;
}

void Displayer::drawPolygon(const vector<Point>& points, const Range xrange, const Range yrange, bool draw_vertices) {
  const int numpoints = points.size();
  c->colour = 0x00000000;

  const int true_plot_height = plot_height - 2 * BORDER_WIDTH;
  const int true_plot_width = plot_width - 2 * BORDER_WIDTH;

  const int true_plot_offset_x = plot_offset_x + BORDER_WIDTH;
  const int true_plot_offset_y = plot_offset_y + BORDER_WIDTH;

  // Lines to go on graph
  c->setClip(true_plot_offset_x, true_plot_offset_y, true_plot_offset_x + true_plot_width, true_plot_offset_y + true_plot_height);
  c->setAntiAliased();
  c->setLineThickness(LINE_WIDTH);
  int x1, x2, y1, y2;
  for (int i = 0; i < numpoints; i++) {
    // Convert to XY from latlon
    x1 = lontox(true_plot_width, 0, xrange.min(), xrange.max(), points[i].x);
    y1 = lattoy(true_plot_height, 0, yrange.min(), yrange.max(), points[i].y);

    x2 = lontox(true_plot_width, 0, xrange.min(), xrange.max(), points[(i + 1) % numpoints].x);
    y2 = lattoy(true_plot_height, 0, yrange.min(), yrange.max(), points[(i + 1) % numpoints].y);

    // Lines
    c->drawLine(x1 + true_plot_offset_x, y1 + true_plot_offset_y, x2 + true_plot_offset_x, y2 + true_plot_offset_y);
  }
  if (draw_vertices) { // May need to draw circles with thicker polygon edges; right now LINE_WIDTH is 3.
    c->setLineThickness(1);
    for (int i = 0; i < numpoints; i++) {
      // Convert to XY from latlon
      x1 = lontox(plot_width, 0, xrange.min(), xrange.max(), points[i].x);
      y1 = lattoy(plot_height, 0, yrange.min(), yrange.max(), points[i].y);

      // Point
      if (points[i].selected) {
        c->fillRectAbs(x1 - (POINT_SIZE - 2) + true_plot_offset_x, y1 - (POINT_SIZE - 2) + true_plot_offset_y, x1 + (POINT_SIZE - 2) + true_plot_offset_x, y1 + (POINT_SIZE - 2) + true_plot_offset_y, 0x00FF0000);

      } else {
        c->fillRectAbs(x1 - 3 + true_plot_offset_x, y1 - 3 + true_plot_offset_y, x1 + 3 + true_plot_offset_x, y1 + 3 + true_plot_offset_y, 0x00000000);
      }

      // Box around point
      c->drawRect(x1 - (POINT_SIZE - 1) + true_plot_offset_x, y1 - (POINT_SIZE - 1) + true_plot_offset_y, 2 * (POINT_SIZE - 1), 2 * (POINT_SIZE - 1), 0x00FFFFFF);
      c->drawRect(x1 - POINT_SIZE + true_plot_offset_x, y1 - POINT_SIZE + true_plot_offset_y, 2 * POINT_SIZE, 2 * POINT_SIZE);
    }
  }
  c->setClip(0, 0, img_width, img_height);
}


Legend* Displayer::getLegend(const Range& datarange) {
  Range r(datarange);
  double leg_range_min = xrange_min;
  double leg_range_max = xrange_max;

  if (range_dynamic) {
    leg_range_min = datarange.min();
    leg_range_max = datarange.max();
  }
  r.setmin(leg_range_min);
  r.setmax(leg_range_max);

  return new Legend(r, colour_map, colour_map_rev);
}

void Displayer::drawScale(Legend& leg_colours) {
  // Generate the colour map
  int num_leg_segments = 10;
  int leg_dec_places;
  int x, y;

  if (range_dynamic && colour_map == STEPWISE) {
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
  for (int i = numcolours - 1; (i + 1); i--) {
    int colour = leg_colours.lookup(i);
    c->fillRectAbs((int)round(i * factor) + leg_left, leg_top, (int)round((i + 1) * factor) + leg_left, leg_top + (int)ceil(height / 2), colour);
    c->fillRectAbs((int)round(i * factor) + leg_left, leg_top + (int)ceil(height / 2), (int)round((i + 1) * factor) + leg_left, leg_bottom, (colour >> 1) & 0x7F7F7F7F);
  }

  double lbl_offset;
  if (range_dynamic) {
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
  for (int i = 0; i <= num_leg_segments; i++) {
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
  y = identify_offset_y + 5;
  c->fontsize = 10;
  c->drawText(credit_text, x, y, Canvas::BOTTOM, Canvas::RIGHT);

  c->colour = 0x00FF0000;
  x = plot_offset_x + 2;
  y = plot_offset_y + 2;
  c->fontsize = 14;
  c->drawText(warning, x, y, Canvas::TOP, Canvas::LEFT);
}

bool Displayer::writePng(string filename) {
  FILE* f;
  f = fopen(filename.c_str(), "wb");
  if (!f) {
    fprintf(stderr, "Could not open output file!\n");
    return false;
  }
  gdImagePng(c->img, f);
  fclose(f);
  return true;
}
