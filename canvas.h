#ifndef __GENIMAGE_CANVAS_H
#define __GENIMAGE_CANVAS_H
#include <gd.h>
#include <string>

#include "legendtoken.h"

#define POLY_SIZE 11

enum STYLES{DASHED};

// BIG WORK SAVER: gdImageSetClip
// Look into it... could clean up code a -lot-

class Canvas {
public:
  enum VAlignment { TOP, MIDDLE, BOTTOM };
  enum HAlignment { LEFT, CENTER, RIGHT };
  Canvas();
  Canvas(Canvas& c);
  Canvas(std::string fontfile, int width, int height, gdImagePtr img = 0);
  ~Canvas();
  bool writeImage(char* filename);
  bool writeImage(FILE* f);
  void setLineThickness(int width);
  void setAntiAliased();
  void setAntiAliased(int colour);
  void drawLine(int x1, int y1, int x2, int y2);
  void drawLine(int x1, int y1, int x2, int y2, int colour);
  void copy(gdImagePtr src, int dx, int dy, int sx, int sy, int width, int height);
  void copy(Canvas& src, int dx, int dy, int sx, int sy, int width, int height);
  void drawRect(int x, int y, int width, int height);
  void drawRect(int x, int y, int width, int height, int colour);
  void drawRectAbs(int x1, int y1, int x2, int y2);
  void drawRectAbs(int x1, int y1, int x2, int y2, int colour);
  void fillRect(int x, int y, int width, int height);
  void fillRect(int x, int y, int width, int height, int colour);
  void fillRectAbs(int x1, int y1, int x2, int y2);
  void fillRectAbs(int x1, int y1, int x2, int y2, int colour);

  void drawText(char* s, int x, int y, VAlignment v, HAlignment h, double angle = 0);
  void drawText(char* s, int x, int y, VAlignment v, HAlignment h, int colour, int size, double angle = 0);

  void drawText(std::string s, int x, int y, VAlignment v, HAlignment h, double angle = 0);
  void drawText(std::string s, int x, int y, VAlignment v, HAlignment h, int colour, int size, double angle = 0);

  void setOffsets(int offset_x, int offset_y);

  void setClip(int x1, int y1, int x2, int y2);

  void setAlpha(int alpha) {
    gdImageAlphaBlending(img, alpha);
  }

  void setStyle(enum STYLES s) {
    int grey = 0x00808080;
    const int dashed[] = { gdTransparent, gdTransparent, grey, grey, gdTransparent, gdTransparent };
    switch(s) {
    case DASHED:
      gdImageSetStyle(img, (int*)dashed, 6);
      break;
    }
  }

  void drawSymbol(enum SYMBOL s, int x, int y);
  void fillSymbol(enum SYMBOL s, int x, int y);
  
  /* Image to draw on */
  gdImagePtr img;
  bool hasParent;

  /* Offset to draw at */
  int offset_x;
  int offset_y;

  /* Duplication of canvas data */
  int width;
  int height;

  /* Colour */
  int colour;

  /* Font options */
  std::string fontfile;
  int fontsize;
};

#endif
