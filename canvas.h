#ifndef __GENIMAGE_CANVAS_H
#define __GENIMAGE_CANVAS_H
#include <gd.h>
#include <string>

using namespace std;

class Canvas {
public:
  enum VAlignment { TOP, MIDDLE, BOTTOM };
  enum HAlignment { LEFT, CENTER, RIGHT };
  Canvas();
  Canvas(Canvas& c);
  Canvas(string fontfile, int width, int height, gdImagePtr img = 0);
  ~Canvas();
  bool writeImage(char* filename);
  bool writeImage(FILE* f);
  void setLineThickness(int width);
  void setAntiAliased();
  void setAntiAliased(int colour);
  void drawLine(int x1, int y1, int x2, int y2);
  void drawLine(int x1, int y1, int x2, int y2, int colour);
  void copy(gdImagePtr src, int dx, int dy, int sx, int sy, int width, int height);
  void drawRect(int x, int y, int width, int height);
  void drawRect(int x, int y, int width, int height, int colour);
  void drawRectAbs(int x1, int y1, int x2, int y2);
  void drawRectAbs(int x1, int y1, int x2, int y2, int colour);
  void fillRect(int x, int y, int width, int height);
  void fillRect(int x, int y, int width, int height, int colour);
  void fillRectAbs(int x1, int y1, int x2, int y2);
  void fillRectAbs(int x1, int y1, int x2, int y2, int colour);
  void drawText(char* s, int x, int y, VAlignment v, HAlignment h);
  void drawText(char* s, int x, int y, VAlignment v, HAlignment h, int colour, int size);
  void drawText(string s, int x, int y, VAlignment v, HAlignment h);
  void drawText(string s, int x, int y, VAlignment v, HAlignment h, int colour, int size);
  void setOffsets(int offset_x, int offset_y);
  
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
  string fontfile;
  int fontsize;
};

#endif
