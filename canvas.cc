#include <iostream>
#include <string>
#include <math.h>
#include "canvas.h"

Canvas::Canvas() {
  width = 0;
  height = 0;
  offset_x = 0;
  offset_y = 0;
  hasParent = false;
}

Canvas::Canvas(Canvas& c) {
  this->img = c.img;
  this->hasParent = c.hasParent;
  this->offset_x = c.offset_x;
  this->offset_y = c.offset_y;
  this->width = c.width;
  this->height = c.height;
  this->colour = c.colour;
  this->fontfile = c.fontfile;
  this->fontsize = c.fontsize;
}

Canvas::Canvas(std::string fontfile, int width, int height, gdImagePtr img) {
  this->width = width;
  this->height = height;
  this->fontfile = fontfile;
  offset_x = 0;
  offset_y = 0;
  hasParent = true;
  if(!img) {
    this->img = gdImageCreateTrueColor(width, height);
    hasParent = false;
  }
}

Canvas::~Canvas() {
  if(!hasParent) {
    gdImageDestroy(img);
  }
}

bool Canvas::writeImage(char* filename) {
  FILE* f = fopen(filename, "wb");
  if(!f) {
    return false;
  }
  gdImagePng(img, f);
  fclose(f);
  return true;
}

bool Canvas::writeImage(FILE* f) {
  if(!f) {
    return false;
  }
  gdImagePng(img, f);
  return true;
}

void Canvas::setLineThickness(int width) {
  gdImageSetThickness(img, width);
}

void Canvas::setAntiAliased() {
  gdImageSetAntiAliased(img, colour);
}

void Canvas::setAntiAliased(int colour) {
  gdImageSetAntiAliased(img, colour);
}

void Canvas::drawLine(int x1, int y1, int x2, int y2) {
  gdImageLine(img, x1 + offset_x, y1 + offset_y, x2 + offset_x, y2 + offset_y, colour);
}

void Canvas::drawLine(int x1, int y1, int x2, int y2, int colour) {
  gdImageLine(img, x1 + offset_x, y1 + offset_y, x2 + offset_x, y2 + offset_y, colour);
}

void Canvas::copy(gdImagePtr src, int dx, int dy, int sx, int sy, int width, int height) {
  gdImageCopy(img, src, dx + offset_x, dy + offset_y, sx, sy, width, height);
}

void Canvas::copy(Canvas& src, int dx, int dy, int sx, int sy, int width, int height) {
  gdImageCopy(img, src.img, dx + offset_x, dy + offset_y, sx, sy, width, height);
}

void Canvas::drawRect(int x, int y, int width, int height) {
  gdImageRectangle(img, x + offset_x, y + offset_y, x + offset_x + width - 1, y + offset_y + height - 1, colour);
}

void Canvas::drawRect(int x, int y, int width, int height, int colour) {
  gdImageRectangle(img, x + offset_x, y + offset_y, x + offset_x + width - 1, y + offset_y + height - 1, colour);
}

void Canvas::drawRectAbs(int x1, int y1, int x2, int y2) {
  gdImageRectangle(img, x1 + offset_x, y1 + offset_y, x2 + offset_x, y2 + offset_y, colour);
}

void Canvas::drawRectAbs(int x1, int y1, int x2, int y2, int colour) {
  gdImageRectangle(img, x1 + offset_x, y1 + offset_y, x2 + offset_x, y2 + offset_y, colour);
}

void Canvas::fillRect(int x, int y, int width, int height) {
  gdImageFilledRectangle(img, x + offset_x, y + offset_y, x + offset_x + width - 1, y + offset_y + height - 1, colour);
}

void Canvas::fillRect(int x, int y, int width, int height, int colour) {
  gdImageFilledRectangle(img, x + offset_x, y + offset_y, x + offset_x + width - 1, y + offset_y + height - 1, colour);
}

void Canvas::fillRectAbs(int x1, int y1, int x2, int y2) {
  gdImageFilledRectangle(img, x1 + offset_x, y1 + offset_y, x2 + offset_x, y2 + offset_y, colour);
}

void Canvas::fillRectAbs(int x1, int y1, int x2, int y2, int colour) {
  gdImageFilledRectangle(img, x1 + offset_x, y1 + offset_y, x2 + offset_x, y2 + offset_y, colour);
}

void Canvas::drawText(std::string s, int x, int y, VAlignment v, HAlignment h, double angle) {
  drawText((char*)s.c_str(), x, y, v, h, colour, fontsize, angle);
}

void Canvas::drawText(std::string s, int x, int y, VAlignment v, HAlignment h, int colour, int size, double angle) {
  drawText((char*)s.c_str(), x, y, v, h, colour, size, angle);
}

void Canvas::drawText(char* s, int x, int y, VAlignment v, HAlignment h, double angle) {
  drawText(s, x, y, v, h, colour, fontsize, angle);
}

enum BRECT_POS { LOWERLEFT_X, LOWERLEFT_Y, LOWERRIGHT_X, LOWERRIGHT_Y, UPPERRIGHT_X, UPPERRIGHT_Y, UPPERLEFT_X, UPPERLEFT_Y };
void Canvas::drawText(char* s, int x, int y, VAlignment v, HAlignment h, int colour, int size, double angle) {
  int brect[8];
  gdImageStringFT(NULL, brect, colour, (char*)fontfile.c_str(), size, angle, 0, 0, s);

  switch(h) {
  case LEFT:
    x -= brect[UPPERLEFT_X];
    break;
    
  case CENTER:
    x -= (brect[LOWERRIGHT_X] + brect[UPPERLEFT_X]) / 2;
    break;
    
  case RIGHT:
    x -= brect[LOWERRIGHT_X];
    break;
  }
  
  switch(v) {
  case TOP:
    y -= brect[UPPERRIGHT_Y];
    break;
    
  case MIDDLE:
    y -= (brect[LOWERLEFT_Y] + brect[UPPERRIGHT_Y]) / 2;
    break;
    
  case BOTTOM:
    y -= brect[LOWERLEFT_Y];
    break;
  }
  x += offset_x;
  y += offset_y;

  gdImageStringFT(img, brect, colour, (char*)fontfile.c_str(), size, angle, x, y, s);
} 

void Canvas::setOffsets(int offset_x, int offset_y) {
  this->offset_x = offset_x;
  this->offset_y = offset_y;
}

void addValue(gdPoint* p, int x, int y, int numpoints) {
  for(int i = 0; i < numpoints; i++) {
    p[i].x += x;
    p[i].y += y;
  }
}

void Canvas::drawSymbol(enum SYMBOL s, int x, int y) {
  gdPoint utri[] = { {5, 1}, {0, 9}, {10, 9} };
  gdPoint ltri[] = { { 1, 0}, {1, 10}, {9, 5} };
  gdPoint rtri[] = { { 9, 0}, {9, 10}, {1, 5} };
  gdPoint dtri[] = { { 5, 9}, {0, 1}, {10, 1} };
  gdPoint dmond[] = { { 5, 1}, {1, 5}, {5, 9}, {9, 5} };
  gdPoint star6[] = { {8, 5}, {10, 2}, {7, 2}, {5, -1}, {4, 2}, {0, 2}, {2, 5}, {0, 9}, {3, 8}, {5, 11}, {7, 8}, {10, 9} };
  switch(s) {
  case UTRIANGLE:
    addValue(utri, x, y, 3);
    gdImagePolygon(img, utri, 3, colour);
    break;
  case LTRIANGLE:
    addValue(ltri, x, y, 3);
    gdImagePolygon(img, ltri, 3, colour);
    break;
  case RTRIANGLE:
    addValue(rtri, x, y, 3);
    gdImagePolygon(img, rtri, 3, colour);
    break;
  case DTRIANGLE:
    addValue(dtri, x, y, 3);
    gdImagePolygon(img, dtri, 3, colour);
    break;
  case DIAMOND:
    addValue(dmond, x, y, 4);
    gdImagePolygon(img, dmond, 4, colour);
    break;
  case STAR6:
    addValue(star6, x, y, 12);
    gdImagePolygon(img, star6, 12, colour);
    break;
  case SQUARE:
    gdImageRectangle(img, x + 1, y + 1, x + 9, y + 9, colour);
    break;
  case CIRCLE:
    gdImageArc(img, x + 5, y + 5, 9, 9, 0, 360, colour);
    break;
  case NONE:
    break;
  }
}

void Canvas::fillSymbol(enum SYMBOL s, int x, int y) {
  gdPoint utri[] = { {5, 1}, {0, 9}, {10, 9} };
  gdPoint ltri[] = { { 1, 0}, {1, 10}, {9, 5} };
  gdPoint rtri[] = { { 9, 0}, {9, 10}, {1, 5} };
  gdPoint dtri[] = { { 5, 9}, {0, 1}, {10, 1} };
  gdPoint dmond[] = { { 5, 1}, {1, 5}, {5, 9}, {9, 5} };
  gdPoint star6[] = { {8, 5}, {10, 2}, {7, 2}, {5, -1}, {4, 2}, {0, 2}, {2, 5}, {0, 9}, {3, 8}, {5, 11}, {7, 8}, {10, 9} };
  switch(s) {
  case UTRIANGLE:
    addValue(utri, x, y, 3);
    gdImageFilledPolygon(img, utri, 3, colour);
    break;
  case LTRIANGLE:
    addValue(ltri, x, y, 3);
    gdImageFilledPolygon(img, ltri, 3, colour);
    break;
  case RTRIANGLE:
    addValue(rtri, x, y, 3);
    gdImageFilledPolygon(img, rtri, 3, colour);
    break;
  case DTRIANGLE:
    addValue(dtri, x, y, 3);
    gdImageFilledPolygon(img, dtri, 3, colour);
    break;
  case DIAMOND:
    addValue(dmond, x, y, 4);
    gdImageFilledPolygon(img, dmond, 4, colour);
    break;
  case STAR6:
    addValue(star6, x, y, 12);
    gdImageFilledPolygon(img, star6, 12, colour);
    break;
  case SQUARE:
    gdImageFilledRectangle(img, x + 1, y + 1, x + 9, y + 9, colour);
    break;
  case CIRCLE:
    gdImageFilledEllipse(img, x + 5, y + 5, 9, 9, colour);
    break;
  case NONE:
    break;
  }
}
