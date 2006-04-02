#include <iostream>
#include <string>
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

Canvas::Canvas(string fontfile, int width, int height, gdImagePtr img) {
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

void Canvas::drawText(string s, int x, int y, VAlignment v, HAlignment h) {
  drawText((char*)s.c_str(), x, y, v, h, colour, fontsize);
}

void Canvas::drawText(string s, int x, int y, VAlignment v, HAlignment h, int colour, int size) {
  drawText((char*)s.c_str(), x, y, v, h, colour, size);
}

void Canvas::drawText(char* s, int x, int y, VAlignment v, HAlignment h) {
  drawText(s, x, y, v, h, colour, fontsize);
}

void Canvas::drawText(char* s, int x, int y, VAlignment v, HAlignment h, int colour, int size) {
  int brect[8];
  gdImageStringFT(NULL, brect, colour, (char*)fontfile.c_str(), size, 0, 0, 0, s);

  switch(h) {
  case LEFT:
    x -= brect[0];
    break;
    
  case CENTER:
    x -= (brect[2] + brect[0]) / 2;
    break;
    
  case RIGHT:
    x -= brect[2];
    break;
  }
  
  switch(v) {
  case TOP:
    y -= brect[5];
    break;
    
  case MIDDLE:
    y -= (brect[1] + brect[5]) / 2;
    break;
    
  case BOTTOM:
    y -= brect[1];
    break;
  }
  x += offset_x;
  y += offset_y;

  gdImageStringFT(img, brect, colour, (char*)fontfile.c_str(), size, 0, x, y, s);
} 

void Canvas::setOffsets(int offset_x, int offset_y) {
  this->offset_x = offset_x;
  this->offset_y = offset_y;
}
