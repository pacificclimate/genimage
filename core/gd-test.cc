#include "genimage.h"
#include "ConfigFile.h"

class Canvas {
public:
  enum VAlignment { TOP, MIDDLE, BOTTOM };
  enum HAlignment { LEFT, CENTER, RIGHT };
  Canvas(int width, int height, gdImagePtr img = 0) {
    this->width = width;
    this->height = height;
    offset_x = 0;
    offset_y = 0;
    hasParent = true;
    if(!img) {
      this->img = gdImageCreateTrueColor(width, height);
      hasParent = false;
    }
  }

  ~Canvas() {
    if(!hasParent) {
      gdImageDestroy(img);
    }
  }

  bool writeImage(char* filename) {
    FILE* f = fopen(filename, "wb");
    if(!f) {
      return false;
    }
    gdImagePng(img, f);
    fclose(f);
    return true;
  }

  bool writeImage(FILE* f) {
    if(!f) {
      return false;
    }
    gdImagePng(img, f);
    return true;
  }

  void setLineThickness(int width) {
    gdImageSetThickness(img, width);
  }

  void setAntiAliased() {
    gdImageSetAntiAliased(img, colour);
  }

  void setAntiAliased(int colour) {
    gdImageSetAntiAliased(img, colour);
  }

  void drawLine(int x1, int y1, int x2, int y2) {
    gdImageLine(img, x1 + offset_x, y1 + offset_y, x2 + offset_x, y2 + offset_y, colour);
  }

  void drawLine(int x1, int y1, int x2, int y2, int colour) {
    gdImageLine(img, x1 + offset_x, y1 + offset_y, x2 + offset_x, y2 + offset_y, colour);
  }

  void copy(gdImagePtr src, int dx, int dy, int sx, int sy, int width, int height) {
    gdImageCopy(img, src, dx + offset_x, dy + offset_y, sx, sy, width, height);
  }

  void fillRect(int x, int y, int width, int height) {
    gdImageFilledRectangle(img, x + offset_x, y + offset_y, width, height, colour);
  }

  void fillRect(int x, int y, int width, int height, int colour) {
    gdImageFilledRectangle(img, x + offset_x, y + offset_y, width, height, colour);
  }

  void drawText(char* s, int x, int y, VAlignment v, HAlignment h) {
    drawText(s, x, y, v, h, colour, fontsize);
  }

  void drawText(char* s, int x, int y, VAlignment v, HAlignment h, int colour, int size) {
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

int main(int argc, char** argv) {
  Canvas c(400, 100);
  char* string1 = "This Corrosion joke";
  char* string2 = "Bool Sheet";
  char* string3 = "Math Emetics";

  FILE* f = fopen("/home/bronaugh/work/cics/scen-basemaps/cgcm2_1_0_cdn.png", "rb");
  gdImagePtr p = gdImageCreateFromPng(f);
  fclose(f);

  printf("%x\n", p->pixels);
  printf("%x\n", p->pixels[0]);
  printf("%x\n", p->pixels[0][0]);

  if(argc < 2) {
    printf("Usage: ./gd-test <config file>\n");
  }

  ConfigFile cf(argv[1]);
  cf.readInto(c.fontfile, "fontfile");

  c.colour = 0x00FFFFFF;
  c.fontsize = 12;

  c.drawText(string1, 0, 0, Canvas::TOP, Canvas::LEFT);
  c.drawText(string2, 200, 50, Canvas::MIDDLE, Canvas::CENTER);
  c.drawText(string3, 400, 100, Canvas::BOTTOM, Canvas::RIGHT);

  gdImagePng(c.img, stdout);
}
