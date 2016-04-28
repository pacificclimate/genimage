#ifndef __GENIMAGE_LEGENDTOKEN_H
#define __GENIMAGE_LEGENDTOKEN_H

enum SYMBOL {UTRIANGLE, LTRIANGLE, RTRIANGLE, DTRIANGLE, DIAMOND, STAR6, SQUARE, CIRCLE, NONE};

class LegendToken {
public:
  LegendToken(std::string name, int colour, enum SYMBOL sym, bool filled) {
    this->name = name;
    this->colour = colour;
    this->sym = sym;
    this->filled = filled;
  }
  std::string name;
  int colour;
  enum SYMBOL sym;
  bool filled;
};

#endif
