const static int continuous[] = {
  0x00006FFF, 0x00007FFF, 0x00008FFF, 0x00009FFF, 0x0000AFFF, 0x0000BFFF, 0x0000CFFF, 0x0000DFFF, 0x0000EFFF, 0x0000FFFF, 0x000FFFFF, 0x001FFFEF, 0x002FFFDF, 0x003FFFCF, 0x004FFFBF, 0x005FFFAF, 0x006FFF9F, 0x007FFF8F, 0x008FFF7F, 0x009FFF6F, 0x00AFFF5F, 0x00BFFF4F, 0x00CFFF3F, 0x00DFFF2F, 0x00EFFF1F, 0x00FFFF0F, 0x00FFFF00, 0x00FFEF00, 0x00FFDF00, 0x00FFCF00, 0x00FFBF00, 0x00FFAF00, 0x00FF9F00, 0x00FF8F00, 0x00FF7F00, 0x00FF6F00, 0x00FF5F00, 0x00FF4F00, 0x00FF3F00, 0x00FF2F00, 0x00FF1F00, 0x00FF0F00, 0x00FF0000, 0x00EF0000, 0x00DF0000, 0x00CF0000, 0x00BF0000, 0x00AF0000, 0x009F0000, 0x008F0000
};

const static int continuous_rev[] = {
  0x008F0000, 0x009F0000, 0x00AF0000, 0x00BF0000, 0x00CF0000, 0x00DF0000, 0x00EF0000, 0x00FF0000, 0x00FF0F00, 0x00FF1F00, 0x00FF2F00, 0x00FF3F00, 0x00FF4F00, 0x00FF5F00, 0x00FF6F00, 0x00FF7F00, 0x00FF8F00, 0x00FF9F00, 0x00FFAF00, 0x00FFBF00, 0x00FFCF00, 0x00FFDF00, 0x00FFEF00, 0x00FFFF00, 0x00FFFF0F, 0x00EFFF1F, 0x00DFFF2F, 0x00CFFF3F, 0x00BFFF4F, 0x00AFFF5F, 0x009FFF6F, 0x008FFF7F, 0x007FFF8F, 0x006FFF9F, 0x005FFFAF, 0x004FFFBF, 0x003FFFCF, 0x002FFFDF, 0x001FFFEF, 0x000FFFFF, 0x0000FFFF, 0x0000EFFF, 0x0000DFFF, 0x0000CFFF, 0x0000BFFF, 0x0000AFFF, 0x00009FFF, 0x00008FFF, 0x00007FFF, 0x00006FFF
};

// Stub
const static int stepwise_rev[] = {
  0x00FF0000, 0x00FF5F2A, 0x00FF9138, 0x00FFCD05, 0x00FFFF6A, 0x00D6FF74, 0x00CECD00, 0x0000FF00, 0x0000FEFC, 0x0000807F, 0x00007FFF, 0x000000FF
};

// Stub
const static int stepwise[] = {
  0x000000FF, 0x00007FFF, 0x0000807F, 0x0000FEFC, 0x0000FF00, 0x00CECD00, 0x00D6FF74, 0x00FFFF6A, 0x00FFCD05, 0x00FF9138, 0x00FF5F2A, 0x00FF0000
};


/*
Names Trevor has associated with colours:
Red:           0x00FF0000
Redderyellow:  0x00FF5F2A
Redyellow:     0x00FF9138
Redyellower:   0x00FFCD05
Yellow:        0x00FFFF9A
Yellowergreen: 0x00FFFF6A
Yellowgreen:   0x00D6FF74
Yellowgreener: 0x00CECD00
Green, Land:   0x0000FF00
Bluegreener:   0x0000FEFC
Bluegreen:     0x0000807F
Bluergreen:    0x00007FFF
Blue, Ocean:   0x000000FF
Black:         0x00000000
White:         0x00000000
*/
