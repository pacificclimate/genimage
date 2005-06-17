#ifndef __POINT_H
#define __POINT_H
#include <iostream>

#ifndef INFINITY
#define INFINITY infinityf()
#endif

using namespace std;

#define FLOAT_PRECISION 18
#define DOUBLE_PRECISION 48

float clip_precision(float number, int num_mantissa_bits);
double clip_precision(double number, int num_mantissa_bits);

int clip(int number);
float clip(float number);
double clip(double number);

template <class T>
class Point {
 public:
  // Constructors
  Point() { x = 0; y = 0; selected = false; l_used = r_used = false; }
  Point(const T x, const T y, bool sel = false) { this->x = x; this->y = y; selected = sel; }
  Point(const Point<T>& in) { x = in.x; y = in.y; selected = in.selected; }

  // Misc
  void clip() { x = ::clip(x); y = ::clip(y); }
  float dist_from(const Point<float>& p) { return hypotf(x - p.x, y - p.y); }
  double dist_from(const Point<double>& p) { return hypot(x - p.x, y - p.y); }
  int dist_from(const Point<int>& p) { return hypotl(x - p.x, y - p.y); }
  T cross(const Point<T>& p) { T temp; temp = (x * p.y) - (p.x * y); return temp; }

  // Operators
  Point<T> operator - (const Point<T>& p) const { Point<T> temp; temp.x = x - p.x; temp.y = y - p.y; return temp; }
  Point<T> operator + (const Point<T>& p) const { Point<T> temp; temp.x = x + p.x; temp.y = y + p.y; return temp; }
  T operator * (const Point<T>& p) const { return x*p.x +  y*p.y; }
  Point<T> operator * (const T& p) const { return Point<T>(x*p + y*p.y); }
  Point<T> operator / (const Point<T>& p) const { return Point<T>(x/p.x, y/p.y); }
  Point<T> operator / (const T& p) const { return Point<T>(x/p, y/p); }
  Point<T>& operator -= (const Point<T>& p) { x -= p.x; y -= p.y; return *this; }
  Point<T>& operator += (const Point<T>& p) { x += p.x; y += p.y; return *this; }
  bool operator == (const Point<T>& p) const { if(p.x == x && p.y == y) { return true; } else { return false; } }
  bool operator != (const Point<T>& p) const { if(p.x != x || p.y != y) { return true; } else { return false; } }
  bool operator <= (const Point<T>& p) const { if((p.x*p.x + p.y*p.y) >= (x*x + y*y)) { return true; } else { return false; } }
  bool operator < (const Point<T>& p) const { if((p.x*p.x + p.y*p.y) > (x*x + y*y)) { return true; } else { return false; } }
  bool operator >= (const Point<T>& p) const { if((p.x*p.x + p.y*p.y) <= (x*x + y*y)) { return true; } else { return false; } }
  bool operator > (const Point<T>& p) const { if((p.x*p.x + p.y*p.y) < (x*x + y*y)) { return true; } else { return false; } }

  friend ostream& operator<< <> (ostream& os, const Point<T>& f);
  T x, y;
  bool l_used; // Left edge used
  bool r_used; // Right edge used
  bool selected;
};


template <class T>
Point<T> clip(const Point<T>& p) {
  Point<T> temp(clip(p.x), clip(p.y));;
  return temp;
}

template <class T>
ostream& 
operator<<(ostream& os, const Point<T>& f) {
  os << "(" << f.x << "," << f.y << ")";
  return os;
}

#endif
