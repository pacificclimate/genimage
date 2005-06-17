#include "point.h"

float clip_precision(float number, int num_mantissa_bits) {
  number = roundf(number * (1 << num_mantissa_bits)) / (1 << num_mantissa_bits);
  return number;
}

double clip_precision(double number, int num_mantissa_bits) {
  number = round(number * (1 << num_mantissa_bits)) / (1 << num_mantissa_bits);
  return number;
}

float clip(float number) {
  return clip_precision(number, FLOAT_PRECISION);
}

double clip(double number) {
  return clip_precision(number, DOUBLE_PRECISION);
}

int clip(int number) {
  return number;
}
