#ifndef __GENIMAGE_GENIMAGE_H
#define __GENIMAGE_GENIMAGE_H
#include <string>
#include <iostream>
#include <stdio.h>
#include <gd.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/file.h>
#include <errno.h>
using namespace std;

#include "range.h"
#include "legend.h"
#include "point.h"
#include "support.h"

#define VALUE_LENGTH 13
#define LL_LENGTH 10
#define SLMASK_LENGTH 2
#define NUM_TIMESLICES 17
#define POINT_SIZE 4
#define LINE_WIDTH 3

#define STATIC_ARGS 29
#define LEG_BOTTOM_TEXT_HEIGHT 55
#define LEG_TOP_TEXT_HEIGHT 35
#define LEG_HEIGHT 120
#define LAT_WIDTH 40
#define LAT_EXTRA_HEIGHT 10
#define LEG_EXTRA_RWIDTH 20
#define BORDER_WIDTH 1

#define IMAGE_TYPE_PNG 0
#define IMAGE_TYPE_JPG 1

// Earth's radius in km
#define EARTH_RADIUS 6378.1

#define JPEG_QUALITY 90

enum PLOT_TYPES{ TYPE_MAP, TYPE_LEGEND, TYPE_LAT, TYPE_ALL, TYPE_TEXT, TYPE_MASK, TYPE_GEOREF, TYPE_PLOTINFO, TYPE_REGIONONLY, TYPE_INVALID };

inline int lontox(int max_x, float off_lon, float min_lon, float max_lon, float lon) {
  float difflon = max_lon - min_lon;

  return (int)roundf(((lon - min_lon - off_lon) / difflon) * max_x);
}

inline int lattoy(int max_y, float off_lat, float min_lat, float max_lat, float lat) {
  float difflat = max_lat - min_lat;

  return (int)roundf(max_y - (((lat - min_lat - off_lat) / difflat) * max_y));
}

inline void hline(int** img, int x1, int x2, int y, int width, int colour) {
  int* imgptr = img[y];
  for(; x1 < x2; x1++) {
    imgptr[x1] = colour;
  }
}

inline void vline(int** img, int x, int y1, int y2, int width, int colour) {
  for(; y1 < y2; y1++) {
    img[y1][x] = colour;
  }
}

inline void seektoword(FILE * in) {
  char c;
  while(!isspace(c = fgetc(in)));
  while(isspace(c = fgetc(in)));
  ungetc(c, in);
}

inline char* seektoword(char * in) {
    while(!isspace(*(in++)) && *in);
    while(isspace(*(in++)) && *in);
    in--;
    return in;
}

inline int readint(FILE * in, char * buf) {
    char * ptr = buf;
    while(!isspace(*(ptr++) = fgetc(in)));
    *ptr = '\0';
    return atoi(buf);
}

inline int readint(char *in, char * buf) {
    char * ptr = buf;
    while(!isspace(*(ptr++) = *(in++)) && *in);
    *ptr = '\0';
    return atoi(buf);
}

inline void seektonextline(FILE * in) {
  while('\n' != fgetc(in));
}

inline int load_grid(FILE *infile, int rows, int cols, char* dataline, int dl_length, float* data, int data_width) {
  char* ptr;
  for(int y = 0; y < rows; y++) {
    if(fgets(dataline, dl_length, infile)) {
      ptr = dataline + 1;
      for(int x = 0; x < cols; x++) {
	*(ptr + (data_width - 1)) = '\0';
	data[y * cols + x] = atof(ptr);
	ptr += data_width;
      }	  
    }
  }
  return 0;
}

inline int load_grid(FILE *infile, int rows, int cols, char* dataline, int dl_length, int* data, int data_width) {
  char* ptr;
  for(int y = 0; y < rows; y++) {
    if(fgets(dataline, dl_length, infile)) {
      ptr = dataline + 1;
      for(int x = 0; x < cols; x++) {
	*(ptr + (data_width - 1)) = '\0';
	data[y * cols + x] = atoi(ptr);
	ptr += data_width;
      }	  
    }
  }
  return 0;
}

inline int load_grid(FILE *infile, int rows, int cols, char* dataline, int dl_length, float* data, int data_width, struct range* drange) {
  char* ptr;
  for(int y = 0; y < rows; y++) {
    if(fgets(dataline, dl_length, infile)) {
      ptr = dataline + 1;
      for(int x = 0; x < cols; x++) {
	*(ptr + (data_width - 1)) = '\0';
	data[y * cols + x] = atof(ptr);
	drange->add(data[y * cols + x]);
	ptr += data_width;
      }	  
    }
  }
  return 0;
}

inline int load_grid(FILE *infile, int rows, int cols, char* dataline, int dl_length, float* data, int data_width, struct range* drange, int* slmask, int allow_value) {
  char* ptr;
  for(int y = 0; y < rows; y++) {
    if(fgets(dataline, dl_length, infile)) {
      ptr = dataline + 1;
      for(int x = 0; x < cols; x++) {
	*(ptr + (data_width - 1)) = '\0';
	data[y * cols + x] = atof(ptr);
	if(slmask[y * cols + x] == allow_value) {
	  drange->add(data[y * cols + x]);
	}
	ptr += data_width;
      }	  
    }
  }
  return 0;
}

inline int parse_data_point(char* datapt, Point<float>& pt) {
  if(datapt && strchr(datapt, ':')) {
    char* nxt = strchr(datapt, ':');
    if(nxt && *nxt) {
      *nxt = '\0';
      nxt++;
    } else {
      return 0;
    }

    pt.x = atof(datapt);
    datapt = nxt;

    nxt = strchr(datapt, ':');
    if(nxt && *nxt) {
      *nxt = '\0';
      nxt++;
    } else {
      return 0;
    }
    pt.y = atof(datapt);
    pt.selected = atoi(nxt);
  } else {
    return 0;
  }
  return 1;
}

inline Point<float>* parse_data_point(char* datapt) {
  Point<float>* dpoint = new Point<float>();
  if(parse_data_point(datapt, *dpoint)) {
    return dpoint;
  } else {
    delete dpoint;
    return NULL;
  }
}

#endif
