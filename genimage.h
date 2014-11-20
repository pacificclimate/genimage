#ifndef __GENIMAGE_GENIMAGE_H
#define __GENIMAGE_GENIMAGE_H
#include <string>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <gd.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/file.h>
#include <errno.h>

#include "range.h"
#include "legend.h"
#include "point.h"
#include "support.h"


#define LTRED 0x00FF0000
#define LTORANGE 0x00FF6600
#define LTYELLOW 0x00FFFF00
#define LTYELLOWGRN 0x00CC00FF
#define DKYELLOWGRN 0x00660080
#define DKGREENYLW 0x0066FF00
#define LTGREENYLW 0x0066FF00
#define BRTGREEN 0x0080FF80
#define LTGREEN 0x0000FF00
#define DKGREEN 0x00008000
#define DKRGREEN 0x00004000
#define BRTBLUEGRN 0x0080FFCC
#define LTBLUEGRN 0x0000FF66
#define DKBLUEGRN 0x00008033
#define DKRBLUEGRN 0x00004018
#define LTGCYAN 0x0000FFCC
#define DKGCYAN 0x00008066
#define DKRGCYAN 0x00004033
#define LTBCYAN 0x0000CCFF
#define LTCBLUE 0x000066FF
#define LTBLUE 0x000000FF
#define LTPURPLE 0x006600FF
#define DKPURPLE 0x00330080
#define LTPMAGENTA 0x00CC00FF
#define LTMAGENTA 0x00FF00CC
#define LTRMAGENTA 0x00FF0080
#define DKRMAGENTA 0x00800040
#define DKRRMAGENTA 0x00400020

#define WHITE 0x00FFFFFF
#define LTGRAY 0x00C0C0C0
#define DKGRAY 0x00808080
#define DKRGRAY 0x00404040
#define BLACK 0x00000000

// MPN: Semi-transparent gray
#define TGRAY 0x60707070

#define VALUE_LENGTH 13
#define LL_LENGTH 10
#define SLMASK_LENGTH 2
#define NUM_TIMESLICES 17
#define GCMINFO_COLS 45

#define STATIC_ARGS 29

#define IMAGE_TYPE_PNG 0
#define IMAGE_TYPE_JPG 1

// Earth's radius in km
#define EARTH_RADIUS 6378.1

#define JPEG_QUALITY 90

enum ANOM_TYPE{DEFAULT=0, ANOMALY, ABSOLUTE};

enum PLOT_TYPES {
  TYPE_MAP = 0, // Use this to plot a map
  TYPE_TEXT = 1,
  TYPE_MASK = 2,
  TYPE_GEOREF = 3,
  TYPE_PLOTINFO = 4,  // unused
  TYPE_REGIONONLY = 5, // image shown on 'region' tab
  TYPE_SCATTER_TIMESLICE = 6, // Use this to plot a single var scatter
  TYPE_SCATTER_TIMESLICE_HIST = 7,
  TYPE_SCATTER_VARIABLE = 8,
  TYPE_SCATTER_TIMESLICE_TEXT = 9,
  TYPE_SCATTER_VARIABLE_TEXT = 10,
  TYPE_MAP_DIFFERENCE = 11, // Map plot of change with respect to timeslice2.
  TYPE_SCENARIO_DATA = 12, // Hackish. Spits out scenarios format data (.dat file)
  TYPE_SLMASK_DATA = 13, // / and so on for SLMASK, LAT,S LONGS
  TYPE_LATS_DATA = 14,
  TYPE_LONGS_DATA = 15,
  TYPE_BOXPLOT_TIMESLICE = 16,
  TYPE_BOXPLOT_TIMESLICE_TEXT = 17,
  TYPE_STICKPLOT = 18,
  TYPE_GEOTIFF = 19,
  TYPE_BANDS_TIMESLICE = 20,
  TYPE_BANDS_TIMESLICE_HIST = 21,
  TYPE_SCENARIO_SET_METADATA = 22,
  TYPE_INVALID };

inline int nearest_offset(int size, double* data, double value) {
  double old_diff = INFINITY;
  int i;
  for(i = 0; i < size; i++) {
    const double diff = fabs(data[i] - value);
    if(old_diff < diff) {
      break;
    }
    old_diff = diff;
  }
  return i - 1;
}

inline double nearest(int size, double* data, double value) {
  return data[nearest_offset(size, data, value)];
}

inline int lontox(int max_x, double off_lon, double min_lon, double max_lon, double lon) {
  double difflon = max_lon - min_lon;

  return (int)roundf(((lon - min_lon - off_lon) / difflon) * max_x);
}

inline int lattoy(int max_y, double off_lat, double min_lat, double max_lat, double lat) {
  double difflat = max_lat - min_lat;

  return (int)roundf(max_y - (((lat - min_lat - off_lat) / difflat) * max_y));
}

inline void hline(int** img, int x1, int x2, int y, int /*width*/, int colour) {
  int* imgptr = img[y];
  for(; x1 < x2; x1++) {
    imgptr[x1] = colour;
  }
}

inline void vline(int** img, int x, int y1, int y2, int /*width*/, int colour) {
  for(; y1 < y2; y1++) {
    img[y1][x] = colour;
  }
}

inline void seektoword(FILE * in) {
  char c;
  while(!isspace(c = fgetc(in))) ;
  while(isspace(c = fgetc(in))) ;
  ungetc(c, in);
}

inline char* seektoword(char * in) {
    while(!isspace(*(in++)) && *in) ;
    while(isspace(*(in++)) && *in) ;
    in--;
    return in;
}

inline int readint(FILE * in, char * buf) {
    char * ptr = buf;
    while(!isspace(*(ptr++) = fgetc(in))) ;
    *ptr = '\0';
    return atoi(buf);
}

inline int readint(char *in, char * buf) {
    char * ptr = buf;
    while(!isspace(*(ptr++) = *(in++)) && *in) ;
    *ptr = '\0';
    return atoi(buf);
}

inline void seektonextline(FILE * in) {
  while('\n' != fgetc(in)) ;
}

inline int load_grid(FILE *infile, int rows, int cols, char* dataline, int dl_length, double* data, int data_width) {
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

inline int load_grid(FILE *infile, int rows, int cols, char* dataline, int dl_length, double* data, int data_width, Range* drange) {
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

inline int load_grid(FILE *infile, int rows, int cols, char* dataline, int dl_length, double* data, int data_width, Range* drange, int* slmask, int allow_value) {
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

inline int parse_data_point(const char* datapt, Point& pt) {
  // "lon:lat:selected"
  if(datapt && strchr(datapt, ':')) {
    char* const dp_thing = new char[strlen(datapt) + 1];
    char* dp = strcpy(dp_thing, datapt);
    char* nxt = strchr(dp, ':');
    if(nxt && *nxt) {
      *nxt = '\0';
      nxt++;
    } else {
      return 0;
    }

    pt.x = atof(dp);
    dp = nxt;

    nxt = strchr(dp, ':');
    if(nxt && *nxt) {
      *nxt = '\0';
      nxt++;
    } else {
      return 0;
    }
    pt.y = atof(dp);
    pt.selected = atoi(nxt);
    delete[] dp_thing;
  } else {
    return 0;
  }
  return 1;
}

inline Point parse_data_point(const char* datapt) {
  Point dpoint;
  parse_data_point(datapt, dpoint);
  return dpoint;
}

#endif
