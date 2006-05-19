#include "genimage.h"
#include "canvas.h"
#include "displayer.h"
#include "datamanager.h"
#include "config.h"
#include "ConfigFile.h"

#include <map>
#include <algorithm>
#include <boost/tokenizer.hpp>

/* Equatorial Cylindrical Equidistant Projection
    * lat = ((-y + center_y) / height) * 180
    * lon = ((x - center_x) / width) * 360
    * x = ((lon / 360) * width) + center_x
    * y = ((-lat / 180) * height) + center_y
*/


// NOTE: CONVERT TO auto_ptr ET AL: http://ootips.org/yonat/4dev/smart-pointers.html
using namespace std;

void parseArgs(Config& c, Displayer& disp, DataManager& dm, int argc, char** argv) {
  char val;
  Point* temp;
  list<Point* > lpoints;
  char* optstring = "a:b:c:e:f:g:h:i:j:k:l:m:n:o:p:r:s:t:u:v:w:x:y:z:";
  static struct option options[] = {
    { "xrange-min", 1, 0, 'a'},
    { "xvariable", 1, 0, 'b'},
    { "colour-map", 1, 0, 'c'},
    { "yvariable", 1, 0, 'd'},
    { "xrange-max", 1, 0, 'e'},
    { "identify-text", 1, 0, 'f'},
    { "ocean-plot", 1, 0, 'g'},
    { "scatter-type", 1, 0, 'h'},
    { "poly-point", 1, 0, 'j'},
    { "y-axis-text", 1, 0, 'k'},
    { "x-axis-text", 1, 0, 'l'},
    { "model", 1, 0, 'm'},
    { "region", 1, 0, 'n'},
    { "output-file", 1, 0, 'o'},
    { "plot-type", 1, 0, 'p'},
    { "scenario-set", 1, 0, 'q'},
    { "resolution", 1, 0, 'r'},
    { "expt", 1, 0, 's'},
    { "timeslice", 1, 0, 't'},
    { "yrange-min", 1, 0, 'u'},
    { "yrange-max", 1, 0, 'v'},
    { "data-point", 1, 0, 'x'},
    { "timeofyear", 1, 0, 'y'},
    { "zoom-factor", 1, 0, 'z'},
    { "show-grid", 0, &disp.grid, 1},
    { "colour-map-inverse", 0, &disp.colour_map_rev, 1},
    { "dynamic-range", 0, &disp.range_dynamic, 1},
    { "percentiles", 0, &c.percentiles, 1},
    { "percent-change-calculations", 0, &c.pct_change_data, 1},
    { 0, 0, 0, 0 }
  };
  while((val = getopt_long(argc, argv, optstring, options, NULL)) != -1) {
    switch(val) {
    case 0:
      break;
    case 'a':
      // xrange-min
      disp.xrange_min = atof(optarg);
      break;
    case 'b':
      // X variable
      c.xvariable = optarg;
      break;
    case 'c':
      // colour-map
      disp.colour_map = atoi(optarg);
      break;
    case 'd':
      // Y variable
      c.yvariable = optarg;
      break;
    case 'e':
      // xrange-max
      disp.xrange_max = atof(optarg);
      break;
    case 'f':
      // identify-text
      disp.identify_text = optarg;
      break;
    case 'g':
      // ocean-plot
      c.plot_over_ocean = atoi(optarg);
      break;
    case 'h':
      // scatter-type
      c.scatter_type = atoi(optarg);
      break;
    case 'j':
      // poly-point
      temp = parse_data_point(optarg);
      if(temp) {
	lpoints.push_back(temp);
      }
      break;
    case 'k':
      // Y axis text
      disp.yaxis_text = optarg;
      break;
    case 'l':
      // X axis text
      disp.xaxis_text = optarg;
      break;
    case 'm':
      // model
      dm.model = optarg;
      break;
    case 'n':
      // region
      dm.region = optarg;
      break;
    case 'o':
      // output-file
      c.outfile = optarg;
      break;
    case 'p':
      // plot-type
      c.plot_type = atoi(optarg);
      break;
    case 'q':
      // scenario-set
      c.scenario_set = optarg;
      break;
    case 'r':
      // resolution
      c.resolution = atoi(optarg);
      break;
    case 's':
      // expt
      dm.expt = optarg;
      break;
    case 't':
      // timeslice
      dm.timeslice = optarg;
      break;
    case 'u':
      // yrange-min
      disp.yrange_min = atof(optarg);
      break;
    case 'v':
      // yrange-max
      disp.yrange_max = atof(optarg);
      break;
    case 'w':
      // plot-width
      disp.plot_width = atoi(optarg);
      break;
    case 'x':
      // data-point
      if(parse_data_point(optarg, c.dpoint)) {
	c.point_present = true;
      }
      break;
    case 'y':
      // timeofyear
      dm.timeofyear = atoi(optarg);
      break;
    case 'z':
      // zoom-factor
      c.zoom_factor = atoi(optarg);
      break;
    }
  }

  // Store up the data points
  c.numpoints = lpoints.size();
  c.points = 0;
  if(c.numpoints) {
    list<Point* >::const_iterator lpoints_iter;
    int i = 0;
    c.points = new Point*[c.numpoints];
    for(lpoints_iter = lpoints.begin(); lpoints_iter != lpoints.end(); lpoints_iter++) {
      c.points[i] = (*lpoints_iter);
      i++;
    }
  }
}

void handleMap(Displayer& disp, DataManager& dm) {
  dm.variable = dm.config.xvariable;
  dm.open_datafile();

  int datasize = dm.data_size();
  double* data = new double[datasize];
  double* lats = new double[dm.lats_size()];
  double* longs = new double[dm.longs_size()];
  double* grid_lats = new double[dm.gridlats_size()];
  double* grid_longs = new double[dm.gridlongs_size()];
  int* slmask = new int[datasize];
  int* drawmask = new int[datasize];
  int* mask = new int[datasize];
  gdImagePtr basemap = dm.get_basemap();

  dm.get_slmask(slmask);
  dm.get_lats(lats);
  dm.get_longs(longs);
  dm.get_gridlats(grid_lats, lats);
  dm.get_gridlongs(grid_longs, longs);
  dm.get_drawmask(drawmask, slmask);
  dm.get_datamask(mask, slmask, grid_lats, grid_longs);
  dm.get_data(data);

  Range datarange(data, datasize);
  Range xrange(grid_longs, dm.gridlongs_size());
  Range yrange(grid_lats, dm.gridlats_size());
  Legend* legend = disp.getLegend(datarange);

  disp.setOffsets(basemap);
  disp.createCanvas(dm.config.fontfile);
  disp.drawMap(basemap, *legend, drawmask, mask, data, grid_lats, grid_longs, dm.numrows(), dm.numcols());
  disp.drawTicks(xrange, yrange);
  disp.drawScale(*legend);
  disp.clearIdentifyArea();
  disp.drawIdentifyText();
  disp.drawCreditText();
  disp.drawPolygon(dm.config.numpoints, dm.config.points, xrange, yrange);
  disp.fillMapGaps();

  if(dm.config.point_present) {
    disp.plotDecal(dm.config.decalfile, dm.config.dpoint, xrange, yrange);
  }

  disp.writePng(dm.config.outfile);

  gdImageDestroy(basemap);
  delete[] lats;
  delete[] longs;
  delete[] grid_lats;
  delete[] grid_longs;
  delete[] slmask;
  delete[] drawmask;
  delete[] mask;
  delete[] data;
}

void handleText(Displayer& disp, DataManager& dm) {
  dm.variable = dm.config.xvariable;
  dm.open_datafile();

  FILE* outfd;

  if(dm.config.outfile.length() && dm.config.outfile != "-") {
    outfd = fopen(dm.config.outfile.c_str(), "w");
    if(!outfd) {
      fprintf(stderr, "Couldn't open file %s\n", dm.config.outfile.c_str());
      exit(1);
    }

  } else {
    outfd = stdout;
  }

  int rows = dm.numrows();
  int cols = dm.numcols();
  int datasize = dm.data_size();
  double* data = new double[datasize];
  double* lats = new double[dm.lats_size()];
  double* grid_lats = new double[dm.gridlats_size()];
  double* grid_longs = new double[dm.gridlongs_size()];
  double* longs = new double[dm.longs_size()];
  int* slmask = new int[datasize];
  int* mask = new int[datasize];
  gdImagePtr basemap = dm.get_basemap();
  disp.setOffsets(basemap);

  dm.get_slmask(slmask);
  dm.get_lats(lats);
  dm.get_longs(longs);
  dm.get_gridlats(grid_lats, lats);
  dm.get_gridlongs(grid_longs, longs);
  dm.get_datamask(mask, slmask, grid_lats, grid_longs);
  dm.get_data(data);

  Range datarange(data, datasize);
  Range longrange(grid_longs, dm.gridlongs_size());
  Range latrange(grid_lats, dm.gridlats_size());
  Legend* legend = disp.getLegend(datarange);
  
  fprintf(outfd, "Lat range: (%0.2f)-(%0.2f)\n", latrange.min(), latrange.max());
  fprintf(outfd, "Lon range: (%0.2f)-(%0.2f)\n", longrange.min(), longrange.max());
  fprintf(outfd, "Data range: (%0.6f)-(%0.6f)\n", datarange.min(), datarange.max());
  fprintf(outfd, "Legend colour range: (%f)-(%f)\n", legend->range.min(), legend->range.max());
  fprintf(outfd, "Size: (%i)-(%i)\n", disp.img_width, disp.img_height);
  fprintf(outfd, "Map size: (%i)-(%i)\n", disp.plot_width, disp.plot_height);
  fprintf(outfd, "Map offset: (%i)-(%i)\n", disp.plot_offset_x, disp.plot_offset_y);

  // Data analysis loop stuff

  list<WPoint> wpoints;
  list<WPoint>::const_iterator p_iter;
  double base_data_sum = 0;
  double base_data_avg = 0;
  double total_weight = 0;
  double total_squared_weight = 0;
  double data_avg = 0;
  double data_sum = 0;
  double data_variance = 0;
  double data_stddev = 0;
  double data_median = 0;
  double data_wmedian = 0;
  int i, j;
  
  double* base_data = 0;
  
  LatLonRange llrange;
  
  double* dataptr = data;
  
  if(dm.config.pct_change_data) {
    base_data = new double[datasize];
    dm.get_basedata(base_data);
    
    // Correct the baseline data grid for the % change
    for(i = 0; i < rows; i++) {
      for(j = 0; j < cols; j++) {
	const int offset = (i * cols) + j;
	data[offset] = base_data[offset] + base_data[offset] * (data[offset] / 100);
      }
    }
  }
  
  for(i = 0; i < rows; i++) {
    for(j = 0; j < cols; j++) {
      const int offset = (i * cols) + j;
      if(mask[offset] == 1) {
	const double area = (M_PI / 180) * squared(EARTH_RADIUS) * fabs(sin(grid_lats[i] * (M_PI / 180)) - sin(grid_lats[i + 1] * (M_PI / 180))) * fabs(grid_longs[j] - grid_longs[j + 1]);
	data_sum += area * dataptr[offset];
	total_squared_weight += squared(area);
	total_weight += area;
	llrange.addPoint(data[offset], longs[j], lats[i]);
	wpoints.push_back(WPoint(data[offset], area));
	if(dm.config.pct_change_data) {
	  base_data_sum += area * base_data[offset];
	}
      }
    }
  }
  
  data_avg = data_sum / total_weight;
  
  if(dm.config.pct_change_data) {
    base_data_avg = base_data_sum / total_weight;
  }
  
  // Go through the list of pts and weights getting variance components
  for(p_iter = wpoints.begin(); p_iter != wpoints.end(); p_iter++) {
    data_variance += squared((*p_iter).p - data_avg) * (*p_iter).w;
  }
  
  if(dm.config.pct_change_data) {
    data_avg = ((data_avg - base_data_avg) / base_data_avg) * 100;
  }
  
  // Calculate variance
  // See http://pygsl.sourceforge.net/reference/pygsl/node35.html
  data_variance *= total_weight / (squared(total_weight) - total_squared_weight);
  
  // Calculate standard deviation
  data_stddev = sqrt(data_variance);
  
  // Work out the median
  wpoints.sort();
  double desired_weight = total_weight / 2;
  WPoint curr_pt(0,0);
  WPoint last_pt(0,0);
  double wsum = 0;
  for(p_iter = wpoints.begin(); p_iter != wpoints.end() && wsum < desired_weight; p_iter++) {
    last_pt = curr_pt;
    curr_pt = *p_iter;
    wsum += curr_pt.w;
  }
  
  // Interpolate and find best fit median
  double wdiff = curr_pt.w;
  data_wmedian = ((desired_weight - (wsum - wdiff)) / wdiff) * last_pt.p + ((wsum - desired_weight) / wdiff) * curr_pt.p;
  
  desired_weight = (double)wpoints.size() / 2;
  wsum = 0;
  for(p_iter = wpoints.begin(); p_iter != wpoints.end() && wsum <= desired_weight; p_iter++) {
    last_pt = curr_pt;
    curr_pt = *p_iter;
    wsum++;
  }
  
  // Finish the median calculation
  if(wsum - desired_weight < 1) {
    data_median = curr_pt.p;
  } else {
    // Average of 2
    data_median = (last_pt.p + curr_pt.p) / 2;
  }
  
  fprintf(outfd, "Selection area weighted mean: (%0.6f)\n", data_avg);
  fprintf(outfd, "Selection area weighted median: (%0.6f)\n", data_wmedian);
  fprintf(outfd, "Selection area weighted standard deviation: (%0.6f)\n", data_stddev);
  fprintf(outfd, "Selection median: (%0.6f)\n", data_median);
  fprintf(outfd, "Selection data min (d, lon, lat): (%0.6f)-(%0.2f)-(%0.2f)\n", llrange.getMin(), llrange.getMinLong(), llrange.getMinLat());
  fprintf(outfd, "Selection data max (d, lon, lat): (%0.6f)-(%0.2f)-(%0.2f)\n", llrange.getMax(), llrange.getMaxLong(), llrange.getMaxLat());
  fprintf(outfd, "Selection area: (%.0f) km<sup>2</sup>\n", total_weight);
  fprintf(outfd, "Selection num grid boxes: (%i)\n", wpoints.size());
  
  if(base_data) {
    delete[] base_data;
  }

  // If the user wants to know about the data box that a specific lat/lon 
  // is in, report back
  if(dm.config.point_present) {
    // Find the column it's in
    int i, j;
    for(i = 0; i < cols; i++) {
      if(grid_longs[i] < dm.config.dpoint.x && grid_longs[i + 1] > dm.config.dpoint.x) {
	break;
      }
    }
    
    // Find the row it's in
    for(j = 0; j < rows; j++) {
      if(grid_lats[j + 1] < dm.config.dpoint.y && grid_lats[j] > dm.config.dpoint.y) {
	break;
      }
    }
    
    if(j != rows && i != cols) {
      fprintf(outfd, "Grid box: (%i)-(%i)\n", i, j);
    }
    fprintf(outfd, "Data point (d, lon, lat): (%0.6f)-(%0.2f)-(%0.2f)\n", data[(j * cols) + i], longs[i], lats[j]);
  }

  if(dm.config.outfile.length() && dm.config.outfile != "-") {
    fclose(outfd);
  }

  delete[] lats;
  delete[] longs;
  delete[] grid_lats;
  delete[] grid_longs;
  delete[] slmask;
  delete[] mask;
  delete[] data;
}

void handleMask(Displayer& disp, DataManager& dm) {
  dm.variable = dm.config.xvariable;
  dm.open_datafile();

  FILE* outfd;

  if(dm.config.outfile.length() && dm.config.outfile != "-") {
    outfd = fopen(dm.config.outfile.c_str(), "w");
  } else {
    outfd = stdout;
  }

  int rows = dm.numrows();
  int cols = dm.numcols();
  int datasize = dm.data_size();
  double* lats = new double[dm.lats_size()];
  double* grid_lats = new double[dm.gridlats_size()];
  double* grid_longs = new double[dm.gridlongs_size()];
  double* longs = new double[dm.longs_size()];
  int* slmask = new int[datasize];
  int* mask = new int[datasize];
  gdImagePtr basemap = dm.get_basemap();
  disp.setOffsets(basemap);
  int i, j;

  dm.get_slmask(slmask);
  dm.get_lats(lats);
  dm.get_longs(longs);
  dm.get_gridlats(grid_lats, lats);
  dm.get_gridlongs(grid_longs, longs);
  dm.get_datamask(mask, slmask, grid_lats, grid_longs);

  for(i = 0; i < rows; i++) {
    for(j = 0; j < cols; j++) {
      fprintf(outfd, "% i", mask[(i * cols) + j]);
    }
    fprintf(outfd, "\r\n");
  }

  if(dm.config.outfile.length() && dm.config.outfile != "-") {
    fclose(outfd);
  }

  delete[] lats;
  delete[] longs;
  delete[] grid_lats;
  delete[] grid_longs;
  delete[] slmask;
  delete[] mask;
}

void handleGeoref(Displayer& disp, DataManager& dm) {
  dm.variable = dm.config.xvariable;
  dm.open_datafile();

  FILE* outfd;

  if(dm.config.outfile.length() && dm.config.outfile != "-") {
    outfd = fopen(dm.config.outfile.c_str(), "w");
  } else {
    outfd = stdout;
  }

  int rows = dm.numrows();
  int cols = dm.numcols();
  int datasize = dm.data_size();
  double* lats = new double[dm.lats_size()];
  double* grid_lats = new double[dm.gridlats_size()];
  double* grid_longs = new double[dm.gridlongs_size()];
  double* longs = new double[dm.longs_size()];
  int* slmask = new int[datasize];
  int* mask = new int[datasize];
  gdImagePtr basemap = dm.get_basemap();
  int i, j;
  double* fulldata[NUM_TIMESLICES];

  disp.setOffsets(basemap);

  dm.get_slmask(slmask);
  dm.get_lats(lats);
  dm.get_longs(longs);
  dm.get_gridlats(grid_lats, lats);
  dm.get_gridlongs(grid_longs, longs);
  dm.get_datamask(mask, slmask, grid_lats, grid_longs);

  for(i = 0; i < NUM_TIMESLICES; i++) {
    fulldata[i] = new double[datasize];
    dm.get_data(fulldata[i]);
  }

  for(i = 0; i < rows; i++) {
    for(j = 0; j < cols; j++) {
      const int offset = (i * cols) + j;
      if(mask[offset]) {
	fprintf(outfd, " % #9.4f % #9.4f", lats[i], longs[j]);
	for(int q = 0; q < NUM_TIMESLICES; q++) {
	  fprintf(outfd, " % E", fulldata[q][offset]);
	}
	fprintf(outfd, "\r\n");
      }
    }
  }

  for(i = 0; i < NUM_TIMESLICES; i++) {
    delete[] fulldata[i];
  }

  if(dm.config.outfile.length() && dm.config.outfile != "-") {
    fclose(outfd);
  }

  delete[] lats;
  delete[] longs;
  delete[] grid_lats;
  delete[] grid_longs;
  delete[] slmask;
  delete[] mask;
}

void handlePlotInfo(Displayer& disp, DataManager& dm) {
  dm.variable = dm.config.xvariable;
  dm.open_datafile();
  FILE* outfd;

  if(dm.config.outfile.length() && dm.config.outfile != "-") {
    outfd = fopen(dm.config.outfile.c_str(), "w");
  } else {
    outfd = stdout;
  }

  gdImagePtr basemap = dm.get_basemap();
  double* lats = new double[dm.lats_size()];
  double* longs = new double[dm.longs_size()];
  double* grid_lats = new double[dm.gridlats_size()];
  double* grid_longs = new double[dm.gridlongs_size()];
  disp.setOffsets(basemap);
  
  dm.get_lats(lats);
  dm.get_longs(longs);
  dm.get_gridlats(grid_lats, lats);
  dm.get_gridlongs(grid_longs, longs);

  Range xrange(grid_longs, dm.gridlongs_size());
  Range yrange(grid_lats, dm.gridlats_size());

  fprintf(outfd, "Lat range: (%0.2f)-(%0.2f)\n", yrange.min(), yrange.max());
  fprintf(outfd, "Lon range: (%0.2f)-(%0.2f)\n", xrange.min(), xrange.max());
  fprintf(outfd, "Size: (%i)-(%i)\n", disp.img_width, disp.img_height);
  fprintf(outfd, "Map size: (%i)-(%i)\n", disp.plot_width, disp.plot_height);
  fprintf(outfd, "Map offset: (%i)-(%i)\n", disp.plot_offset_x, disp.plot_offset_y);

  if(dm.config.outfile.length() && dm.config.outfile != "-") {
    fclose(outfd);
  }

  gdImageDestroy(basemap);
  delete[] lats;
  delete[] longs;
  delete[] grid_lats;
  delete[] grid_longs;
}

void handleRegionOnly(Displayer& disp, DataManager& dm) {
  dm.variable = dm.config.xvariable;
  dm.open_datafile();
  FILE* outfd;

  if(dm.config.outfile.length() && dm.config.outfile != "-") {
    outfd = fopen(dm.config.outfile.c_str(), "w");
  } else {
    outfd = stdout;
  }

  gdImagePtr basemap = dm.get_basemap();

  double* lats = new double[dm.lats_size()];
  double* longs = new double[dm.longs_size()];
  double* grid_lats = new double[dm.gridlats_size()];
  double* grid_longs = new double[dm.gridlongs_size()];

  dm.get_lats(lats);
  dm.get_longs(longs);
  dm.get_gridlats(grid_lats, lats);
  dm.get_gridlongs(grid_longs, longs);

  Range xrange(grid_longs, dm.gridlongs_size());
  Range yrange(grid_lats, dm.gridlats_size());

  disp.setOffsets(basemap);
  disp.createCanvas(dm.config.fontfile);
  disp.copyMap(basemap);
  disp.drawTicks(xrange, yrange);
  disp.fillMapGaps();
  disp.fillGapsRegionMap();
  disp.clearIdentifyArea();
  disp.drawCreditText();
  disp.drawPolygon(dm.config.numpoints, dm.config.points, xrange, yrange);

  if(dm.config.point_present) {
    disp.plotDecal(dm.config.decalfile, dm.config.dpoint, xrange, yrange);
  }

  disp.writePng(dm.config.outfile);

  if(dm.config.outfile.length() && dm.config.outfile != "-") {
    fclose(outfd);
  }

  gdImageDestroy(basemap);
  delete[] lats;
  delete[] longs;
  delete[] grid_lats;
  delete[] grid_longs;
}

inline bool is_available(string token) {
  return !(token == "o" || token == "0");
}

void readGCMInfo(DataManager& dm, list<string*>& list) {
  ifstream f(dm.config.gcminfofile.c_str());
  if(!f.is_open()) {
    cerr << "Could not read gcminfo file!" << endl;
    exit(3);
  }

  std::string filedata((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());

  boost::tokenizer<boost::escaped_list_separator<char> > tok(filedata);
  boost::tokenizer<boost::escaped_list_separator<char> >::iterator beg=tok.begin();

  int (*lowercase)(int) = tolower;
  while(beg != tok.end()) {
    int count;
    string* bits = new string[GCMINFO_COLS];
    // Read a line of vars
    string model = *beg;
    std::transform(model.begin(), model.end(), model.begin(), lowercase);
    model.erase(0, model.find_first_not_of("\n\r"));
    bits[0] = model;
    ++beg;
    for(count = 1; count < GCMINFO_COLS && beg != tok.end(); ++beg, ++count)
      bits[count] = *beg;
    if(count == GCMINFO_COLS) {
      list.push_back(bits);
    } else {
      delete[] bits;
    }
  }
}

#define DKRED 0x00800000
#define LTRED 0x00FF0000
#define YELLOW 0x00FFFF00
#define DKGREEN 0x00008000
#define LTGREEN 0x0000FF00
#define CYAN 0x0000FFFF
#define DKBLUE 0x00000080
#define LTBLUE 0x000000FF
#define MAGENTA 0x00FF00FF

#define WHITE 0x00FFFFFF
#define DKGRAY 0x00808080
#define LTGRAY 0x00C0C0C0
#define BLACK 0x00000000

void loadData(DataManager& dm, list<ScatterVars* >& vars, list<LegendToken* >& legend_bits, Range& xrange, Range& yrange) {
  int datasize = 0;
  double* xdata = 0;
  double* ydata = 0;
  double* lats = 0;
  double* longs = 0;
  double* grid_lats = 0;
  double* grid_longs = 0;
  int* slmask = 0;
  int* mask = 0;
  int x = 0, y = 0;

  map<const string, int> colours;
  map<const string, enum SYMBOL> symbols;

  colours["ccsrnies"] = LTGREEN;
  colours["cgcm1"] = BLACK;
  colours["cgcm2"] = DKGRAY;
  colours["csiromk2b"] = MAGENTA;
  colours["echam4"] = LTRED;
  colours["hadcm2"] = LTBLUE;
  colours["hadcm3"] = DKBLUE;
  colours["gfdlr15"] = CYAN;
  colours["gfdlr30"] = CYAN;
  colours["ncarpcm"] = YELLOW;
  colours["hadcm3"] = DKBLUE;
  colours["hadcm3"] = DKBLUE;

  symbols["A1"] = LTRIANGLE;
  symbols["A2"] = UTRIANGLE;
  symbols["B1"] = RTRIANGLE;
  symbols["B2"] = DIAMOND;
  symbols["A1FI"] = STAR6;
  symbols["A1T"] = DTRIANGLE;
  symbols["ga"] = CIRCLE;
  symbols["gg"] = SQUARE;

  string last_model = "";
  string last_expt = "";
  string last_nmodel = "";
  string last_nexpt = "";
  LegendToken* tok = 0;
  LegendToken* special = 0;
  LegendToken* ensemble = 0;
  list<ScatterVars*>::iterator vars_iter;

  for(vars_iter = vars.begin(); vars_iter != vars.end(); vars_iter++) {
    // Set symbols to use on map and legend
    if((*vars_iter)->model != last_model || (*vars_iter)->expt != last_expt) {
      int colour;
      enum SYMBOL symbol;

      if(dm.config.percentiles) {
	symbol = NONE;
	colour = LTGRAY;
      } else {
	string id = (*vars_iter)->expt.substr(0, 2);
	if(symbols.find(id) != symbols.end()) {
	  symbol = symbols[id];
	} else {
	  symbol = NONE;
	}
	colour = colours[(*vars_iter)->model];
      }
      if(symbols.find((*vars_iter)->expt) != symbols.end()) {
	// Special case: If this is A1FI or such...
	if(dm.config.percentiles) {
	  symbol = NONE;
	} else {
	  symbol = symbols[(*vars_iter)->expt];
	}
	special = new LegendToken((*vars_iter)->model + " " + (*vars_iter)->expt, colour, symbol, true);
	legend_bits.push_back(special);
      } else if((*vars_iter)->expt.substr(2, 1) == "x") {
	// Ensemble run
	ensemble = new LegendToken((*vars_iter)->model + " " + (*vars_iter)->expt, colour, symbol, false);
	legend_bits.push_back(ensemble);
      } else {
	// Normal...

	// If the first 2 digits of the experiment match, tack the last digit onto the end of the name
	if((*vars_iter)->model == last_nmodel && (*vars_iter)->expt.substr(0, 2) == last_nexpt.substr(0,2)) {
	  tok->name += "," + (*vars_iter)->expt.substr(2, 1);
	} else {
	  // If they aren't the same, it's time to push a new token up
	  tok = new LegendToken((*vars_iter)->model + " " + (*vars_iter)->expt, colour, symbol, true);
	  legend_bits.push_back(tok);
	}

	last_nmodel = (*vars_iter)->model;
	last_nexpt = (*vars_iter)->expt;
      }
    }

    if(symbols.find((*vars_iter)->expt) != symbols.end()) {
      (*vars_iter)->symbol = special;
    } else if((*vars_iter)->expt.substr(2, 1) == "x") {
      (*vars_iter)->symbol = ensemble;
    } else {
      (*vars_iter)->symbol = tok;
    }

    // Set the appropriate data file bits
    dm.loadScatterVars(*vars_iter);

    // If model change...
    if((*vars_iter)->model != last_model) {
      dm.open_datafile();
      datasize = dm.data_size();
      
      if(xdata)
	delete[] xdata;
      if(ydata)
	delete[] ydata;
      if(lats)
	delete[] lats;
      if(longs)
	delete[] longs;
      
      xdata = new double[datasize];
      ydata = new double[datasize];
      lats = new double[dm.lats_size()];
      longs = new double[dm.longs_size()];

      dm.get_lats(lats);
      dm.get_longs(longs);

      if(dm.config.scatter_type == ST_POINT) {
	x = nearest_offset(dm.longs_size(), longs, dm.config.dpoint.x);
	y = nearest_offset(dm.lats_size(), lats, dm.config.dpoint.y);
      } else if(dm.config.scatter_type == ST_REGION) {
	if(slmask)
	  delete[] slmask;
	if(mask)
	  delete[] mask;
	if(grid_lats)
	  delete[] grid_lats;
	if(grid_longs)
	  delete[] grid_longs;

	slmask = new int[datasize];
	mask = new int[datasize];
	grid_lats = new double[dm.gridlats_size()];
	grid_longs = new double[dm.gridlongs_size()];

	dm.get_slmask(slmask);
	dm.get_gridlats(grid_lats, lats);
	dm.get_gridlongs(grid_longs, longs);
	dm.get_datamask(mask, slmask, grid_lats, grid_longs);
      }
    }

    // Load the data in
    dm.get_data(ydata);
    if((*vars_iter)->xvariable != "") {
      dm.loadScatterVars(*vars_iter, true);
      dm.get_data(xdata);
    }

    // Grab the data we want
    if(dm.config.scatter_type == ST_POINT) {
      // Load data at grid point
      (*vars_iter)->setCoord(longs[x], lats[y]);
      if((*vars_iter)->xvariable != "") {
	(*vars_iter)->setXData(xdata[(y * dm.numcols()) + x]);
	xrange.add(xdata[(y * dm.numcols()) + x]);
      }
      (*vars_iter)->setYData(ydata[(y * dm.numcols()) + x]);
      yrange.add(ydata[(y * dm.numcols()) + x]);
      //cout << longs[x] << ", " << lats[y] << ": " << ydata[(y * dm.numcols()) + x] << endl;
    } else if(dm.config.scatter_type == ST_REGION) {
      // Load in required data
      const int rows = dm.lats_size();
      const int cols = dm.longs_size();
      double xdata_avg;
      double ydata_avg;
      double total_weight = 0;
      double xdata_sum = 0;
      double ydata_sum = 0;

      // Calculate weighted average over area; use that as data point
      if((*vars_iter)->xvariable == "") {
	// If we don't have an X variable...
	for(int i = 0; i < rows; i++) {
	  for(int j = 0; j < cols; j++) {
	    const int offset = (i * cols) + j;
	    if(mask[offset] == 1) {
	      const double area = (M_PI / 180) * squared(EARTH_RADIUS) * fabs(sin(grid_lats[i] * (M_PI / 180)) - sin(grid_lats[i + 1] * (M_PI / 180))) * fabs(grid_longs[j] - grid_longs[j + 1]);
	      ydata_sum += area * ydata[offset];
	      total_weight += area;
	    }
	  }
	}
	ydata_avg = ydata_sum / total_weight;
	yrange.add(ydata_avg);
	(*vars_iter)->setYData(ydata_avg);
      } else {
	// If we do have an X variable...
	for(int i = 0; i < rows; i++) {
	  for(int j = 0; j < cols; j++) {
	    const int offset = (i * cols) + j;
	    if(mask[offset] == 1) {
	      const double area = (M_PI / 180) * squared(EARTH_RADIUS) * fabs(sin(grid_lats[i] * (M_PI / 180)) - sin(grid_lats[i + 1] * (M_PI / 180))) * fabs(grid_longs[j] - grid_longs[j + 1]);
	      xdata_sum += area * xdata[offset];
	      ydata_sum += area * ydata[offset];
	      total_weight += area;
	    }
	  }
	}
	xdata_avg = xdata_sum / total_weight;
	ydata_avg = ydata_sum / total_weight;
	xrange.add(xdata_avg);
	yrange.add(ydata_avg);
	(*vars_iter)->setXData(xdata_avg);
	(*vars_iter)->setYData(ydata_avg);
      }
      (*vars_iter)->setCoord(0, 0);
    }

    last_expt = (*vars_iter)->expt;
    last_model = (*vars_iter)->model;
  }
  if(xdata)
    delete[] xdata;
  if(ydata)
    delete[] ydata;
  if(lats)
    delete[] lats;
  if(longs)
    delete[] longs;
  if(slmask)
    delete[] slmask;
  if(mask)
    delete[] mask;
  if(grid_lats)
    delete[] grid_lats;
  if(grid_longs)
    delete[] grid_longs;
}

void setRanges(Range& xrange, Range& yrange) {
  double xtick_spacing = tick_spacing(xrange, DESIRED_XTICKS);
  xrange.setmin(floor(xrange.min() / xtick_spacing) * xtick_spacing);
  xrange.setmax(ceil(xrange.max() / xtick_spacing) * xtick_spacing);
  if(xrange.min() > 0) xrange.setmin(0);
  double ytick_spacing = tick_spacing(yrange, DESIRED_YTICKS);
  yrange.setmin(floor(yrange.min() / ytick_spacing) * ytick_spacing);
  yrange.setmax(ceil(yrange.max() / ytick_spacing) * ytick_spacing);
  if(yrange.min() > 0) yrange.setmin(0);
}

bool compareScatterYVar(const ScatterVars* a, const ScatterVars* b) {
  return a->daty < b->daty;
}

void pctile(list<ScatterVars*>& vars, double pct, ScatterVars* v) {
  int element = (int)floor(pct * (vars.size() - 1));
  list<ScatterVars*>::const_iterator i;
  int j;
  for(j = 0, i = vars.begin(); j < element; ++i, ++j);
  v->datx = (*i)->datx;
  v->daty = (*i)->daty;
  v->lon = (*i)->lon;
  v->lat = (*i)->lat;
}

void calcPercentiles(list<ScatterVars*>& vars, list<LegendToken* >& leg_tokens) {
  const double percentiles[] = { 0.1, 0.5, 0.9 };
  const string years[] = { "2020", "2050", "2080" };
  const string desc[] = { "10th percentile", "Median", "90th percentile" };
  const enum SYMBOL symbols[] = { DTRIANGLE, DIAMOND, UTRIANGLE };

  // Set up the list-map combination
  map<const string, list<ScatterVars*> > vars_by_ts;
  list<ScatterVars*> l[3];
  for(int i = 0; i < 3; i++) {
    vars_by_ts[years[i]] = l[i];
  }

  // Populate the lists in the map
  list<ScatterVars*>::const_iterator i;
  for(i = vars.begin(); i != vars.end(); ++i) {
    vars_by_ts[(*i)->timeslice].push_back(*i);
  }

  // Sort the lists
  for(int i = 0; i < 3; i++) {
    vars_by_ts[years[i]].sort(compareScatterYVar);
  }

  // Create new symbols and give them data
  for(int i = 0; i < 3; i++) {
    LegendToken* tok = new LegendToken(desc[i], 0x00000000, symbols[i], true);
    leg_tokens.push_back(tok);
    for(int j = 0; j < 3; j++) {
      ScatterVars* v = new ScatterVars(desc[i], "", years[j], "", "");
      pctile(vars_by_ts[years[j]], percentiles[i], v);
      v->symbol = tok;
      vars.push_back(v);
    }
  }
}

enum GCMINFO_OFFSETS{MODEL_OFFSET, MODEL_COUNT, JUNK1, SERIES_OFFSET, MODELNAME_OFFSET, EXPT_OFFSET, S2020_OFFSET, S2050_OFFSET, S2080_OFFSET};
void handleScatterTimeslice(Displayer& disp, DataManager& dm, bool textOnly = false) {
  int count = 0;
  int var_offset = -1;

  if(dm.config.scatter_type == ST_POINT && !dm.config.point_present) {
    cerr << "Must specify point when doing a scatter plot using a data point" << endl;
    exit(1);
  }

  // Read in GCMINFO file
  list<string*> gcminfo;
  readGCMInfo(dm, gcminfo);
  list<string*>::iterator beg = gcminfo.begin();

  // Read header and identify variable
  for(count = 0; count < GCMINFO_COLS; ++count) {
    if((*beg)[count] == dm.config.yvariable)
      var_offset = count;
  }
  delete[] (*beg);
  ++beg;

  // Create a list of the available expts
  list<ScatterVars* > vars;
  for(; beg != gcminfo.end(); ++beg) {
    string model = (*beg)[MODEL_OFFSET];
    string expt = (*beg)[EXPT_OFFSET];
    bool var_available = is_available((*beg)[var_offset]);

    if(var_available && (dm.config.scenario_set == "" || dm.config.scenario_set == (*beg)[SERIES_OFFSET])) {
      if(is_available((*beg)[S2020_OFFSET])) {
	ScatterVars* v = new ScatterVars(model, expt, "2020", "", dm.config.yvariable);
	v->setXData(2020);
	vars.push_back(v);
      }
      if(is_available((*beg)[S2050_OFFSET])) {
	ScatterVars* v = new ScatterVars(model, expt, "2050", "", dm.config.yvariable);
	v->setXData(2050);
	vars.push_back(v);
      }
      if(is_available((*beg)[S2080_OFFSET])) {
	ScatterVars* v = new ScatterVars(model, expt, "2080", "", dm.config.yvariable);
	v->setXData(2080);
	vars.push_back(v);
      }
    }
    delete[] (*beg);
  }

  // Load in the data
  Range xrange;
  Range yrange;
  list<LegendToken* > leg_tokens;

  loadData(dm, vars, leg_tokens, xrange, yrange);

  if(disp.range_dynamic) {
    setRanges(xrange, yrange);
  } else {
    yrange.setmin(disp.yrange_min);
    yrange.setmax(disp.yrange_max);
  }
  xrange.setmin(2010);
  xrange.setmax(2090);
  
  if(dm.config.percentiles) {
    calcPercentiles(vars, leg_tokens);
  }
  
  if(textOnly) {
    FILE* out = fopen(dm.config.outfile.c_str(), "w");

    if(out) {
      list<ScatterVars*>::iterator i = vars.begin();
      for(; i != vars.end(); i++) {
	ScatterVars* s = *i;
	fprintf(out, "%s %s,%0.2f,%0.2f,%s,%0.6f,\n", s->model.c_str(), s->expt.c_str(), s->lon, s->lat, s->timeslice.c_str(), s->daty);
      }      
      fclose(out);
    }
  } else {
    // Set up the plot
    disp.setScatterOffsets();
    disp.setTicks(xrange, yrange);
    disp.createCanvas(dm.config.fontfile);
    
    // Clear the plot area, draw what we want to...
    disp.clearPlot();
    disp.drawScatterGrid(xrange, yrange);
    disp.drawLines(vars, xrange, yrange);
    disp.drawScatter(vars, xrange, yrange);
    
    // Draw tick marks, text labels, and title up X and Y axes
    disp.drawTicks(xrange, yrange);
    disp.drawAxisTitles();
    
    // Draw identifying text
    disp.clearIdentifyArea();
    disp.drawCreditText();
    disp.drawIdentifyText();
    
    // Draw the legend
    disp.drawLegend(leg_tokens);
    
    // Fill in gaps in the map
    disp.fillMapGaps();
    
    // Write out the file
    disp.writePng(dm.config.outfile);
  }

  // Clean up
  list<ScatterVars* >::iterator vars_iter;
  for(vars_iter = vars.begin(); vars_iter != vars.end(); ++vars_iter) {
    delete *vars_iter;
  }

  list<LegendToken* >::iterator lt_iter;
  for(lt_iter = leg_tokens.begin(); lt_iter != leg_tokens.end(); ++lt_iter) {
    delete *lt_iter;
  }
}

void handleScatterVariable(Displayer& disp, DataManager& dm, bool textOnly = false) {
  int count = 0;
  int xvar_offset = -1;
  int yvar_offset = -1;

  map<const string, int> ts;

  ts["2020"] = S2020_OFFSET;
  ts["2050"] = S2050_OFFSET;
  ts["2080"] = S2080_OFFSET;

  if(dm.config.scatter_type == ST_POINT && !dm.config.point_present) {
    cerr << "Must specify point when doing a scatter plot using a data point" << endl;
    exit(1);
  }
  
  if(dm.timeslice == "" || dm.region == "" || dm.config.xvariable == "" || dm.config.yvariable == "") {
    cerr << "Must specify timeslice, region, xvariable, and yvariable" << endl;
    exit(1);
  }

  // Read in GCMINFO file
  list<string*> gcminfo;
  readGCMInfo(dm, gcminfo);
  list<string*>::iterator beg = gcminfo.begin();

  // Read header and identify variables
  for(count = 0; count < GCMINFO_COLS; ++count) {
    if((*beg)[count] == dm.config.xvariable)
      xvar_offset = count;
    if((*beg)[count] == dm.config.yvariable)
      yvar_offset = count;
  }
  delete[] (*beg);
  ++beg;

  // Create a list of the available expts
  list<ScatterVars* > vars;
  for(; beg != gcminfo.end(); ++beg) {
    string model = (*beg)[MODEL_OFFSET];
    string expt = (*beg)[EXPT_OFFSET];
    bool xvar_available = is_available((*beg)[xvar_offset]);
    bool yvar_available = is_available((*beg)[yvar_offset]);
    bool ts_available = dm.timeslice == "1961_1990" || (ts.end() != ts.find(dm.timeslice) && is_available((*beg)[ts[dm.timeslice]]));

    if(xvar_available && yvar_available && ts_available && (dm.config.scenario_set == "" || dm.config.scenario_set == (*beg)[SERIES_OFFSET])) {
      ScatterVars* v = new ScatterVars(model, expt, dm.timeslice, dm.config.xvariable, dm.config.yvariable);
      vars.push_back(v);
    }
    delete[] (*beg);
  }

  // Load in the data
  Range xrange;
  Range yrange;
  list<LegendToken* > leg_tokens;

  loadData(dm, vars, leg_tokens, xrange, yrange);

  if(disp.range_dynamic) {
    setRanges(xrange, yrange);
  } else {
    xrange.setmin(disp.xrange_min);
    xrange.setmax(disp.xrange_max);
    yrange.setmin(disp.yrange_min);
    yrange.setmax(disp.yrange_max);
  }

  if(textOnly) {
    FILE* out = fopen(dm.config.outfile.c_str(), "w");

    if(out) {
      list<ScatterVars*>::iterator i = vars.begin();
      for(; i != vars.end(); i++) {
	ScatterVars* s = *i;
	fprintf(out, "%s %s,%0.2f,%0.2f,%0.6f,%0.6f,\n", s->model.c_str(), s->expt.c_str(), s->lon, s->lat, s->datx, s->daty);
      }      
      fclose(out);
    }
  } else {
    // Set up the plot
    disp.setScatterOffsets();
    disp.setTicks(xrange, yrange);
    disp.createCanvas(dm.config.fontfile);
    
    // Clear the plot area, draw what we want to...
    disp.clearPlot();
    disp.drawScatterGrid(xrange, yrange);
    disp.drawScatter(vars, xrange, yrange);
    
    // Draw tick marks, text labels, and title up X and Y axes
    disp.drawTicks(xrange, yrange);
    disp.drawAxisTitles();
    
    // Draw identifying text
    disp.clearIdentifyArea();
    disp.drawCreditText();
    disp.drawIdentifyText();
    
    // Draw the legend
    disp.drawLegend(leg_tokens);
    
    // Fill in gaps in the map
    disp.fillMapGaps();
    
    // Write out the file
    disp.writePng(dm.config.outfile);
  }

  // Clean up
  list<ScatterVars* >::iterator vars_iter;
  for(vars_iter = vars.begin(); vars_iter != vars.end(); ++vars_iter) {
    delete *vars_iter;
  }

  list<LegendToken* >::iterator lt_iter;
  for(lt_iter = leg_tokens.begin(); lt_iter != leg_tokens.end(); ++lt_iter) {
    delete *lt_iter;
  }
}

int main(int argc, char ** argv) {
  Config c;
  Displayer disp;
  DataManager dm(c);

  ConfigFile cf("/usr/local/etc/genimage.cfg");

  cf.readInto(c.gcminfofile, "gcminfofile");
  cf.readInto(c.fontfile, "fontfile");
  cf.readInto(c.decalfile, "decalfile");
  cf.readInto(disp.credit_text, "credit_text");
  cf.readInto(c.data_dir, "data_dir");
  cf.readInto(c.map_dir, "map_dir");

  parseArgs(c, disp, dm, argc, argv);
    
  switch(c.plot_type) {
  case TYPE_ALL:
    handleMap(disp, dm);
    break;
    
  case TYPE_TEXT:
    handleText(disp, dm);
    break;
    
  case TYPE_MASK:
    handleMask(disp, dm);
    break;
    
  case TYPE_GEOREF:
    handleGeoref(disp, dm);
    break;
    
  case TYPE_PLOTINFO:
    handlePlotInfo(disp, dm);
    break;
    
  case TYPE_REGIONONLY:
    handleRegionOnly(disp, dm);
    break;

  case TYPE_SCATTER_TIMESLICE:
    handleScatterTimeslice(disp, dm);
    break;
    
  case TYPE_SCATTER_VARIABLE:
    handleScatterVariable(disp, dm);
    break;
    
  case TYPE_SCATTER_TIMESLICE_TEXT:
    handleScatterTimeslice(disp, dm, true);
    break;
    
  case TYPE_SCATTER_VARIABLE_TEXT:
    handleScatterVariable(disp, dm, true);
    break;
    
  default:
    fprintf(stderr, "Invalid plot type!\n");
    break;
  }

  if(c.points) {
    for(int i = 0; i < c.numpoints; i++)
      delete c.points[i];
    delete[] c.points;
  }
}
