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

using namespace std;

void parseArgs(Config& c, Displayer& disp, DataManager& dm, int argc, char** argv) {
  char val;
  Point* temp;
  list<Point* > lpoints;
  char* optstring = "a:b:c:e:f:g:h:i:j:k:l:m:n:o:p:r:s:t:u:v:w:x:y:z:";
  static struct option options[] = {
    { "range-min", 1, 0, 'a'},
    { "variable", 1, 0, 'b'},
    { "colour-map", 1, 0, 'c'},
    { "range-max", 1, 0, 'e'},
    { "identify-text", 1, 0, 'f'},
    { "ocean-plot", 1, 0, 'g'},
    { "scatter-type", 1, 0, 'h'},
    { "legend-decimal-places", 1, 0, 'i'},
    { "poly-point", 1, 0, 'j'},
    { "y-axis-text", 1, 0, 'k'},
    { "x-axis-text", 1, 0, 'l'},
    { "model", 1, 0, 'm'},
    { "region", 1, 0, 'n'},
    { "output-file", 1, 0, 'o'},
    { "plot-type", 1, 0, 'p'},
    { "resolution", 1, 0, 'r'},
    { "expt", 1, 0, 's'},
    { "timeslice", 1, 0, 't'},
    //{ "lon-offset", 1, 0, 'u'},
    //{ "lat-offset", 1, 0, 'v'},
    { "data-point", 1, 0, 'x'},
    { "timeofyear", 1, 0, 'y'},
    { "zoom-factor", 1, 0, 'z'},
    { "show-grid", 0, &disp.grid, 1},
    { "colour-map-inverse", 0, &disp.colour_map_rev, 1},
    { "dynamic-range", 0, &disp.range_dynamic, 1},
    { "percent-change-calculations", 0, &c.pct_change_data, 1},
    { 0, 0, 0, 0 }
  };
  while((val = getopt_long(argc, argv, optstring, options, NULL)) != -1) {
    switch(val) {
    case 0:
      break;
    case 'a':
      // range-min
      disp.range_min = atof(optarg);
      break;
    case 'b':
      // variable
      dm.variable = optarg;
      break;
    case 'c':
      // colour-map
      disp.colour_map = atoi(optarg);
      break;
    case 'e':
      // range-max
      disp.range_max = atof(optarg);
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
    case 'i':
      // legend-decimal-places
      disp.leg_dec_places = atoi(optarg);
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
      // lat-offset
      c.offset_lat = atof(optarg);
      break;
    case 'v':
      // lon-offset
      c.offset_lon = atof(optarg);
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
  if(dm.config.numpoints > 2) {
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
    }
  }

  cout << "Done reading GCM info" << endl;
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

void loadData(DataManager& dm, list<ScatterVars* >& vars, list<LegendToken* >& legend_bits) {
  int datasize = 0;
  double* data = 0;
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

  cout << "Started loading data" << endl;
  
  string last_model = "";
  string last_expt = "";
  string last_nmodel = "";
  string last_nexpt = "";
  LegendToken* tok = 0;
  LegendToken* special = 0;
  LegendToken* ensemble = 0;
  list<ScatterVars*>::iterator vars_iter;
  for(vars_iter = vars.begin(); vars_iter != vars.end(); vars_iter++) {
    //cout << "Data: '" << (*vars_iter)->model << "' '" << (*vars_iter)->expt << "' '" << (*vars_iter)->timeslice << "' '" << (*vars_iter)->variable << "'\n";
    
    // Set symbols to use on map and legend
    if((*vars_iter)->model != last_model || (*vars_iter)->expt != last_expt) {
      if(symbols.find((*vars_iter)->expt) != symbols.end()) {
	// Special case: If this is A1FI or such...
	special = new LegendToken((*vars_iter)->model + " " + (*vars_iter)->expt, colours[(*vars_iter)->model], symbols[(*vars_iter)->expt], true);
	legend_bits.push_back(special);
	
      } else if((*vars_iter)->expt.substr(2, 1) == "x") {
	// Ensemble run
	ensemble = new LegendToken((*vars_iter)->model + " " + (*vars_iter)->expt, colours[(*vars_iter)->model], symbols[(*vars_iter)->expt.substr(0, 2)], false);
	legend_bits.push_back(ensemble);
      } else {
	// Normal...
	if((*vars_iter)->model == last_nmodel && (*vars_iter)->expt.substr(0, 2) == last_nexpt.substr(0,2)) {
	  tok->name += "," + (*vars_iter)->expt.substr(2, 1);
	} else {
	  // Reset tok
	  tok = new LegendToken((*vars_iter)->model + " " + (*vars_iter)->expt, colours[(*vars_iter)->model], symbols[(*vars_iter)->expt.substr(0, 2)], true);
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
      
      if(data)
	delete[] data;
      if(lats)
	delete[] lats;
      if(longs)
	delete[] longs;
      
      data = new double[datasize];
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
    dm.get_data(data);

    // Grab the data we want
    if(dm.config.scatter_type == ST_POINT) {
      // Load data at grid point
      //cout << longs[x] << ", " << lats[y] << ": " << data[(y * dm.numcols()) + x] << endl;
      (*vars_iter)->setCoord(longs[x], lats[y]);
      (*vars_iter)->setYData(data[(y * dm.numcols()) + x]);
    } else if(dm.config.scatter_type == ST_REGION) {
      // Load in required data
      const int rows = dm.lats_size();
      const int cols = dm.longs_size();
      double data_avg;
      double total_weight = 0;
      double data_sum = 0;

      // Calculate weighted average over area; use that as data point
      for(int i = 0; i < rows; i++) {
	for(int j = 0; j < cols; j++) {
	  const int offset = (i * cols) + j;
	  if(mask[offset] == 1) {
	    const double area = (M_PI / 180) * squared(EARTH_RADIUS) * fabs(sin(grid_lats[i] * (M_PI / 180)) - sin(grid_lats[i + 1] * (M_PI / 180))) * fabs(grid_longs[j] - grid_longs[j + 1]);
	    data_sum += area * data[offset];
	    total_weight += area;
	  }
	}
      }
      data_avg = data_sum / total_weight;
      //      cout << data_avg << endl;
      (*vars_iter)->setCoord(0, 0);
      (*vars_iter)->setYData(data_avg);
    }

    last_expt = (*vars_iter)->expt;
    last_model = (*vars_iter)->model;
  }
  if(data)
    delete[] data;
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

  cout << "Done loading data" << endl;
}

enum GCMINFO_OFFSETS{MODEL_OFFSET, MODEL_COUNT, JUNK1, SERIES_OFFSET, MODELNAME_OFFSET, EXPT_OFFSET, S2020_OFFSET, S2050_OFFSET, S2080_OFFSET};
void handleScatterTimeslice(Displayer& disp, DataManager& dm) {
  int count = 0;
  int var_offset = -1;

  if(dm.config.scatter_type == ST_POINT && !dm.config.point_present) {
    cerr << "Must specify point when doing a scatter plot using a data point" << endl;
  }

  // Read in GCMINFO file
  list<string*> gcminfo;
  readGCMInfo(dm, gcminfo);
  list<string*>::iterator beg = gcminfo.begin();

  // Read header and identify variable
  for(count = 0; count < GCMINFO_COLS; ++count) {
    if((*beg)[count] == dm.variable)
      var_offset = count;
    cout << (*beg)[count] << " " << count << endl;
  }
  ++beg;

  // Create a list of the available expts
  list<ScatterVars* > vars;
  for(; beg != gcminfo.end(); ++beg) {
    string model = (*beg)[MODEL_OFFSET];
    string expt = (*beg)[EXPT_OFFSET];
    bool var_available = is_available((*beg)[var_offset]);

    if(var_available) {
      if(is_available((*beg)[S2020_OFFSET])) {
	ScatterVars* v = new ScatterVars(model, expt, "2020", dm.variable);
	v->setXData(2020);
	vars.push_back(v);
      }
      if(is_available((*beg)[S2050_OFFSET])) {
	ScatterVars* v = new ScatterVars(model, expt, "2050", dm.variable);
	v->setXData(2050);
	vars.push_back(v);
      }
      if(is_available((*beg)[S2080_OFFSET])) {
	ScatterVars* v = new ScatterVars(model, expt, "2080", dm.variable);
	v->setXData(2080);
	vars.push_back(v);
      }
    }
  }

  // Load in the data
  Range xrange(2010, 2090);
  Range yrange(disp.range_min, disp.range_max);
  list<LegendToken* > leg_tokens;

  loadData(dm, vars, leg_tokens);

  // Set up the plot
  disp.setScatterOffsets();
  disp.setTicks(xrange, yrange);
  disp.createCanvas(dm.config.fontfile);

  // Clear the plot area, draw what we want to...
  disp.clearPlot();
  disp.drawScatterGrid(xrange, yrange);
  disp.drawScatter(vars, xrange, yrange);
  disp.drawLines(vars, xrange, yrange);

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

int main(int argc, char ** argv) {
  Config c;
  Displayer disp;
  DataManager dm(c);

  ConfigFile cf("genimage.cfg");

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
    
  default:
    fprintf(stderr, "Invalid plot type!\n");
    break;
  }
}
