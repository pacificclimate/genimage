#include "genimage.h"
#include "canvas.h"
#include "datamanager.h"
#include "displayer.h"
#include "config.h"
#include "ConfigFile.h"

/* Equatorial Cylindrical Equidistant Projection
    * lat = ((-y + center_y) / height) * 180
    * lon = ((x - center_x) / width) * 360
    * x = ((lon / 360) * width) + center_x
    * y = ((-lat / 180) * height) + center_y
*/

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
    //{ "map-height", 1, 0, 'h'},
    { "legend-decimal-places", 1, 0, 'i'},
    { "poly-point", 1, 0, 'j'},
    //{ "lon-text-spacing", 1, 0, 'k'},
    { "legend-text", 1, 0, 'l'},
    { "model", 1, 0, 'm'},
    { "region", 1, 0, 'n'},
    { "output-file", 1, 0, 'o'},
    { "plot-type", 1, 0, 'p'},
    { "resolution", 1, 0, 'r'},
    { "scenario", 1, 0, 's'},
    { "timeslice", 1, 0, 't'},
    //{ "lon-offset", 1, 0, 'u'},
    //{ "lat-offset", 1, 0, 'v'},
    //{ "map-width", 1, 0, 'w'},
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
      // map-height
      disp.map_height = disp.img_height = atoi(optarg);
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
      // lon-text-spacing
      disp.lon_text_spacing = atoi(optarg);
      break;
    case 'l':
      // legend-text
      disp.leg_text = optarg;
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
      // lat-text-spacing
      c.resolution = atoi(optarg);
      break;
    case 's':
      // scenario
      dm.scenario = optarg;
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
      // map-width
      disp.map_width = disp.img_width = atoi(optarg);
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
  disp.drawIdentifyText();
  disp.drawPolygon(dm.config.numpoints, dm.config.points, xrange, yrange);
  disp.fillGaps();

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
  fprintf(outfd, "Map size: (%i)-(%i)\n", disp.map_width, disp.map_height);
  fprintf(outfd, "Map offset: (%i)-(%i)\n", disp.map_offset_x, disp.map_offset_y);

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
  fprintf(outfd, "Map size: (%i)-(%i)\n", disp.map_width, disp.map_height);
  fprintf(outfd, "Map offset: (%i)-(%i)\n", disp.map_offset_x, disp.map_offset_y);

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
  disp.fillGaps();
  disp.fillGapsRegionMap();
  disp.drawIdentifyText();
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

int main(int argc, char ** argv) {
  Config c;
  Displayer disp;
  DataManager dm(c);

  ConfigFile cf("genimage.cfg");

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
    
  default:
    fprintf(stderr, "Invalid plot type!\n");
    break;
  }
}
