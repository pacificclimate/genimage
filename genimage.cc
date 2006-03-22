#include "genimage.h"

/* Equatorial Cylindrical Equidistant Projection
    * lat = ((-y + center_y) / height) * 180
    * lon = ((x - center_x) / width) * 360
    * x = ((lon / 360) * width) + center_x
    * y = ((-lat / 180) * height) + center_y
*/

/* Need dynamic and user-specified legend limit */

/* Todo: Need to implement all-on-one-image output */

int main(int argc, char ** argv) {
  // Input variables
  // --------------
  // Filenames
  char *outfile = NULL, *datafile = NULL, *basedatafile = NULL, *mapfile = NULL, *slmaskfile = NULL, *latfile = NULL, *lonfile = NULL;
  // --------------
  // Misc text
  char * leg_text = "", *decalfile = NULL, *credit_text="Created by genimage v1.0", *identify_text="";
  // --------------
  // Flags
  int grid = 0, plot_over_ocean = 0, colour_map = 0, colour_map_rev = 0, range_dynamic = 0, image_type = 0;
  // --------------
  // Map data
  int img_width = -1, img_height = -1;
  double offset_lon = 0, offset_lat = 0;
  int zoom_factor = 1;
  // --------------
  // Range data
  double range_min = 0, range_max = 0;
  // --------------
  // Dec places for Legend data; and lat/lon data
  int dat_dec_places = 1, leg_dec_places = 1, lon_text_spacing = 20, lat_text_spacing = 10;
  // --------------
  // Other data
  int pct_change_data = 0;
  int timeofyear = 16, plot_type = 0;
  // --------------
  // List of points (relative to offset)
  Point* temp;
  Point** points;
  list<Point* > lpoints;
  list<Point* >::const_iterator lpoints_iter;
  
  FILE* infile, *infile2;
  
  // Point to get specific data about
  bool point_present = false;
  Point dpoint;

  int num_leg_blocks;
  int num_leg_segments;
  double leg_range_min;
  double leg_range_max;
  int colour;
  int i, j;
  int x, y;
  //double xf, yf;
  int xend, yend;
  int cols, rows;
  int dl_length, headerlen;
  int numpoints;
  double maxlat, minlat;
  double maxlon, minlon;
  double difflat, difflon, diffrange;
  double max_data, min_data;
  
  double* data;
  double* lats;
  double* lons;
  double* grid_lats;
  double* grid_lons;
  double* r_grid;
  int* xpoint;
  int* ypoint;
  
  char * dataline;
  char buf[1024];
  char buf2[128];
  int brect[8];
  char * ptr = buf;
  legend* leg_colours;
  range data_range;
  int imgsize;
  int** imgptr;
  int* lineptr;
  unsigned int* data_colours, *data_colours_ptr;
  int* slmask;
  int* data_mask, *draw_mask;
  char* font;
  gdImagePtr img;

  int xoffset, yoffset;
  int map_width, map_height;
  int leg_width, leg_height;
  int lat_width, lat_height;
  int map_offset_x, map_offset_y;
  int leg_offset_x, leg_offset_y;
  int lat_offset_x, lat_offset_y;
  int val = 0;

  char* optstring = "o:f:d:b:m:s:t:n:l:s:g:c:w:h:u:v:z:a:e:q:i:r:k:y:p:x:j:";
  static struct option options[] = {
    { "show-grid", 0, &grid, 1},
    { "colour-map-inverse", 0, &colour_map_rev, 1},
    { "dynamic-range", 0, &range_dynamic, 1},
    { "percent-change-calculations", 0, &pct_change_data, 1},
    { "image-type", 1, 0, '0'},
    { "decal-file", 1, 0, '1'},
    { "output-file", 1, 0, 'o'},
    { "font-file", 1, 0, 'f'},
    { "data-file", 1, 0, 'd'},
    { "base-data-file", 1, 0, 'b'},
    { "map-file", 1, 0, 'm'},
    { "slmask-file", 1, 0, 's'},
    { "lat-file", 1, 0, 't'},
    { "lon-file", 1, 0, 'n'},
    { "legend-text", 1, 0, 'l'},
    { "slmask-file", 1, 0, 's'},
    { "ocean-plot", 1, 0, 'g'},
    { "colour-map", 1, 0, 'c'},
    { "color-map", 1, 0, 'c'},
    { "map-width", 1, 0, 'w'},
    { "map-height", 1, 0, 'h'},
    { "lon-offset", 1, 0, 'u'},
    { "lat-offset", 1, 0, 'v'},
    { "zoom-factor", 1, 0, 'z'},
    { "range-min", 1, 0, 'a'},
    { "range-max", 1, 0, 'e'},
    { "data-decimal-places", 1, 0, 'q'},
    { "legend-decimal-places", 1, 0, 'i'},
    { "lat-text-spacing", 1, 0, 'r'},
    { "lon-text-spacing", 1, 0, 'k'},
    { "timeofyear", 1, 0, 'y'},
    { "plot-type", 1, 0, 'p'},
    { "data-point", 1, 0, 'x'},
    { "poly-point", 1, 0, 'j'},
    { "credit-text", 1, 0, '2'},
    { "identify-text", 1, 0, '3'},
    { 0, 0, 0, 0 }
  };
  while((val = getopt_long(argc, argv, optstring, options, NULL)) != -1) {
    switch(val) {
    case 0:
      break;
    case '0':
      // image-type
      image_type = atoi(optarg);
      break;
    case '1':
      // decal-file
      decalfile = optarg;
      break;
    case '2':
      // credit-text
      credit_text = optarg;
      break;
    case '3':
      // identify-text
      identify_text = optarg;
      break;
    case 'a':
      // range-min
      range_min = atof(optarg);
      break;
    case 'b':
      // base-data-file
      basedatafile = optarg;
      break;
    case 'c':
      // colour-map
      colour_map = atoi(optarg);
      break;
    case 'd':
      // data-file
      datafile = optarg;
      break;
    case 'e':
      // range-max
      range_max = atof(optarg);
      break;
    case 'f':
      // font-file
      font = optarg;
      break;
    case 'g':
      // ocean-plot
      plot_over_ocean = atoi(optarg);
      break;
    case 'h':
      // map-height
      map_height = img_height = atoi(optarg);
      break;
    case 'i':
      // legend-decimal-places
      leg_dec_places = atoi(optarg);
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
      lon_text_spacing = atoi(optarg);
      break;
    case 'l':
      // legend-text
      leg_text = optarg;
      break;
    case 'm':
      // map-file
      mapfile = optarg;
      break;
    case 'n':
      // lon-file
      lonfile = optarg;
      break;
    case 'o':
      // output-file
      outfile = optarg;
      break;
    case 'p':
      // plot-type
      plot_type = atoi(optarg);
      break;
    case 'q':
      // data-decimal-places
      dat_dec_places = atoi(optarg);
      break;
    case 'r':
      // lat-text-spacing
      lat_text_spacing = atoi(optarg);
      break;
    case 's':
      // slmask-file
      slmaskfile = optarg;
      break;
    case 't':
      // lat-file
      latfile = optarg;
      break;
    case 'u':
      // lat-offset
      offset_lat = atof(optarg);
      break;
    case 'v':
      // lon-offset
      offset_lon = atof(optarg);
      break;
    case 'w':
      // lon-offset
      map_width = img_width = atoi(optarg);
      break;
    case 'x':
      // data-point
      if(parse_data_point(optarg, dpoint)) {
	point_present = 1;
      }
      break;
    case 'y':
      // timeofyear
      timeofyear = atoi(optarg);
      break;
    case 'z':
      // zoom-factor
      zoom_factor = atoi(optarg);
      break;
    }
  }

  // Make sure we have enough data to continue
  if(!(outfile && datafile && (basedatafile || !pct_change_data) && mapfile && slmaskfile && latfile && lonfile && (range_dynamic || (range_min || range_max)) && (!point_present || decalfile))) {
    fprintf(stderr, "Insufficient arguments\n");
    exit(1);
  }
  
  if(plot_type >= TYPE_INVALID) {
    fprintf(stderr, "Invalid plot type\n");
  }

  // Store up the data points
  numpoints = lpoints.size();
  if(numpoints) {
    points = new Point*[numpoints];
    i = 0;
    for(lpoints_iter = lpoints.begin(); lpoints_iter != lpoints.end(); lpoints_iter++) {
      points[i] = (*lpoints_iter);
      i++;
    }
  }
  
  infile = fopen(mapfile, "r");
  if(!infile) {
    exit(2);
  }
  gdImagePtr image2 = gdImageCreateFromPng(infile);
  if(!image2) {
    exit(2);
  }
  fclose(infile);
  
  if(img_width == -1 || img_height == -1) {
    img_width = image2->sx + 2;
    img_height = image2->sy + 2;
  }
  
  // Size adjustments for different image types
  map_width = img_width;
  map_height = img_height;
  leg_width = map_width + LAT_WIDTH + LEG_EXTRA_RWIDTH;
  leg_height = LEG_HEIGHT;
  lat_width = LAT_WIDTH;
  lat_height = map_height + LAT_EXTRA_HEIGHT;

  // Special case of TYPE_ALL plot (integrates 3 images into 1)
  if(plot_type == TYPE_ALL || plot_type == TYPE_REGIONONLY || plot_type == TYPE_TEXT) {
    img_width = leg_width;
    img_height = leg_height + map_height;
    map_offset_x = lat_width;
    map_offset_y = 0;
    leg_offset_x = 0;
    leg_offset_y = map_height;
    lat_offset_x = 0;
    lat_offset_y = 0;

    img = gdImageCreateTrueColor(img_width, img_height);

    // Clear a particular region of the map
    gdImageFilledRectangle(img, map_offset_x + map_width, map_offset_y, img_width, map_height, 0x00FFFFFF);

  } else {
    map_offset_x = 0;
    map_offset_y = 0;
    leg_offset_x = 0;
    leg_offset_y = 0;
    lat_offset_x = 0;
    lat_offset_y = 0;
  }

  // Set up some vars etc
  imgsize = img_width * img_height;

  // Open data file
  infile = fopen(datafile, "r");
  
  if(!infile) {
    exit(3);
  }

  if(fgets(buf, 1024, infile)) {
    headerlen = strlen(buf);
    ptr = buf;
    // Assume first line is parameters -- valid.
    for(i = 0; i < 4; i++) {
      ptr = seektoword(ptr);
    }
    
    // At this point, we're at the # of cols
    cols = readint(ptr, buf2);
    
    // Now seek to the next non-space
    ptr = seektoword(ptr);
    
    // Now load in # of rows
    rows = readint(ptr, buf);
  } else {
    exit(4);
  }

  // Now, set up the arrays
  data = new double[rows * cols];
  lats = new double[rows];
  lons = new double[cols + 1];
  grid_lats = new double[rows + 1];
  grid_lons = new double[cols + 1];
  r_grid = new double[rows * cols];
  draw_mask = new int[rows * cols];
  data_mask = new int[rows * cols];
  ypoint = new int[rows + 1];
  xpoint = new int[cols + 1];
  slmask = new int[rows * cols];

  // Data line length
  dl_length = cols * VALUE_LENGTH + 3;
  dataline = new char[dl_length];

  // Load in sea-land mask file
  infile2 = fopen(slmaskfile, "r");
  if(!infile2) {
    exit(9);
  }
  //fgets(dataline, dl_length, infile2);
  load_grid(infile2, rows, cols, dataline, dl_length, slmask, SLMASK_LENGTH);
  fclose(infile2);

  // If we're getting rid of the land
  if(plot_over_ocean == 2) {
    for(i = 0; i < rows * cols; i++) {
      draw_mask[i] = !slmask[i];
      data_mask[i] = !slmask[i];
    }
  } else if(plot_over_ocean == 1) {
    for(i = 0; i < rows * cols; i++) {
      draw_mask[i] = 1;
      data_mask[i] = 1;
    }
  } else if(plot_over_ocean == 0) {
    for(i = 0; i < rows * cols; i++) {
      draw_mask[i] = slmask[i];
      data_mask[i] = slmask[i];
    }
  }
  
  // Seek to the data block we want in the data file
  if(timeofyear > 0) {
    if(fgets(dataline, dl_length, infile)) {
      // Calculate offset to beginning of data
      i = (dl_length - 1) * rows * timeofyear + (timeofyear + 1)  * headerlen;
      // Seek to offset (not relative)
      if(fseek(infile, i, SEEK_SET) == -1) {
	exit(5);
      }
    }
  }

  // Load up the data grid
  load_grid(infile, rows, cols, dataline, dl_length, data, VALUE_LENGTH, &data_range, data_mask, 1);
  fclose(infile);

  max_data = data_range.max();
  min_data = data_range.min();
  
  // Dynamic ranging switch
  // Do dynamic ranging before selection made

  // Generate the colour map
  leg_colours = new legend(range_min, range_max, colour_map, colour_map_rev);
  num_leg_segments = 10;
  num_leg_blocks = leg_colours->numcolours();
  diffrange = range_max - range_min;

  if(range_dynamic) {
    if(colour_map == 1) {
      num_leg_segments = 12;
    }
    //data_range.round(range_round_factor);
    range_min = data_range.min();
    range_max = data_range.max();
    leg_range_min = range_min;
    leg_range_max = range_max;
    leg_colours->setmin(range_min);
    leg_colours->setmax(range_max);
  } else {
    leg_range_max = range_max;
    leg_range_min = range_min;
    leg_colours->setmax(range_max);
    leg_colours->setmin(range_min);
    range_max += diffrange / (double)(num_leg_blocks - 2);
    range_min -= diffrange / (double)(num_leg_blocks - 2);
  }

  // Auto set the legend's decimal places
  leg_dec_places = range_dynamic;
  double tmp = (leg_range_max - leg_range_min) / (num_leg_segments);
  while(tmp < 1) {
    leg_dec_places++;
    tmp *= 10;
  }

  // Load the lat-lon grids (err, well, part of them. more about that below)
  infile = fopen(latfile, "r");
  if(!infile) {
    exit(6);
  }
  load_grid(infile, rows, 1, dataline, dl_length, lats, LL_LENGTH);
  fclose(infile);
  
  infile = fopen(lonfile, "r");
  if(!infile) {
    exit(7);
  }
  load_grid(infile, 1, cols, dataline, dl_length, lons, LL_LENGTH);
  fclose(infile);

  // Assumption! BEWARE!
  // This code assumes that the grid is essentially rectangular
  // This code also needs to take into account zoomed images

  grid_lats[0] = lats[0] - ((lats[1] - lats[0]) / 2);
  for(i = 1; i < rows; i++) {
    grid_lats[i] = lats[i] - ((lats[i] - lats[i - 1]) / 2);
  }
  grid_lats[rows] = lats[rows - 1] - ((lats[rows - 2] - lats[rows - 1]) / 2);

  maxlat = lats[0] - ((lats[1] - lats[0]) / 2);
  minlat = lats[rows - 1] + ((lats[rows - 1] - lats[rows - 2]) / 2);

  grid_lons[0] = lons[0] - ((lons[1] - lons[0]) / 2);
  for(i = 1; i < cols; i++) {
    grid_lons[i] = lons[i] - ((lons[i] - lons[i - 1]) / 2);
  }
  grid_lons[cols] = lons[cols - 1] + ((lons[cols - 1] - lons[cols - 2]) / 2);

  maxlon = lons[cols - 1] + ((lons[cols - 1] - lons[cols - 2]) / 2);
  minlon = lons[0] - ((lons[1] - lons[0]) / 2);

  // Make it all positive
  difflat = maxlat - minlat;
  difflon = maxlon - minlon;

  // This is bizarre. FIXME
  xoffset = 1;
  yoffset = 1;

  // Manually generate the top edge
  ypoint[0] = yoffset;
  // Now generate the body
  for(y = 1; y < rows; y++) {
    ypoint[y] = (int)(((maxlat - (lats[y] - ((lats[y] - lats[y - 1]) / 2))) / difflat) * (double)(map_height - BORDER_WIDTH*2) + 0.5) + yoffset;
  }
  // Manually generate the bottom edge
  ypoint[rows] = map_height - (BORDER_WIDTH*2) + yoffset;

  // Manually generate the left edge
  xpoint[0] = xoffset;
  // Now generate the body
  for(x = 1; x < cols; x++) {
    xpoint[x] = (int)(((lons[x] - ((lons[x] - lons[x - 1]) / 2) - minlon) / difflon) * (double)(map_width - BORDER_WIDTH*2) + 0.5) + xoffset;
  }
  // Manually generate the right edge
  xpoint[cols] = map_width - (BORDER_WIDTH*2) + xoffset;
  
  // Zero the grid
  for(i = 0; i < rows*cols; i++) {
    r_grid[i] = 0;
  }
  draw_polygon(rows, cols, r_grid, grid_lats, grid_lons, points, numpoints);

  // Mask off all boxes which are outside the poly in the data (stipple) mask
  // Need to allow most-covered grid box to be selected if no others available
  for(i = 0; i < rows; i++) {
    for(j = 0; j < cols; j++) {
      double area = (grid_lats[i] - grid_lats[i + 1]) * (grid_lons[j + 1] - grid_lons[j]);

      // If coverage is less than 30%, throw this box out
      //if(r_grid[(i * cols) + j] / area < 0.30) {
      if(r_grid[(i * cols) + j] / area <= 0.10) {
	data_mask[(i * cols) + j] = 0;
      } else {
	
      }
    }
  }

  if(plot_type == TYPE_MAP || plot_type == TYPE_ALL) {
    if(plot_type == TYPE_MAP) {
      img = gdImageCreateTrueColor(map_width, map_height);
    }
    imgptr = img->tpixels;
    data_colours = new unsigned int[rows * cols];

    unsigned char** map = image2->pixels;
    unsigned char* mlineptr;

    // Preprocess the data array with mask file
    for(i = 0; i < (rows * cols); i++) {
      // See comment below -- it applies here, in a slightly different form
      data_colours[i] = leg_colours->lookup(data[i]) | (~(((int)draw_mask[i] ^ 1) - 1) & 0x00FFFFFF);
    }
    
    // There is slower and easier to understand code for this in snippets.h
    data_colours_ptr = data_colours;
    y = ypoint[0] + map_offset_y;
    for(j = 0; j < rows; j++) {
      yend = ypoint[j + 1] + map_offset_y;
      for(; y < yend; y++) {
	lineptr = imgptr[y];
	mlineptr = map[y - yoffset - map_offset_y];
	data_colours_ptr = &data_colours[j * cols];
	x = xpoint[0] + map_offset_x;
	for(i = 0; i < cols; i++) {
	  // Look up colour
	  colour = data_colours_ptr[i];
	  xend = xpoint[i + 1] + map_offset_x;
	  for(; x < xend; x++) {
	    // Since map is either 0 or 1, when we subtract 1 we get either
	    // -1 or 0. 2's complement of -1 is 0xFFFFFFFF
	    // Now, this isn't quite what we want, so we invert the mask, 
	    // causing (for zero case) a mask of 0x0000000. Then we AND it 
	    // against the colour. This will zero any bits which are 1; aka, 
	    // make black.
	    if(numpoints > 2 && !data_mask[(j * cols) + i] && ((y + x) & 1)) {
	    //if(numpoints > 2 && !data_mask[(j * cols) + i] && (y & 1) && (x & 1)) {
	      lineptr[x] = 0;
	    } else {
	      lineptr[x] = colour & ~((int)mlineptr[x - xoffset - map_offset_x] - 1);
	    }
	  }
	}
      }
    }
    delete[] data_colours;
  }
  
  if(plot_type == TYPE_MAP || plot_type == TYPE_ALL || plot_type == TYPE_REGIONONLY) {
    if(plot_type == TYPE_REGIONONLY) {
      gdImageCopy(img, image2, map_offset_x + 1, map_offset_y + 1, 0, 0, image2->sx, image2->sy);
      imgptr = img->tpixels;
    }
    // After we're done with the map, we put the grid on it
    colour = 0x00A0A0A0;
    if(grid) {
      for(j = 0; j < rows; j++) {
	hline(imgptr, 0 + map_offset_x, map_width + map_offset_x, ypoint[j] + map_offset_y, map_width, colour);
      }
      hline(imgptr, 0 + map_offset_x, map_width + map_offset_x, ypoint[rows] - 1 + map_offset_y, map_width, colour);
      for(j = 0; j < cols; j++) {
	vline(imgptr, xpoint[j] + map_offset_x, 0 + map_offset_y, map_height + map_offset_y, map_width, colour);
      }
      vline(imgptr, xpoint[cols] - 1 + map_offset_x, 0 + map_offset_y, map_height + map_offset_y, map_width, colour);
    }
    
    colour = 0x00000000;

    // Lines to go on graph
    // This code is acceptable now
    gdImageSetAntiAliased(img, colour);
    gdImageSetThickness(img, LINE_WIDTH);
    int x1, x2, y1, y2;
    for(i = 0; i < numpoints; i++) {
      // Convert to XY from latlon
      x1 = lontox(map_width, offset_lon, minlon, maxlon, points[i]->x);
      y1 = lattoy(map_height, offset_lat, minlat, maxlat, points[i]->y);

      x2 = lontox(map_width, offset_lon, minlon, maxlon, points[(i + 1) % numpoints]->x);
      y2 = lattoy(map_height, offset_lat, minlat, maxlat, points[(i + 1) % numpoints]->y);

      // Lines
      gdImageLine(img, x1 + map_offset_x, y1 + map_offset_y, x2 + map_offset_x, y2 + map_offset_y, colour);
    }
    gdImageSetThickness(img, 1);
    for(i = 0; i < numpoints; i++) {
      // Convert to XY from latlon
      x1 = lontox(map_width, offset_lon, minlon, maxlon, points[i]->x);
      y1 = lattoy(map_height, offset_lat, minlat, maxlat, points[i]->y);

      // Point
      if(points[i]->selected) {
	gdImageFilledRectangle(img, x1 - (POINT_SIZE - 2) + map_offset_x, y1 - (POINT_SIZE - 2) + map_offset_y, x1 + (POINT_SIZE - 2) + map_offset_x, y1 + (POINT_SIZE - 2) + map_offset_y, 0x00FF0000);
	
      } else {
	gdImageFilledRectangle(img, x1 - 3 + map_offset_x, y1 - 3 + map_offset_y, x1 + 3 + map_offset_x, y1 + 3 + map_offset_y, 0x00000000);
      }
      
      // Box around point
      gdImageRectangle(img, x1 - (POINT_SIZE - 1) + map_offset_x, y1 - (POINT_SIZE - 1) + map_offset_y, x1 + (POINT_SIZE - 1) + map_offset_x, y1 + (POINT_SIZE - 1) + map_offset_y, 0x00FFFFFF);
      gdImageRectangle(img, x1 - POINT_SIZE + map_offset_x, y1 - POINT_SIZE + map_offset_y, x1 + POINT_SIZE + map_offset_x, y1 + POINT_SIZE + map_offset_y, 0x00000000);
    }      

    if(point_present) {
      // Display the data point
      // First, load the decal image (erroring if DNE)
      infile = fopen(decalfile, "r");
      if(!infile) {
	exit(66);
      }
      gdImagePtr decalimg = gdImageCreateFromPng(infile);
      if(!image2) {
	exit(66);
      }
      fclose(infile);
      
      // Next apply the decal image to the right spot
      gdImageCopy(img, decalimg, 
		  lontox(map_width, offset_lon, minlon, maxlon, dpoint.x) - 
		  ((decalimg->sx - 1) >> 1) + map_offset_x, 
		  lattoy(map_height, offset_lat, minlat, maxlat, dpoint.y) - 
		  ((decalimg->sy - 1) >> 1) + map_offset_y, 
		  0, 0, decalimg->sx, decalimg->sy);
    }
  }

  if(plot_type == TYPE_LEGEND || plot_type == TYPE_ALL || plot_type == TYPE_REGIONONLY) {
    if(plot_type == TYPE_LEGEND) {
      img = gdImageCreateTrueColor(leg_width, leg_height);
    }
    imgptr = img->tpixels;
    
    int leg_left, leg_right;
    int leg_top, leg_bottom;
    int width;
    double lbl_offset;
    double factor;
    double temp;
    double temp2;

    leg_left = LAT_WIDTH + BORDER_WIDTH + leg_offset_x;
    leg_right = leg_width - LEG_EXTRA_RWIDTH - BORDER_WIDTH + leg_offset_x;

    leg_top = LEG_TOP_TEXT_HEIGHT + leg_offset_y;
    leg_bottom = leg_height - LEG_BOTTOM_TEXT_HEIGHT + leg_offset_y;

    width = leg_right - leg_left;

    if(plot_type == TYPE_REGIONONLY) {
      gdImageFilledRectangle(img, 0 + leg_offset_x, 0 + leg_offset_y, leg_width + leg_offset_x, leg_height + leg_offset_y, 0x00FFFFFF);
    } else {
      i = leg_colours->numcolours();
      
      factor = (double)width / (double)i;
      
      // White out the white areas
      gdImageFilledRectangle(img, 0 + leg_offset_x, 0 + leg_offset_y, leg_width + leg_offset_x, leg_top - 2, 0x00FFFFFF);
      gdImageFilledRectangle(img, 0 + leg_offset_x, leg_top - 2, leg_left - 2, leg_bottom + 2, 0x00FFFFFF);
      gdImageFilledRectangle(img, leg_right + 2, leg_top - 2, leg_width, leg_bottom + 2, 0x00FFFFFF);
      gdImageFilledRectangle(img, 0 + leg_offset_x, leg_bottom + 2, leg_width + leg_offset_x, leg_height + leg_offset_y, 0x00FFFFFF);
      
      // Draw the surrounding box
      gdImageLine(img, leg_left - 1, leg_top - 1, leg_right + 1, leg_top - 1, 0x00000000);
      gdImageLine(img, leg_right + 1, leg_top - 1, leg_right + 1, leg_bottom + 1, 0x00000000);
      gdImageLine(img, leg_right + 1, leg_bottom + 1, leg_left - 1, leg_bottom + 1, 0x00000000);
      gdImageLine(img, leg_left - 1, leg_bottom + 1, leg_left - 1, leg_top - 1, 0x00000000);

      // Legend text
      // Get ready to render text
      colour = 0x00000000;
      x = leg_left + width / 2;
      y = leg_offset_y + leg_height;
      gdImageStringFT(NULL, brect, colour, font, 12, 0, x, y, leg_text);
      
      // Center the text horizontally and align it upwards vertically
      x -= (brect[2] - brect[0]) >> 1;
      y -= (brect[3] - (leg_offset_y + leg_height)) + 15;
      
      // Render
      gdImageStringFT(img, brect, colour, font, 12, 0, x, y, leg_text);
      
      // Draw the legend
      for(i = leg_colours->numcolours() - 1; (i + 1); i--) {
	colour = leg_colours->lookup(i);
	gdImageFilledRectangle(img, (int)round(i * factor) + leg_left, leg_top, (int)round((i + 1) * factor) + leg_left, leg_bottom, colour);
      }
      
      // Draw on the legend and give it text
      if(range_dynamic) {
	lbl_offset = 0;
	factor = (double)(width) / (num_leg_segments);
      } else {
	lbl_offset = factor;
	factor = (double)(width - (2 * lbl_offset)) / num_leg_segments;
      } 
      
      // Create the tick marks in the legend
      colour = 0x00000000;
      for(i = 0; i <= num_leg_segments; i++) {
	x = leg_left + (int)round(i * factor + lbl_offset);
	gdImageLine(img, x, leg_top, x, leg_top + 4, colour);
	gdImageLine(img, x, leg_bottom - 4, x, leg_bottom, colour);
      }
      
      // Label the ticks
      temp2 = ((leg_range_max - leg_range_min) / num_leg_segments);
      for(i = 0;i <= num_leg_segments; i++) {
	x = leg_left + (int)round(i * factor + lbl_offset);
	y = leg_bottom;
	temp = temp2 * i + leg_range_min;
	
	sprintf(buf2, "%%.%if", leg_dec_places);
	sprintf(buf, buf2, temp);
	
	gdImageStringFT(NULL, brect, colour, font, 10, 0, x, y, buf);
	
	// Center the text horizontally and align it downwards vertically
	x -= (brect[2] - brect[0]) >> 1;
	y += (leg_bottom - brect[5]) + 5;
	
	// Render
	gdImageStringFT(img, brect, colour, font, 10, 0, x, y, buf);
      }

      // Identify text
      // Get ready to render text
      colour = 0x00000000;
      x = leg_offset_x;
      y = leg_offset_y + leg_height;
      gdImageStringFT(NULL, brect, colour, font, 8, 0, x, y, identify_text);
      
      // Align text left and up
      x -= brect[0];
      y -= (brect[3] - (leg_offset_y + leg_height));
      
      //fprintf(stderr, "x: %i, y: %i\n", x, y);
      
      // Render
      gdImageStringFT(img, brect, colour, font, 8, 0, x, y, identify_text);
    }

    // Credit text
    // Get ready to render text
    colour = 0x00000000;
    x = leg_offset_x + leg_width;
    y = leg_offset_y + leg_height;
    gdImageStringFT(NULL, brect, colour, font, 10, 0, x, y, credit_text);
    
    // Align text right and up
    x -= (brect[2] - (leg_offset_x + leg_width));
    y -= (brect[3] - (leg_offset_y + leg_height));
    
    //fprintf(stderr, "x: %i, y: %i\n", x, y);
    
    // Render
    gdImageStringFT(img, brect, colour, font, 10, 0, x, y, credit_text);
    
    // Lon grid
    factor = width / difflon;
    
    j = lon_text_spacing;
    for(i = (int)floorf(maxlon / j) * j; i > minlon; i -= j) {
      x = leg_left + (int)round((i - minlon) * factor);
      y = leg_offset_y;

      // Put bottom tick marks on
      gdImageLine(img, x, y, x, y + 3, 0x00000000);

      // Format the lon string nicely
      if(i < 0) {
	sprintf(buf, "%iW", abs(i));
      } else if(i == 0) {
	sprintf(buf, "%i", i);
      } else if(i > 0) {
	sprintf(buf, "%iE", i);
      }

      gdImageStringFT(NULL, brect, colour, font, 10, 0, x, y, buf);
      
      // Center the text horizontally and align it downwards vertically
      x -= (brect[2] - brect[0]) >> 1;
      y -= (brect[5] - leg_offset_y) - 6;
      
      // Render
      gdImageStringFT(img, brect, colour, font, 10, 0, x, y, buf);

    }
  }

  if(plot_type == TYPE_LAT || plot_type == TYPE_ALL || plot_type == TYPE_REGIONONLY) {
    double factor;
    if(plot_type == TYPE_LAT) {
      img = gdImageCreateTrueColor(lat_width, lat_height);
    }
    imgptr = img->tpixels;
    gdImageFilledRectangle(img, 0 + lat_offset_x, 0 + lat_offset_y, lat_width + lat_offset_x - 1, lat_height + lat_offset_y - 1, 0x00FFFFFF);

    
    // Lat points
    factor = lat_height / difflat;
    colour = 0x00000000;

    j = lat_text_spacing;
    for(i = (int)floorf(maxlat / j) * j; i > minlat; i -= j) {
      x = lat_offset_x;
      y = lat_offset_y + (lat_height - (int)round((i - minlat) * factor));

      gdImageLine(img, lat_width - 4, y, lat_width - 1, y, 0x00000000);

      // Format the lon string nicely
      if(i < 0) {
	sprintf(buf, "%iS", abs(i));
      } else if(i == 0) {
	sprintf(buf, "%i", i);
      } else if(i > 0) {
	sprintf(buf, "%iN", i);
      }

      gdImageStringFT(NULL, brect, colour, font, 10, 0, x, y, buf);
      
      // Center the text horizontally and align it downwards vertically
      x = lat_width - brect[2] - 6;
      y -= (brect[5] - brect[1]) >> 1;
      
      if(x >= 0 && x < lat_width && y >= 0 && y < lat_height) {
	// Render
	gdImageStringFT(img, brect, colour, font, 10, 0, x, y, buf);
      }
    }
  }
  FILE * outfd;

  if(strcmp(outfile, "-") == 0) {
    outfd = stdout;
  } else {
    int error;
    outfd = fopen(outfile, "wb");

    if(!outfd) {
      fprintf(stderr, "Couldn't open file %s for writing!\n", outfile);
      exit(12);
    }

    // Try to lock the file
    error = flock(fileno(outfd), LOCK_EX|LOCK_NB);

    // Error handling
    if(error == EBADF || error == EINTR || error == EINVAL) {
      fprintf(stderr, "Something bad happened while trying to lock file %s!\n", outfile);
      exit(12);
    } else if(error == EWOULDBLOCK) {
      // File is open by somebody else
      fclose(outfd);
      exit(99);
    }
  }

  if(plot_type == TYPE_PLOTINFO) {
    fprintf(outfd, "Lat range: (%0.2f)-(%0.2f)\n", minlat, maxlat);
    fprintf(outfd, "Lon range: (%0.2f)-(%0.2f)\n", minlon, maxlon);
    fprintf(outfd, "Size: (%i)-(%i)\n", img_width, img_height);
    fprintf(outfd, "Map size: (%i)-(%i)\n", map_width, map_height);
    fprintf(outfd, "Map offset: (%i)-(%i)\n", map_offset_x, map_offset_y);
  }

  if(plot_type == TYPE_TEXT) {
    double* pred_data = new double[rows * cols];
    double* base_data = new double[rows * cols];

    double map_data_max = -INFINITY;
    double map_data_max_lat = 0;
    double map_data_max_lon = 0;
    double map_data_min = INFINITY;
    double map_data_min_lat = 0;
    double map_data_min_lon = 0;
    for(i = 0; i < rows; i++) {
      for(j = 0; j < cols; j++) {
	if(draw_mask[(i * cols) + j]) {
	  if(data[(i * cols) + j] > map_data_max) {
	    map_data_max = data[(i * cols) + j];
	    map_data_max_lon = lons[j];
	    map_data_max_lat = lats[i];
	  } 
	  if(data[(i * cols) + j] < map_data_min) {
	    map_data_min = data[(i * cols) + j];
	    map_data_min_lon = lons[j];
	    map_data_min_lat = lats[i];
	  }
	}
      }
    }
    // Data analysis loop stuff
    if(numpoints > 2) {
      list<WPoint> wpoints;
      list<WPoint> pwpoints;
      list<WPoint>::const_iterator p_iter;
      double base_data_sum = 0;
      double base_data_avg = 0;
      double total_weight = 0;
      double total_squared_weight = 0;
      double data_avg = 0;
      double data_sum = 0;
      double data_min = INFINITY;
      double data_min_lon = 0;
      double data_min_lat = 0;
      double data_max = -INFINITY;
      double data_max_lon = 0;
      double data_max_lat = 0;
      double data_variance = 0;
      double data_stddev = 0;
      double data_median = 0;
      double data_wmedian = 0;

      double* dataptr = data;

      if(pct_change_data) {
	// Open baseline data file
	infile = fopen(basedatafile, "r");
	
	if(!infile) {
	  exit(13);
	}
	
	if(fgets(buf, 1024, infile)) {
	  headerlen = strlen(buf);
	  ptr = buf;
	  // Assume first line is parameters -- valid.
	  for(i = 0; i < 4; i++) {
	    ptr = seektoword(ptr);
	  }
	  
	  // At this point, we're at the # of cols
	  // Check and make sure that baseline and prediction have same # cols
	  if(cols != readint(ptr, buf2)) {
	    exit(14);
	  }
	  
	  // Now seek to the next non-space
	  ptr = seektoword(ptr);
	  
	  // Now load in # of rows
	  // Check and make sure that baseline and prediction have same # rows
	  if(rows != readint(ptr, buf2)) {
	    exit(15);
	  }
	} else {
	  exit(16);
	}
	
	// Seek to the data block we want in the data file
	if(timeofyear > 0) {
	  if(fgets(dataline, dl_length, infile)) {
	    // Calculate offset to beginning of data
	    i = (dl_length - 1) * rows * timeofyear + (timeofyear + 1)  * headerlen;
	    // Seek to offset (not relative)
	    if(fseek(infile, i, SEEK_SET) == -1) {
	      exit(17);
	    }
	  }
	}
	
	// Load up the baseline data grid
	load_grid(infile, rows, cols, dataline, dl_length, base_data, VALUE_LENGTH);
	fclose(infile);

	// Correct the baseline data grid for the % change
	for(i = 0; i < rows; i++) {
	  for(j = 0; j < cols; j++) {
	    pred_data[(i * cols) + j] = base_data[(i * cols) + j] + base_data[(i * cols) + j] * (data[(i * cols) + j] / 100);
	  }
	}
	dataptr = pred_data;
      }

      for(i = 0; i < rows; i++) {
	for(j = 0; j < cols; j++) {
	  if(data_mask[(i * cols) + j] == 1) {
	    double area = (M_PI / 180) * squared(EARTH_RADIUS) * fabs(sin(grid_lats[i] * (M_PI / 180)) - sin(grid_lats[i + 1] * (M_PI / 180))) * fabs(grid_lons[j] - grid_lons[j + 1]);

	    //double area = cos(lats[i] * (M_PI / 180));
	    data_sum += area * dataptr[(i * cols) + j];
	    total_squared_weight += squared(area);
	    total_weight += area;
	    if(data[(i * cols) + j] > data_max) {
	      data_max = data[(i * cols) + j];
	      data_max_lon = lons[j];
	      data_max_lat = lats[i];
	    } 
	    if(data[(i * cols) + j] < data_min) {
	      data_min = data[(i * cols) + j];
	      data_min_lon = lons[j];
	      data_min_lat = lats[i];
	    }
	    wpoints.push_back(WPoint(data[(i * cols) + j], area));
	    if(pct_change_data) {
	      base_data_sum += area * base_data[(i * cols) + j];
	      pwpoints.push_back(WPoint(dataptr[(i * cols) + j], area));
	    }
	  }
	}
      }

      data_avg = data_sum / total_weight;
      
      if(pct_change_data) {
	base_data_avg = base_data_sum / total_weight;
      }


      // EVIL HACK: CHANGE THIS BUG BUG BUG
      if(pct_change_data && 0) {
	// Go through the list of pts and weights getting variance components
	for(p_iter = pwpoints.begin(); p_iter != pwpoints.end(); p_iter++) {
	  data_variance += squared((*p_iter).p - data_avg) * (*p_iter).w;
	}
      } else {
	// Go through the list of pts and weights getting variance components
	for(p_iter = wpoints.begin(); p_iter != wpoints.end(); p_iter++) {
	  data_variance += squared((*p_iter).p - data_avg) * (*p_iter).w;
	}
      }
      // Calculate variance
      // See http://pygsl.sourceforge.net/reference/pygsl/node35.html
      data_variance *= total_weight / (squared(total_weight) - total_squared_weight);

      // Calculate standard deviation
      data_stddev = sqrtf(data_variance);

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
      //desired_weight = wpoints.size() / 2;
      wsum = 0;
      for(p_iter = wpoints.begin(); p_iter != wpoints.end() && wsum <= desired_weight; p_iter++) {
	last_pt = curr_pt;
	curr_pt = *p_iter;
	wsum++;
      }

      if(wsum - desired_weight < 1) {
	data_median = curr_pt.p;
      } else {
	// Average
	data_median = (last_pt.p + curr_pt.p) / 2;
      }

      if(pct_change_data) {
	data_avg = ((data_avg - base_data_avg) / base_data_avg) * 100;
	//data_stddev = (data_stddev / base_data_avg) * 100;
      }

      sprintf(buf2, "Selection area weighted mean: (%%0.%if)\n", dat_dec_places);
      fprintf(outfd, buf2, data_avg);

      sprintf(buf2, "Selection area weighted median: (%%0.%if)\n", dat_dec_places);
      fprintf(outfd, buf2, data_wmedian);


      sprintf(buf2, "Selection area weighted standard deviation: (%%0.%if)\n", dat_dec_places);
      fprintf(outfd, buf2, data_stddev);

      sprintf(buf2, "Selection median: (%%0.%if)\n", dat_dec_places);
      fprintf(outfd, buf2, data_median);

      sprintf(buf2, "Selection data min (d, lon, lat): (%%0.%if)-(%%0.2f)-(%%0.2f)\n", dat_dec_places);
      fprintf(outfd, buf2, data_min, data_min_lon, data_min_lat);

      sprintf(buf2, "Selection data max (d, lon, lat): (%%0.%if)-(%%0.2f)-(%%0.2f)\n", dat_dec_places);
      fprintf(outfd, buf2, data_max, data_max_lon, data_max_lat);

      fprintf(outfd, "Selection area: (%.0f) km<sup>2</sup>\n", total_weight);
      fprintf(outfd, "Selection num grid boxes: (%i)\n", wpoints.size());
    }

    fprintf(outfd, "Lat range: (%0.2f)-(%0.2f)\n", minlat, maxlat);
    fprintf(outfd, "Lon range: (%0.2f)-(%0.2f)\n", minlon, maxlon);
    sprintf(buf2, "Data range: (%%0.%if)-(%%0.%if)\n", dat_dec_places, dat_dec_places);
    fprintf(outfd, buf2, min_data, max_data);
    fprintf(outfd, "Legend colour range: (%f)-(%f)\n", range_min, range_max);
    fprintf(outfd, "Size: (%i)-(%i)\n", img_width, img_height);
    fprintf(outfd, "Map size: (%i)-(%i)\n", map_width, map_height);
    fprintf(outfd, "Map offset: (%i)-(%i)\n", map_offset_x, map_offset_y);

    sprintf(buf2, "Map data min (d, lon, lat): (%%0.%if)-(%%0.2f)-(%%0.2f)\n", dat_dec_places);
    fprintf(outfd, buf2, map_data_min, map_data_min_lon, map_data_min_lat);
    
    sprintf(buf2, "Map data max (d, lon, lat): (%%0.%if)-(%%0.2f)-(%%0.2f)\n", dat_dec_places);
    fprintf(outfd, buf2, map_data_max, map_data_max_lon, map_data_max_lat);

    // If the user wants to know about the data box that a specific lat/lon 
    // is in, report back
    if(point_present) {
      // Find the column it's in
      for(i = 0; i < cols; i++) {
	if(grid_lons[i] < dpoint.x && grid_lons[i + 1] > dpoint.x) {
	  break;
	}
      }

      // Find the row it's in
      for(j = 0; j < rows; j++) {
	if(grid_lats[j + 1] < dpoint.y && grid_lats[j] > dpoint.y) {
	  break;
	}
      }

      if(j != rows && i != cols) {
	fprintf(outfd, "Grid box: (%i)-(%i)\n", i, j);
      }
      sprintf(buf2, "Data point (d, lon, lat): (%%0.%if)-(%%0.2f)-(%%0.2f)\n", dat_dec_places);
      fprintf(outfd, buf2, data[(j * cols) + i], lons[i], lats[j]);
      //fprintf(outfd, buf2, data[(j * cols) + i], (double)i, (double)j);
    }
    delete[] pred_data;
    delete[] base_data;
  } else if(plot_type == TYPE_MASK) {
    for(i = 0; i < rows; i++) {
      for(j = 0; j < cols; j++) {
	fprintf(outfd, "% i", data_mask[(i * cols) + j]);
      }
      fprintf(outfd, "\r\n");
    }
  } else if(plot_type == TYPE_GEOREF) {
    // Time for georef. Ugh.
    // Open data file
    double* fulldata[NUM_TIMESLICES];

    infile = fopen(datafile, "r");
    
    if(!infile) {
      exit(18);
    }

    for(i = 0; i < NUM_TIMESLICES; i++) {
      fulldata[i] = new double[rows * cols];
      if(!fgets(buf, 1024, infile)) {
	exit(19);
      }
    
      load_grid(infile, rows, cols, dataline, dl_length, fulldata[i], VALUE_LENGTH);
    }
    fclose(infile);

    for(i = 0; i < rows; i++) {
      for(j = 0; j < cols; j++) {
	if(data_mask[(i * cols) + j]) {
	  fprintf(outfd, " % #9.4f % #9.4f", lats[i], lons[j]);
	  for(int q = 0; q < NUM_TIMESLICES; q++) {
	    fprintf(outfd, " % E", fulldata[q][(i * cols) + j]);
	  }
	  fprintf(outfd, "\r\n");
	}
      }
    }

    for(i = 0; i < NUM_TIMESLICES; i++) {
      delete[] fulldata[i];
    }
  } else {
    if(image_type == IMAGE_TYPE_PNG) {
      gdImagePng(img, outfd);
    } else {
      gdImageJpeg(img, outfd, JPEG_QUALITY);
    }
    gdImageDestroy(img);
  }  

  if(strcmp(outfile, "-") != 0) {
    flock(fileno(outfd), LOCK_UN);
    fclose(outfd);
  }

  delete leg_colours;

  for(i = 0; i < numpoints; i++) {
    delete points[i];
  }
  if(numpoints) {
    delete[] points;
  }
  gdImageDestroy(image2);
  delete[] slmask;
  delete[] data_mask;
  delete[] draw_mask;
  delete[] xpoint;
  delete[] ypoint;
  delete[] data;
  delete[] lats;
  delete[] lons;
  delete[] grid_lats;
  delete[] grid_lons;
  delete[] r_grid;

  delete[] dataline;

  return 0;
}
