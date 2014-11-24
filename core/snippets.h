/*
bool sameside(const Point& p1, const Point& p2, const Line& l) {
  double crossproduct = (l.to - l.from).cross(p1 - l.from) * (l.to - l.from).cross(p2 - l.from);
  if(crossproduct >= 0) {
    return true;
  } else {
    return false;
  }
}
// Internal
// Swaps things in the swap array if they're not ordered
void swap_if_unordered(int* swap_array, 
		       Point points[], 
		       Point from, 
		       int to1, 
		       int to2)
{
  if(points[swap_array[to1]].valid() && points[swap_array[to2]].valid()) {
    if(from.dist_from(points[swap_array[to1]]) > from.dist_from(points[swap_array[to2]])) {
      swap(swap_array[to1], swap_array[to2]);
    }
  } else {
    if(!points[swap_array[to1]].valid()) {
      // Try to move NULLs to end
      swap(swap_array[to1], swap_array[to2]);
    }
  }
}

// Returns true if above line, false otherwise
bool right_of_line(const Point& p, const Line& l) {
  if(l.slope == 0) {
    return (((l.to.x - l.from.x) * (l.from.y - p.y)) > 0);
  } else if(l.slope == INFINITY || l.slope == -INFINITY) {
    return (((l.from.x - p.x) * (l.to.y - l.from.y)) < 0);
  } else {
    return (((((p.x - l.from.x) * l.slope) + l.from.y - p.y) * (l.to.y - l.from.y) * (l.to.x - l.from.x)) > 0);
  }
}

// Assumes rectangle
bool within_bbox(const Point& p, const Point& topright, const Point& bottomleft) {
  return ((bottomleft.x <= p.x) && (topright.x > p.x) && (bottomleft.y <= p.y) && (topright.y > p.y));
}

*/

/*
// Assumes points come in in clockwise order - THIS IS REQUIRED
void draw_triangle(int rows, int cols, double* grid, const double* lats, const double* lons, Point p1, Point p2, Point p3) {
  int i, j, k, l;
  //double maxlat = lats[0], minlat = lats[rows];
  //double maxlon = lons[cols], minlon = lons[0];
  Point topright, bottomleft;

  //cout << p1 << p2 << p3 << endl;

  Line line[] = { 
    Line(p1, p2, true), 
    Line(p2, p3, true), 
    Line(p3, p1, true)
  };

  int numlines = sizeof(line) / sizeof(Line);
  Point ipoints[numlines];
  Line* box;

  for(i = 0; i < rows; i++) {
    // Reverse the lats -- arg. I want them to go min -> max
    double tlat = lats[i * 2];
    double blat = lats[i * 2 + 1];

    //cout.setf(ios::fixed);
    //cout << tlat << "-" << blat << ": ";


    for(j = 0; j < cols; j++) {
      double llon = lons[j * 2];
      double rlon = lons[j * 2 + 1];
      topright = Point(rlon, tlat);
      bottomleft = Point(llon, blat);

      list<Point > outline;
      list<Point >::iterator ol_iter;

      box = get_box_outline(topright, bottomleft);
      
      // Find out if top left corner is inside
      bool inside = true;
      for(l = 0; l < numlines; l++) {
	//inside = inside && right_of_line(box[0].from, line[l]);
	inside = inside && sameside(box[0].from, line[l].from, line[(l + 1) % numlines]);
      }

      for(k = 0; k < QUAD_SIDES; k++) {
	int sar[] = { 0, 1, 2};
	for(l = 0; l < numlines; l++) {
	  ipoints[l] = ls_intersect(box[k], line[l]);
	}
	
	// Sort the intersects the right way using a translation array
	swap_if_unordered(sar, ipoints, box[k].from, 0, 1);
	swap_if_unordered(sar, ipoints, box[k].from, 0, 2);
	swap_if_unordered(sar, ipoints, box[k].from, 1, 2);
	
	for(l = 0; ipoints[sar[l]].valid(); l++) {
	  //for(l = 0; l < numlines; l++) {
	  if(ipoints[sar[l]].valid()) {
	    // Flip 'inside' bit
	    inside = !inside;
	    if(within_bbox(line[sar[l]].from, topright, bottomleft)) {
	      // Going out
	      outline.push_back(line[sar[l]].from);
	    }
	    outline.push_back(ipoints[sar[l]]);
	  }
	}
	if(inside) {
	  outline.push_back(box[k].to);
	}
      }

      if(outline.size() > 0) {
	if(outline.size() > 2) {
	  // Now we've supposedly generated the outline... get the area
	  // Triangle fan approach
	  ol_iter = outline.begin();
	  Point pivot, prev, cur;
	  pivot = *ol_iter;
	  
	  ol_iter++;
	  prev = *ol_iter;
	  
	  ol_iter++;
	  cur = *ol_iter;
	  while(ol_iter != outline.end()) {
	    // Add area of this triangle to the grid box
	    grid[i * cols + j] += triangle_area(pivot, prev, cur);
	    prev = cur;
	    ol_iter++;
	    cur = *ol_iter;
	  }
	} else {
	  fprintf(stderr, "%s:%i: Error: Can't have a poly with less than 3 sides\n", __FILE__, __LINE__);
	}
      }
      //cout.setf(ios::left);
      //cout.precision(3);
      //cout << "" << grid[i * cols + j];

      delete[] box;
      //cout << inside << " ";
      //cout << "Max: Lat point: " << blat << "\tLon point: " << rlon << endl;
    }
    //cout << endl;
  }
}


bool clockwise(int piv, int nxt, int prv, const vector<Point>& points) {
  return ((points[nxt] - points[piv]).cross(points[piv] - points[prv]) > 0);
}

bool line_crosses(int from, int to, const vector<Point>& points) {
  int nxt;
  Point intersect;
  const int numlines = points.size();
  for(int i = 0; i < numlines; i++) {
    nxt = (i + 1) % numlines;
    // Make sure that the lines do not include endpoints
    if(from == i || from == nxt) {
      // Do no checks
    } else {
      // Do not check endpoints
      if((from != i) && (from != nxt) && (to != i) && (to != nxt)) {
	// Check if the line segments intersect
	intersect = ls_intersect(points[from], points[to], points[i], points[nxt]);
	if(intersect.valid()) {
	  // If so, clean up and say they do
	  return true;
	}
      }
    }
  }
  return false;
}

// Tries to find a triangle originating at <whatever>
bool solve_connects(int idx, int points[], int connects[], int numpoints, int i_initial) {
  // If we are at the end of our run (base case for our recursion)
  if(idx == 2) {
    // Do final check
    if(connects[(points[0] * numpoints) + points[2]] >= LINE) {
      // If all is good, we have a solution
      return true;
    } else {
      return false;
    }
  }
  // If we aren't at end of run, go through the connect list looking for
  // a connection to this point
  for(int i = i_initial; i < numpoints; i++) {
    if(connects[(points[idx] * numpoints) + i] >= LINE) {
      points[idx + 1] = i;
      if(solve_connects(idx + 1, points, connects, numpoints, i + 1)) {
	return true;
      }
    }
  } // end for

  // If we reach here, the possibilities have been exhausted. There is no
  // solution.
  return false; 
}

// Tries to find a triangle originating at <whatever>
bool solve_connects(int idx, int points[], int connects[], int numpoints) {
  return solve_connects(idx, points, connects, numpoints, 0);
}


void remove_intercepts(int p1, int p2, const vector<Point>& points, int connects[]) {
  // If temp line...
  const int numpoints = points.size();
  if(connects[(p1 * numpoints) + p2] == TEMP) {
    for(int i = 0; i < numpoints; i++) {
      if(i == p1 || i == p2)
	continue;

      for(int j = 0; j < numpoints; j++) {
	if(j == p1 || j == p2)
	  continue;
	
	if(connects[(i * numpoints) + j] == TEMP && ls_intersect(points[p1], points[p2], points[i], points[j]).valid())
	  connects[(i * numpoints) + j] = USED;
      } // end for
    } // end for
  }
}

// Decrements this intersect and its equivalent at other end of line
void dec(int p1, int p2, int connects[], int numpoints) {
  if(p1 == p2) {
    // CANT HAPPEN
    assert(false);
  } else {
    connects[(p1 * numpoints) + p2]--;
    connects[(p2 * numpoints) + p1]--;
  }
}


// Adaptation of algorithm found at:
// http://www.ecse.rpi.edu/Homepages/wrf/research/geom/pnpoly.html
// Copyright (c) 1970-2003, Wm. Randolph Franklin
bool is_inside(Point p, const vector<Point>& points) {
  bool inside = false;
  int i, j;
  const int numpoints = points.size();
  for(i = 0, j = numpoints - 1; i < numpoints; j = i++) {
    if ((((points[i].y <= p.y) && (p.y < points[j].y)) ||
         ((points[j].y <= p.y) && (p.y < points[i].y))) &&
        (p.x < (points[j].x - points[i].x) * (p.y - points[i].y) / (points[j].y - points[i].y) + points[i].x)) {
      inside = !inside;
    }
  }
  return inside;
}


*/

// Assumes incoming shapes are POLYGONS -- NOT overlapping
// This code sucks.
/*
  Instead of splitting up by triangles in some fucked up way then computing overlap, why not:
   - Only check boxes within the bounding box
   - Within the bounding box, "draw" the lines of the polygon on the grid
     - Store a list of grid boxes by ID, along with which line touched them; O(length of polygon path)
     - Sort the final list by the ID and, if more than one line touched a grid box, merge the entries up into one entry; O(n log n)
     - For each box in this final list, intersect it with the lines that touch it and compute the coverage: O(num lines touching box * length of polygon path)
     - Simple way: For all other grid boxes, test if upper left corner of polygon is within polygon or not using is_inside
       - Continue using this value of "inside" without retesting until you reach an edge
 */
/*
void draw_polygon(int rows, int cols, double* grid, const double* lats, const double* lons, const vector<Point>& points) {
  const int numpoints = points.size();
  int connects[numpoints * numpoints];
  int* conn;
  int i, j;

  // Special case poly within gb
  Point center(0, 0);
  for(i = 0; i < numpoints; i++) {
    center.x += points[i].x;
    center.y += points[i].y;
  }
  center.x /= numpoints;
  center.y /= numpoints;

  int xgb, ygb;
  if(find_in_grid_var(center.x, cols, lons, xgb, true) &&
     find_in_grid_var(center.y, rows, lats, ygb, false)) {
    bool all_inside = true;
    for(i = 0; i < numpoints; i++) {
      bool inside = ((lons[xgb * 2 + 1] >= points[i].x && points[i].x > lons[xgb * 2]) &&
		     (lats[ygb * 2] >= points[i].y && points[i].y > lats[ygb * 2 + 1]));
      all_inside &= inside;
    }

    if(all_inside) {
      grid[(ygb * cols) + xgb] = 1.0;
      return;
    }
  }
     

  // Clear table
  memset(connects, 0, numpoints*numpoints*sizeof(int));

  // Look for intersects
  for(i = 0; i < numpoints; i++) {
    conn = &connects[i * numpoints];
    for(j = 0; j < numpoints; j++) {
      // Don't include self-connects or initialized connects
      if(conn[j] != UNDEF) {
	continue;
      } 
      if(i == j) {
	connects[(i * numpoints) + j] = USED;
	continue;
      }

      // If this is an edge, it can be used once (we set type)
      if(j == (i + 1) % numpoints || i == (j + 1) % numpoints) {
	connects[(i * numpoints) + j] = LINE;
	connects[(j * numpoints) + i] = LINE;
	continue;
      }
      // If there are no intersecting real lines on this...
      if(!line_crosses(i, j, points)) {
	// If this is inside...
	if(is_inside((points[i] + points[j]) / 2, points)) {
	  // It may be an internal line. Can be used twice.
	  connects[(i * numpoints) + j] = TEMP;
	  connects[(j * numpoints) + i] = TEMP;
	  continue;
	}
      } else {
	// Something crossed this would-be internal line -- throw it out
	connects[(i * numpoints) + j] = USED;
	connects[(j * numpoints) + i] = USED;
      }
    }
  }

  // Look for intersects
  // WARNING POSSIBLE INFINITE LOOP
  int pts[3];
  pts[0] = 0;
  while(pts[0] < numpoints) {
    // Tries to find _some_ connection between points
    if(solve_connects(0, pts, connects, numpoints)) {
      // Do intersect checking on TEMP lines for other TEMP lines
      remove_intercepts(pts[0], pts[1], points, connects);
      remove_intercepts(pts[1], pts[2], points, connects);
      remove_intercepts(pts[2], pts[0], points, connects);

      // Mark lines taken
      dec(pts[0], pts[1], connects, numpoints);
      dec(pts[1], pts[2], connects, numpoints);
      dec(pts[2], pts[0], connects, numpoints);
      
      //cerr << pts[0] << " " << pts[1] << " " << pts[2] << endl;
      //cerr << "Triangle at " << points[pts[0]] << "-" << points[pts[1]] << "-" << points[pts[2]] << endl;

      draw_triangle(rows, cols, grid, lats, lons, points[pts[0]], points[pts[1]], points[pts[2]]);
    } else {
      // Done with this point
      pts[0]++;
    }
  } // end while

} // end draw_polygon
double triangle_area(Point p1, Point p2, Point p3) {
  double a, b, c;

  // Lengths of sides
  a = sqrt(squared(p2.x - p1.x) + squared(p2.y - p1.y));
  b = sqrt(squared(p3.x - p2.x) + squared(p3.y - p2.y));
  c = sqrt(squared(p1.x - p3.x) + squared(p1.y - p3.y));

  // Eqn for area of triangle
  return sqrt((a + (b + c) ) * (c - (a - b) ) * (c + (a - b) ) * (a + (b - c) )) / 4;
}
*/



// Main render loop:
    /* Slow code but works
       for(j = 0; j < rows; j++) {
         for(y = ypoint[j]; y < ypoint[j + 1]; y++) {
           for(i = 0; i < cols; i++) {
             // Look up colour
             colour = leg_colours->lookup(data[j * cols + i]);
             for(x = xpoint[i]; x < xpoint[i + 1]; x++) {
               img[y * img_width + x] = colour * map[y * img_width + x];
	     }
	   }
	 }
       }
    */

  /*
  topright.y = p1.y;
  bottomleft.y = p3.y;
  // Figure out bounding box
  if(p1.x > p2.x) {
    if(p1.x > p3.x) {
      bottomright.x = p1.x;
      // Figure out which way p2 and p3 are ordered
      if(p2.x < p3.x) {
	topleft.x = p3.x;
      }
    } else { // if p1.x < p3.x
      // We know p1 is in the middle and p3.x is on bottom
      topleft.x = p2.x;
      bottomright.x = p3.x;
    }
  } else { // if p1.x < p2.x
    if(p1.x > p3.x) {
      // We know p1 is in the middle and p2.x is on bottom
      topleft.x = p3.x;
      bottomright.x = p2.x;
    } else { // if p1.x < p3.x
      topleft.x = p1.x;
      if(p2.x < p3.x) {
	bottomright.x = p3.x;
      }
    }
  }    
  */
  
  /*
  Line<Point<float> > l1, l2;
  if(p1.y - p2.y == 0) {
    l1 = Line<Point<float> >(p2, p3, true);
  } else {
    l1 = Line<Point<float> >(p1, p2, true);
  }
  l2 = Line<Point<float> >(p1, p3, true);

  list<Line<Point<float> > > polygon;
  */

/*
bool compute_outline_change(list<Point<float> >& outline,
			       const Line<Point<float> >& polyline, 
			       const Line<Point<float> >& boxline, 
			       Point<float> topright, 
			       Point<float> bottomleft, 
			       bool& inside) {
  
  Point<float>* intersect = ls_intersect(boxline, polyline);

  if(intersect) {
    if(within_bbox(polyline.from, topright, bottomleft)) {
      // Exiting, since the 'from' point is within the box
      inside = false;
      outline.push_back(Point<float>(*intersect));
    } else if(within_bbox(polyline.to, topright, bottomleft)) {
      // Entering, since the 'from' point is within the box
      inside = true;
      outline.push_back(Point<float>(*intersect));
      outline.push_back(Point<float>(polyline.to));
    } else {
      // Passing through
      // Need to figure out direction
      
      // This is a left or right side
      //if((sign(intersect->y - polyline.from.y) * sign(intersect->x - polyline.from.x) * sign(boxline.to.x - boxline.from.x) * sign(boxline.to.y - boxline.from.y)) > 0) {
      //if((sign(intersect->y - polyline.from.y) * sign(intersect->x - polyline.from.x) * sign(boxline.slope)) > 0) {
      //if(sign(polyline.slope) * sign(boxline.slope) * sign(normal_slope(boxline)) > 0) {
      
      // This doesn't work -- further analysis needed
      // 
      //if(sign(polyline.slope) > sign(normal_slope(boxline))) {
      //if(sign(polyline.slope) > (sign(intersect.x - boxline.from.x) * sign(intersect.y - boxline.from.y))) {

      //if(sign(polyline.slope) <= sign(boxline.slope)) {
      if(sign(polyline.slope) <= sign(boxline.slope)) {
	// Going to outside
	inside = false;
      } else {
	// Going to inside
	inside = true;
      }
      outline.push_back(Point<float>(*intersect));
    }
  } else {
    return false;
  }
  return true;
}  

*/

/*
      // If no intersects...
      // Change a few things
      if(outline.size() == 0) {
	// If one point is within and no intersects, whole thing is within
	if((llon < p1.x) && (rlon > p1.x) && (blat < p1.y) && (tlat > p1.y)) {
	  for(pl_iter = polylines.begin(); pl_iter != polylines.end(); pl_iter++) {
	    outline.push_back((*pl_iter).from);
	  }
	} else {
	  // Either it's completely out or completely in
	  // If out... figure out if it's within the triange
	  // How?
	}
      } else {
	// We have a triangle

	if(outline.size() < 3) {
	  cout << topright << " " << bottomleft << ": ";
	  cout << "Error: Can't have a polygon made up of less than 3 points" << endl;
	  grid[i * cols + j] = 0;
	} else {
	  // We have a valid poly - calculate area of polys making this one up
	  // Triangle fan approach
	  ol_iter = outline.begin();
	  Point<float> pivot, prev, cur;
	  pivot = *ol_iter;

	  pl_iter++;
	  prev = *ol_iter;

	  pl_iter++;
	  cur = *ol_iter;

	  // Zero out the grid box
	  grid[i * cols + j] = 0;
	  while(ol_iter != outline.end()) {
	    // Add area of this triangle to the grid box
	    grid[i * cols + j] += triangle_area(pivot, prev, cur);
	    prev = cur;
	    ol_iter++;
	    cur = *ol_iter;
	  }
	  cout << topright << "-" << bottomleft << ": " << grid[i * cols + j] << endl;
	  
	}
      }
      
*/

/*
float normal_slope(Line<Point<float> > line) {
  float slope;
  if(line.to.y == line.from.y) {
    slope = INFINITY * -(line.to.x - line.from.x);
  } else {
    slope = -(line.to.x - line.from.x) / (line.to.y - line.from.y);
  }
  return slope;
}

int sign(float number) {
  return (number >= -0) ? 1 : -1;
}

*/

/*
// Assumes incoming shapes are POLYGONS -- NOT overlapping
void draw_polygon(int rows, int cols, float* grid, float* lats, float* lons, Point<float>** points, int numpoints) {
  int i, j;

  if(numpoints < 3) {
    // Case 0: Not a polygon
    cout << __FILE__ << ":" << __LINE__ << ": Error: Can't have a poly with less than 3 sides" << endl;
  } else if(numpoints == 3) {
    // Case 1: Triangle (ready to go)
    draw_triangle(rows, cols, grid, lats, lons, *(points[0]), *(points[1]), *(points[2]));
  } else {
    Line<Point<float> > lines[numpoints];
  
    // Initialize lines array
    for(i = 0; i < numpoints; i++) {
      lines[i].from = *points[i];
      lines[i].to = *points[(i + 1) % numpoints];
    }

    int piv, nxt, prv;
    int numtris = 0;
    // Case 2: Polygon
    piv = 0;
    next_point_set(piv, nxt, prv, lines, points, numpoints);
    // Is this a flawed test?
    j = 0;
    while(numtris < (numpoints - 2)) {
      j++;
      bool points_within_triangle = false;
      if(nxt == -1) {
	next_point_set(piv, prv, nxt, lines, points, numpoints);
      }

      cout << piv << " " << prv << " " << nxt << endl;
      for(i = 0; i < numpoints; i++) {
	if(i == piv || i == prv || i == nxt) {
	  continue;
	} else {
	  bool ss1, ss2, ss3;
	  ss1 = sameside(*points[i], *points[piv], Line<Point<float> >(*points[prv], *points[nxt]));
	  ss2 = sameside(*points[i], *points[nxt], Line<Point<float> >(*points[piv], *points[prv]));
	  ss3 = sameside(*points[i], *points[prv], Line<Point<float> >(*points[prv], *points[nxt]));
	  // The IDEA of this test (not the implementation here) is to simply check if each point except the points making up the bounding triangle are within the triangle. Simple. Or so I thought.
	  if(ss1 && ss2 && ss3) {
	    points_within_triangle = true;
	    break;
	  }
	  // Replace this with angle summing? Sum angles of angles attached to points; if less than 180 there's a point inside
	}
      }
      if(points_within_triangle) {
	// Poly has something sticking into it which isn't part of it
	next_point_set(piv, prv, nxt, lines, points, numpoints);
      } else if(has_internal_concave(piv, prv, nxt, points, numpoints)) {
	// Points are not correctly ordered -- triangle is not part of poly
	next_point_set(piv, prv, nxt, lines, points, numpoints);
      } else if(((nxt + 1) % numpoints) != piv && line_crosses(piv, nxt, points, numpoints)) {
	// Triangle's lines cross other polys
	next_point_set(piv, prv, nxt, lines, points, numpoints);
      } else if(((prv + (numpoints - 1)) % numpoints) != piv && line_crosses(piv, prv, points, numpoints)) {
	// Triangle's lines cross other polys
	next_point_set(piv, prv, nxt, lines, points, numpoints);
      } else if(((prv + (numpoints - 1)) % numpoints) != nxt && line_crosses(prv, nxt, points, numpoints)) {
	// Triangle's lines cross other polys
	next_point_set(piv, prv, nxt, lines, points, numpoints);
      } else {
	// Triangle is all proper
	// draw_triangle(rows, cols, grid, lats, lots, piv, prv, nxt);
	cout << "Would draw triangle " << *(points[piv]) << ", " << *(points[prv]) << ", " << *(points[nxt]) << endl;
	set_lines_used(piv, prv, nxt, lines, numpoints);
	prv = nxt;
	nxt = next_point_from(nxt, lines, points, numpoints);
	numtris++;
      }
    }
  }
}
*/
  
/*
void set_lines_used(int p1, int p2, Line<Point<float> >* lines, int numlines) {
  if(p1 == (p2 + 1) % numlines) {
    lines[p2].used = true;
  } else if(p2 == (p1 + 1) % numlines) {
    lines[p1].used = true;
  }
}
void set_lines_used(int p1, int p2, int p3, Line<Point<float> >* lines, int numlines) {
  set_lines_used(p1, p2, lines, numlines);
  set_lines_used(p2, p3, lines, numlines);
  set_lines_used(p3, p1, lines, numlines);
}

// Returns true if angle at point p is greater than 180 degrees or 2pi radians
bool is_concave(int p, Point<float> **points, int numpoints) {
  int l, r;
  double costheta, sintheta;
  double llength, rlength;
  l = (p + 1) % numpoints;
  r = (p + (numpoints - 1)) % numpoints;
  llength = sqrt((*(points[l]) - *(points[p])) * (*(points[l]) - *(points[p])));
  rlength = sqrt((*(points[r]) - *(points[p])) * (*(points[r]) - *(points[p])));
  costheta = (*(points[l]) - *(points[p])) * (*(points[r]) - *(points[p]));
  costheta /= llength * rlength;
  sintheta = sqrt(1 - costheta * costheta);
  if(points[p]->y - sintheta * llength >= points[l]->y) {
    return false;
  } else {
    return true;
  }
}

// Gets the next point set. Caveats from above apply.
void next_point_set(int& piv, int& nxt, int& prv, Line<Point<float> >* lines, Point<float> **points, int numpoints) {
  piv = next_point_from(piv, lines, points, numpoints);
  prv = next_point_from(piv, lines, points, numpoints);
  nxt = next_point_from(prv, lines, points, numpoints);
}

bool has_internal_concave(int p1, int p2, int p3, Point<float> **points, int numpoints) {
  if(p3 < p1) {
    return has_internal_concave(p3, p2, p1, points, numpoints);
  }
  if(p3 < p2) {
    return has_internal_concave(p1, p3, p2, points, numpoints);
  }
  if(p2 < p1) {
    return has_internal_concave(p2, p1, p3, points, numpoints);
  }
  if(p2 == (p1 + 1) % numpoints && p3 == (p2 + 1) % numpoints) {
    bool cwpiv = clockwise(p2, (p2 + (numpoints - 1)) % numpoints, (p2 + 1) % numpoints, points, numpoints);
    if(cwpiv) {
      return true;
    } else {
      return false;
    }
  }
  return false;
}

int next_point_from(int p, Line<Point<float> >* lines, Point<float> **points, int numpoints) {
  int np;
  if(p == -1) {
    return -1;
  }
  for(np = (p + 1) % numpoints; np != p; np = (np + 1) % numpoints) {
    if(lines[(np + (numpoints - 1)) % numpoints].used == false || lines[np].used == false) {
      return np;
    }
  }
  return -1;
}


//Clockwise tests
  //return ((*(points[nxt]) - *(points[piv])).cross(*(points[prv]) - *(points[piv])) > 0);
  //return ((*(points[prv]) - *(points[piv])).cross(*(points[nxt]) - *(points[piv])) > 0);

*/
