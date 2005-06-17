#include "support.h"

enum CONN_STATES { UNDEF, USED, LINE, TEMP };

bool sameside(const Point<float>& p1, const Point<float>& p2, const Line<Point<float> >& l) {
  float crossproduct = (l.to - l.from).cross(p1 - l.from) * (l.to - l.from).cross(p2 - l.from);
  if(crossproduct >= 0) {
    return true;
  } else {
    return false;
  }
}

int sign(float number) {
	  return (number >= -0) ? 1 : -1;
}

// Internal
// Swaps things in the swap array if they're not ordered
void swap_if_unordered(int* swap_array, 
		       Point<float>* points[], 
		       Point<float> from, 
		       int to1, 
		       int to2)
{
  if(points[swap_array[to1]] && points[swap_array[to2]]) {
    if(from.dist_from(*(points[swap_array[to1]])) > from.dist_from(*(points[swap_array[to2]]))) {
      swap(swap_array[to1], swap_array[to2]);
    }
  } else {
    if(!points[swap_array[to1]]) {
      // Try to move NULLs to end
      swap(swap_array[to1], swap_array[to2]);
    }
  }
}

// Returns true if above line, false otherwise
bool right_of_line(const Point<float>& p, const Line<Point<float> >& l) {
  if(l.slope == 0) {
    return (((l.to.x - l.from.x) * (l.from.y - p.y)) > 0);
  } else if(l.slope == INFINITY || l.slope == -INFINITY) {
    return (((l.from.x - p.x) * (l.to.y - l.from.y)) < 0);
  } else {
    return (((((p.x - l.from.x) * l.slope) + l.from.y - p.y) * (l.to.y - l.from.y) * (l.to.x - l.from.x)) > 0);
  }
}

// Assumes rectangle
bool within_bbox(Point<float> p, Point<float> topright, Point<float> bottomleft) {
  return ((bottomleft.x < p.x) && (topright.x > p.x) && (bottomleft.y < p.y) && (topright.y > p.y));
}

// Returns list of points composing bounding box in clockwise order
Line<Point<float> >* get_box_outline(Point<float> topright, Point<float> bottomleft) {
  Line<Point<float> >* lines = new Line<Point<float> >[QUAD_SIDES];
  
  lines[0] = Line<Point<float> >(Point<float>(bottomleft.x, topright.y), Point<float>(topright.x, topright.y), true);
  lines[1] = Line<Point<float> >(Point<float>(topright.x, topright.y), Point<float>(topright.x, bottomleft.y), true);
  lines[2] = Line<Point<float> >(Point<float>(topright.x, bottomleft.y), Point<float>(bottomleft.x, bottomleft.y), true);
  lines[3] = Line<Point<float> >(Point<float>(bottomleft.x, bottomleft.y), Point<float>(bottomleft.x, topright.y), true);
  
  return lines;
}

Point<float>* ls_intersect(const Line< Point<float> >& l1, const Line< Point<float> >& l2, bool ignore_first_range) {
  return ls_intersect(l1.from, l1.to, l2.from, l2.to, ignore_first_range);
}

Point<float>* ls_intersect(const Point<float>& p00, const Point<float>& p01, const Point<float>& p10, const Point<float>& p11, bool ignore_first_range) {
  float mA, mB, x, y, s, t;
  int vlA = 0, vlB = 0;
  Point<float>* point = 0;

  // Check for vertical lines

  if((p00.x - p01.x) == 0) {
    // A is vertical line
    vlA = 1;
    mA = INFINITY;
  } else {
    mA = (p00.y - p01.y) / (p00.x - p01.x);
  }

  if((p10.x - p11.x) == 0) {
    // B is vertical line
    vlB = 1;
    mB = INFINITY;
  } else {
    mB = (p10.y - p11.y) / (p10.x - p11.x);
  }

  //mA = clip(mA);
  //mB = clip(mB);

  // If parallel or antiparallel...
  if(mA == mB) {
    return 0;
  }

  if(vlA) {
    //cout << "(" << p00.x << " - " << p11.x << ") = " << (p00.x - p11.x) << endl;

    // Line A (p00 - p01) is a vertical line
    x = p00.x;
    y = ((p00.x - p11.x) * mB + p11.y);
    s = ((y - p00.y) / (p01.y - p00.y));
    t = ((x - p10.x) / (p11.x - p10.x));
  } else if(vlB) {
    // Line B (p10 - p11) is a vertical line
    x = p10.x;
    y = ((p10.x - p01.x) * mA + p01.y);
    s = ((x - p00.x) / (p01.x - p00.x));
    t = ((y - p10.y) / (p11.y - p10.y));
  } else {
    // Neither line is a vertical line
    x = ( - mB * p10.x + p10.y + mA * p00.x - p00.y ) / ( mA - mB );
    y = (mA * ( x - p00.x ) + p00.y);
    s = ((x - p00.x) / (p01.x - p00.x));
    t = ((x - p10.x) / (p11.x - p10.x));
  }

  // Check range
  if((ignore_first_range || (s >= 0 && s <= 1)) && t >= 0 && t <= 1) {
    point = new Point<float>(x, y);
  }

  // Return whatever we got, or 0 (NULL) if no intersect
  return point;
}

// Assumes points come in in clockwise order - THIS IS REQUIRED
void draw_triangle(int rows, int cols, float* grid, float* lats, float* lons, Point<float> p1, Point<float> p2, Point<float> p3) {
  int i, j, k, l;
  float maxlat = lats[0], minlat = lats[rows];
  float maxlon = lons[cols], minlon = lons[0];
  Point<float> topright, bottomleft;

  //cout << p1 << p2 << p3 << endl;

  Line<Point<float> > line[] = { 
    Line<Point<float> >(p1, p2, true), 
    Line<Point<float> >(p2, p3, true), 
    Line<Point<float> >(p3, p1, true)
  };

  int numlines = sizeof(line) / sizeof(Line<Point<float> >);
  Point<float>* ipoints[numlines];
  Line<Point<float> >* box;

  for(i = 0; i < rows; i++) {
    // Reverse the lats -- arg. I want them to go min -> max
    float tlat = lats[i];
    float blat = lats[i + 1];

    //cout.setf(ios::fixed);
    //cout << tlat << "-" << blat << ": ";


    for(j = 0; j < cols; j++) {
      float llon = lons[j];
      float rlon = lons[j + 1];
      topright = Point<float>(rlon, tlat);
      bottomleft = Point<float>(llon, blat);

      list<Point<float> > outline;
      list<Point<float> >::iterator ol_iter;

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
	
	for(l = 0; ipoints[sar[l]]; l++) {
	  //for(l = 0; l < numlines; l++) {
	  if(ipoints[sar[l]]) {
	    // Flip 'inside' bit
	    inside = !inside;
	    if(within_bbox(line[sar[l]].from, topright, bottomleft)) {
	      // Going out
	      outline.push_back(line[sar[l]].from);
	    }
	    outline.push_back(*(ipoints[sar[l]]));
	    delete ipoints[sar[l]];
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
	  Point<float> pivot, prev, cur;
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
	  cout << __FILE__ << ":" << __LINE__ << ": Error: Can't have a poly with less than 3 sides" << endl;
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

bool clockwise(int piv, int nxt, int prv, Point<float> **points, int numpoints) {
  return ((*(points[nxt]) - *(points[piv])).cross(*(points[piv]) - *(points[prv])) > 0);
}

bool line_crosses(int from, int to, Point<float> **points, int numlines) {
  int nxt;
  Point<float>* intersect;
  for(int i = 0; i < numlines; i++) {
    nxt = (i + 1) % numlines;
    // Make sure that the lines do not include endpoints
    if(from == i || from == nxt) {
      // Do no checks
    } else {
      // Do not check endpoints
      if((from != i) && (from != nxt) && (to != i) && (to != nxt)) {
	// Check if the line segments intersect
	intersect = ls_intersect(*(points[from]), *(points[to]), *(points[i]), *(points[nxt]));
	if(intersect) {
	  // If so, clean up and say they do
	  delete intersect;
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

void remove_intercepts(int p1, int p2, Point<float>** points, int connects[], int numpoints) {
  // If temp line...
  if(connects[(p1 * numpoints) + p2] == TEMP) {
    for(int i = 0; i < numpoints; i++) {
      if(i == p1 || i == p2) {
	continue;
      }
      for(int j = 0; j < numpoints; j++) {
	if(j == p1 || j == p2) {
	  continue;
	}
	
	if(connects[(i * numpoints) + j] == TEMP) {
	  // Do intersect checking
	  Point<float> *intersect;
	  intersect = ls_intersect(*points[p1], *points[p2], *points[i], *points[j]);
	  if(intersect) {
	    connects[(i * numpoints) + j] = USED;
	    delete intersect;
	  }
	}
      } // end for
    } // end for
  }
}

// Decrements this intersect and its equivalent at other end of line
void dec(int p1, int p2, int connects[], int numpoints) {
  if(p1 == p2) {
    // CANT HAPPEN
  } else {
    connects[(p1 * numpoints) + p2]--;
    connects[(p2 * numpoints) + p1]--;
  }
}

// Adaptation of algorithm found at:
// http://www.ecse.rpi.edu/Homepages/wrf/research/geom/pnpoly.html
// Copyright (c) 1970-2003, Wm. Randolph Franklin
bool is_inside(Point<float> p, Point<float> **points, int numpoints) {
  bool inside = false;
  int i, j;
  for(i = 0, j = numpoints - 1; i < numpoints; j = i++) {
    if ((((points[i]->y <= p.y) && (p.y < points[j]->y)) ||
	 ((points[j]->y <= p.y) && (p.y < points[i]->y))) &&
	(p.x < (points[j]->x - points[i]->x) * (p.y - points[i]->y) / (points[j]->y - points[i]->y) + points[i]->x)) {
      inside = !inside;
    }
  }
  return inside;
}

// Assumes incoming shapes are POLYGONS -- NOT overlapping
void draw_polygon(int rows, int cols, float* grid, float* lats, float* lons, Point<float>** points, int numpoints) {
  int connects[numpoints * numpoints];
  int* conn;
  int i, j;

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
      if(!line_crosses(i, j, points, numpoints)) {
	// If this is inside...
	if(is_inside((*points[i] + *points[j]) / 2, points, numpoints)) {
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
      remove_intercepts(pts[0], pts[1], points, connects, numpoints);
      remove_intercepts(pts[1], pts[2], points, connects, numpoints);
      remove_intercepts(pts[2], pts[0], points, connects, numpoints);

      // Mark lines taken
      dec(pts[0], pts[1], connects, numpoints);
      dec(pts[1], pts[2], connects, numpoints);
      dec(pts[2], pts[0], connects, numpoints);
      
      //cerr << pts[0] << " " << pts[1] << " " << pts[2] << endl;
      //cerr << "Triangle at " << *points[pts[0]] << "-" << *points[pts[1]] << "-" << *points[pts[2]] << endl;

      draw_triangle(rows, cols, grid, lats, lons, *points[pts[0]], *points[pts[1]], *points[pts[2]]);
    } else {
      // Done with this point
      pts[0]++;
    }
  } // end while

} // end draw_polygon

// Replace with spherical tri area
float triangle_area(Point<float> p1, Point<float> p2, Point<float> p3) {
  double a, b, c;

  // Lengths of sides
  a = sqrt(squared(p2.x - p1.x) + squared(p2.y - p1.y));
  b = sqrt(squared(p3.x - p2.x) + squared(p3.y - p2.y));
  c = sqrt(squared(p1.x - p3.x) + squared(p1.y - p3.y));

  // Eqn for area of triangle
  return sqrt((a + (b + c) ) * (c - (a - b) ) * (c + (a - b) ) * (a + (b - c) )) / 4;
}
