#include "support.h"
#include "line.h"

enum CONN_STATES { UNDEF, USED, LINE, TEMP };

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
		       Point* points[], 
		       Point from, 
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
bool within_bbox(Point p, Point topright, Point bottomleft) {
  return ((bottomleft.x < p.x) && (topright.x > p.x) && (bottomleft.y < p.y) && (topright.y > p.y));
}

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
  Point* ipoints[numlines];
  Line* box;

  for(i = 0; i < rows; i++) {
    // Reverse the lats -- arg. I want them to go min -> max
    double tlat = lats[i];
    double blat = lats[i + 1];

    //cout.setf(ios::fixed);
    //cout << tlat << "-" << blat << ": ";


    for(j = 0; j < cols; j++) {
      double llon = lons[j];
      double rlon = lons[j + 1];
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

bool clockwise(int piv, int nxt, int prv, Point **points, int numpoints) {
  return ((*(points[nxt]) - *(points[piv])).cross(*(points[piv]) - *(points[prv])) > 0);
}

bool line_crosses(int from, int to, Point **points, int numlines) {
  int nxt;
  Point* intersect;
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

void remove_intercepts(int p1, int p2, Point** points, int connects[], int numpoints) {
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
	  Point *intersect;
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
bool is_inside(Point p, Point **points, int numpoints) {
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
void draw_polygon(int rows, int cols, double* grid, const double* lats, const double* lons, Point** points, int numpoints) {
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
double triangle_area(Point p1, Point p2, Point p3) {
  double a, b, c;

  // Lengths of sides
  a = sqrt(squared(p2.x - p1.x) + squared(p2.y - p1.y));
  b = sqrt(squared(p3.x - p2.x) + squared(p3.y - p2.y));
  c = sqrt(squared(p1.x - p3.x) + squared(p1.y - p3.y));

  // Eqn for area of triangle
  return sqrt((a + (b + c) ) * (c - (a - b) ) * (c + (a - b) ) * (a + (b - c) )) / 4;
}
