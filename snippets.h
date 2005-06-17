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
