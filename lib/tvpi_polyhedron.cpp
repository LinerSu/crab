/* Derived from the TVPI library (original file: polyhedron.cpp).
 * Original author: Axel Simon <A.Simon@kent.ac.uk>
 * Original source: https://github.com/axelsimon/tvpi
 *
 * This file is distributed under the GNU General Public License version 2.
 * See include/crab/domains/tvpi_impl/COPYING for the full license text.
 *
 * Modifications for Crab (2026): Yusen Su <yusen.su@uwaterloo.ca>
 *   - Rewrote include paths to use bundled tvpi_impl/ headers
 *   - Fixed NDEBUG-mode bugs: polyhedron macro parens, dir variable scope
 *   - Fixed axis-aligned Inequality construction (mpq_class -> mpz_class)
 *   - Added #define isIntegral isZ (originally from affine.hh)
 */
// polyhedron.cpp
// A planar polyhedron represented by a finite number of inequalities.
//

#include <crab/domains/tvpi_impl/polyhedron.hh>
#include <crab/domains/tvpi_impl/planar.hh>
#include <crab/domains/tvpi_impl/interval.hh>
#include <crab/domains/tvpi_impl/memory.hh>
#include <algorithm>

#undef DEBUG_HULL
#undef DEBUG_RESULTANTS
#undef DEBUG_QUASISYN
#undef DEBUG_REDUNDANT
#define DEBUG_INTEGRALHULL
#undef DEBUG_ENTAILMENT
#undef DEBUG_ENFORCE
#undef DEBUG_WIDEN
#undef DEBUG_LINOPT

#ifdef DEBUG_MEMORY
#include <iostream>
#include <string>

using std::cerr;
using std::endl;

#endif

#ifdef NDEBUG

#undef DEBUG_HULL
#undef DEBUG_RESULTANTS
#undef DEBUG_QUASISYN
#undef DEBUG_REDUNDANT
#undef DEBUG_INTEGRALHULL
#undef DEBUG_ENTAILMENT
#undef DEBUG_ENFORCE

#else // NDEBUG

#include <iostream>

using std::cerr;
using std::endl;

#endif

using std::vector;
using std::sort;

// TODO: 

// - most of these functions take the size of an array and an array of
// inequalities. Although syntactically a vector looks exactly the
// same, these types are not exchangable. I use the hack &blah[0] to
// pass a vector as an array which obviously assumes that a vector is
// implemented as an array. Duplicating the code for this with
// templates just seems to be nonsense.

// - in case that adding only redundant inequalities to a shared
// polyhedron, the polyhedron and all its inequalities is copied
// although it is no different from the previous one. Keeping an
// additional "twin" pointer could prevent copying the polyhedron the
// next time inequalities are added.

// - in insertInequalities: Inequalities that are two quadrants apart
// but intersect will tighten both bounds by the intersection
// point. However, if one of the bounds is already tighter than the
// intersection point, the other inequality will induce a tighter
// bound on the other axis, too. There is a test to be
// avoided. Somewhere.

// necessity to use the C sorting functions
extern "C" typedef int (*compareFunc) (const void*, const void*);

namespace Tvpi {

// Define a function for checking that a list of inequalities is in order.
#ifdef NDEBUG

#define inequalitiesSorted(size, facet) true

#else // NDEBUG

#define inequalitiesSorted(size, facet) \
  (size>0 ? inequalitiesSorted_(size, facet) : true)

bool inequalitiesSorted_(size_t size, Inequality* facet[]) {
  assert(size>0);
  bool res=true;
  for(size_t i=0; i<size-1; i++) {
    if (*facet[i]>=*facet[i+1]) {
      cerr << i << "th inequality " << *facet[i] << " >= "
	   << i+1 << "th inequality " << *facet[i+1] << endl;
      res=false;
    };
  };
  return res;
}

#endif // NDEBUG

// Partitioning class

// Create partitioning information for a set of inequalities.
Partitioning::Partitioning(size_t facetSize, Inequality** facet) {

  // The first and the last index of the partitioning array are
  // obvious and only there as a convenience.
  dirStart[east]=0;
  dirStart[total]=facetSize;

  // Set the entry west by using entries east and total.
  binarySearch(facet,east,total);

  // Set the entry north by using entries east and west.
  binarySearch(facet,east,west);

  // Set the entry south by using entries west and total.
  binarySearch(facet,west,total);
};

// Create partitioning information for a vector of inequalities.
Partitioning::Partitioning(const vector<Inequality* > facet) {
  if (facet.empty()) for(int d=east; d<numDirs; d++) dirStart[d]=0;
  else {
    // The first and the last index of the partitioning array are
    // obvious and only there as a convenience.
    dirStart[east]=0;
    dirStart[total]=facet.size();
    
    // Set the entry west by using entries east and total.
    binarySearch(&facet.front(),east,total);
    
    // Set the entry north by using entries east and west.
    binarySearch(&facet.front(),east,west);
    
    // Set the entry south by using entries west and total.
    binarySearch(&facet.front(),west,total);
  };
}

// Determine the indices where the inequalities change the direction.
void Partitioning::binarySearch(Inequality* const facet[],
				Direction dirLower,
				Direction dirUpper) {
  // Do a binary search to determine where the direction of the
  // inequalities change.
  size_t lower=dirStart[dirLower];
  size_t upper=dirStart[dirUpper];
  Direction goalDir=Direction((dirLower+dirUpper)/2);
  while (lower<upper) {
    size_t middle=(lower+upper)/2;
    if (facet[middle]->calcDirection()<goalDir) {
      lower=middle+1;
    } else {
      upper=middle;
    };
  };
  // The search found the first entry which belongs to the next
  // partition.
  dirStart[goalDir]=lower;
}

void Partitioning::getNumOfInequalitiesWithSign(bool secondCoeff,
						size_t* plusCoeff,
						size_t* minusCoeff) {
  // Setup the indices for accessing inequalities of specific
  // signs. The last... Variables are the first directions which are
  // not part of the segment anymore.
  Direction firstPlus= (secondCoeff ? east  : south);
  Direction lastPlus=  (secondCoeff ? west  : north);
  Direction firstMinus=(secondCoeff ? west  : north);
  Direction lastMinus= (secondCoeff ? total : south);
  if (plusCoeff) {
    if (firstPlus>lastPlus) {
      *plusCoeff=dirStart[total]-dirStart[firstPlus]+dirStart[lastPlus]; 
      //-dirStart[east]; (but that term is always 0)
    } else {
      *plusCoeff=dirStart[lastPlus]-dirStart[firstPlus];
    }
  };
  if (minusCoeff) *minusCoeff=dirStart[lastMinus]-dirStart[firstMinus];
}

#define incrF1(num) (num+1==f1Size ? 0 : num+1)
#define incrF2(num) (num+1==f2Size ? 0 : num+1)

// Calculate the resultants of two sets of inequalities.
template<bool isZ>
void calcResultantsOfSets(const Inequality* const facet1[],
			  const Partitioning& f1Part,
			  bool f1CommonCoeff,
			  const Interval<isZ>& f1Bound,
			  const Inequality* const facet2[],
			  const Partitioning& f2Part,
			  bool f2CommonCoeff,
			  const Interval<isZ>& f2Bound,
			  vector<Inequality*>& ineqOut) {
  size_t f1Size = f1Part[total];
  if (f1Size==0) return;
  size_t f2Size = f2Part[total];
  if (f2Size==0) return;

  // Combine the positive coefficients in facet1 with the negative
  // coefficients in facet2.
  size_t f1Start=f1Part[(f1CommonCoeff ? east : south)];
  size_t f1End=  (f1CommonCoeff ? f1Part[west] : f1Part[north]+f1Size);
  size_t f2Start=f2Part[(f2CommonCoeff ? west : north)];
  size_t f2End=  f2Part[(f2CommonCoeff ? total : south)];

#ifdef DEBUG_RESULTANTS
  cerr << "First set, " << f1Part << ", from #" << f1Start << " to #"
       << f1End << " with " << f2Part << ", #" << f2Start << " to #" 
       << f2End << endl;
#endif // DEBUG_RESULTANTS

  for(size_t f1Idx=f1Start; f1Idx<f1End; f1Idx++) {
    for(size_t f2Idx=f2Start; f2Idx<f2End; f2Idx++) {
#ifdef DEBUG_RESULTANTS
      cerr << "Combining (" << f1CommonCoeff << "+) "
	   << *facet1[f1Idx % f1Size] << " with ("
	   << f2CommonCoeff << "-) "
	   << *facet2[f2Idx % f2Size];
#endif // DEBUG_RESULTANTS
      ineqOut.push_back(new Inequality(isZ,
				       *facet1[f1Idx % f1Size],
				       f1CommonCoeff,
				       *facet2[f2Idx % f2Size], 
				       f2CommonCoeff));
#ifdef DEBUG_RESULTANTS
      cerr << " to " << *ineqOut.back() << endl;
#endif // DEBUG_RESULTANTS
    }
  };

  // Combine the negative coefficients in facet1 with the positive
  // coefficients in facet2.
  f1Start=f1Part[(f1CommonCoeff ? west : north)];
  f1End=  f1Part[(f1CommonCoeff ? total : south)];
  f2Start=f2Part[(f2CommonCoeff ? east : south)];
  f2End=  (f2CommonCoeff ? f2Part[west] : f2Part[north]+f2Size);

#ifdef DEBUG_RESULTANTS
  cerr << "Second set from #" << f1Start << " to #"
       << f1End << " with #" << f2Start << " to #" 
       << f2End << endl;
#endif // DEBUG_RESULTANTS

  for(size_t f1Idx=f1Start; f1Idx<f1End; f1Idx++) {
    for(size_t f2Idx=f2Start; f2Idx<f2End; f2Idx++) {
#ifdef DEBUG_RESULTANTS
      cerr << "Combining (" << f1CommonCoeff << "-) "
	   << *facet1[f1Idx % f1Size] << " with ("
	   << f2CommonCoeff << "+) "
	   << *facet2[f2Idx % f2Size];
#endif // DEBUG_RESULTANTS
      ineqOut.push_back(new Inequality(isZ,
				       *facet1[f1Idx % f1Size],
				       f1CommonCoeff,
				       *facet2[f2Idx % f2Size],
				       f2CommonCoeff));
#ifdef DEBUG_RESULTANTS
      cerr << " to " << *ineqOut.back() << endl;
#endif // DEBUG_RESULTANTS
    }
  }
};

// Sort a set of constraints and remove all quasi syntactic redundant
// inequalities.
void sortAndRemoveQuasiSyntacticRed(size_t ineqStart,
				    vector<Inequality*>& ineq) {
  size_t ineqSize=ineq.size()-ineqStart;
  if (ineqSize==0) return;

  // The constraints from calcResultantsOfSets are not ordered. It is
  // possible to create resultants as eight sequences of inequalities
  // where each sequence is increasing in angle. These could then be
  // more efficiently sorted by mergesort. Right now the result are
  // several sequences that are increasing or decreasing. Just use
  // quicksort.
  qsort(&ineq[ineqStart], ineqSize, sizeof(Inequality*),
	(compareFunc) &Inequality::compare);

  
  // Remove inequalities with the same angle (quasi-syntactic
  // redundancies, as Jaffar calls them).
  assert(ineq[ineqStart]->isTVPI());
  size_t ineqTop=ineq.size();
  // Track the last valid index.
  size_t lastIneq=ineqStart+1;
  for(size_t i=ineqStart+1; i<ineqTop; i++) {
    assert(ineq[i]->isTVPI());
    ineq[lastIneq]=ineq[i];
    if (*ineq[lastIneq-1]==*ineq[lastIneq]) {
#ifdef DEBUG_QUASISYN
      cerr << "one of " << *ineq[lastIneq-1] << " and "
	   << *ineq[lastIneq] << " is redundant: ";
#endif // DEBUG_QUASISYN
      // The inequalities have the same angle, delete the redundant one.
      if (ineq[lastIneq-1]->includes(*ineq[lastIneq])) {
#ifdef DEBUG_QUASISYN
	cerr << " the former" << endl;
#endif // DEBUG_QUASISYN
	delete ineq[lastIneq-1];
	ineq[lastIneq-1]=ineq[lastIneq];
      } else {
#ifdef DEBUG_QUASISYN
	cerr << " the latter" << endl;
#endif // DEBUG_QUASISYN
	delete ineq[lastIneq];
      };
    } else lastIneq++;
  };
  assert(lastIneq<=ineqTop);
  ineq.resize(lastIneq);

  assert(inequalitiesSorted(lastIneq-ineqStart, &ineq[ineqStart]));
};


// The invalid index, only used internally.
const size_t invalidIndex = (size_t) -1;

  // Create an empty polyhedron.
  template<bool isZ>
  PolyhedronImpl<isZ>::PolyhedronImpl() : refCount(1), twin(0) {};

  // Create a deep copy of the facet array.
  template<bool isZ>
  PolyhedronImpl<isZ>::PolyhedronImpl(const PolyhedronImpl& p) {
    facet.reserve(p.facet.size());
    for(size_t current=0; current<p.facet.size(); current++)
      facet.push_back(new Inequality(*p.facet[current]));
    refCount=1;
    part=p.part;
    twin = 0;
  };


// Calculate the convex hull of two non-empty polyhedra.
template<bool isZ>
void PolyhedronImpl<isZ>::convexHull(const PolyhedronImpl& p1,
				const Interval<isZ>& p1x,
				const Interval<isZ>& p1y,
				const PolyhedronImpl& p2,
				const Interval<isZ>& p2x,
				const Interval<isZ>& p2y,
				bool flipFirst) {
  // The number of points generated by extreme never exceeds the
  // number of inequalities. The set of inequalities in pX.facet does
  // not include the 4 inequalities that are generated on-the-fly for
  // the four the bounds.
  Point points[p1.facet.size()+p2.facet.size()+8];
  size_t lastPoint = 0;
  // The extreme function may create two colinear pairs of rays in
  // case the polyhedron constitues an equality or stripe.
  Ray rays[8];
  size_t lastRay = 0;

  // Setup the values of the class.
  refCount = 1;

  // This debug output is done in DenseTvpi::join()
#ifdef dont_DEBUG_HULL
  cerr << "convexHull: joining ";
  if (flipFirst) cerr << "flipped ";
  cerr << " first x=" << p1x << ", y=" << p1y << endl;
  output(cerr, p1, 0, 0);
  cerr << endl << "with second x=" << p2x << ", y= " << p2y << endl;
  output(cerr, p2, 0, 0);
  cerr << endl;
#endif // DEBUG_HULL

  // Insert all points and rays into the two arrays.
  p1.extreme(p1x, p1y, lastPoint, points, lastRay, rays);
#ifdef DEBUG_HULL
  if (flipFirst) {
    cerr << lastPoint << " points of p1 before flipping are: ";
    for(size_t current=0; current<lastPoint; current++)
      cerr << points[current] << ", ";
    cerr << endl << lastRay << " rays of p1 before flipping are: ";
    for(size_t current=0; current<lastRay; current++)
      cerr << rays[current] << ", ";
    cerr << endl;
    cerr << "second polyhedron x = " << p2x << ", y = " << p2y << endl;
  };
#endif // DEBUG_HULL
  if (flipFirst) {
    for (size_t p = 0; p<lastPoint; p++) points[p].swap();
    for (size_t r = 0; r<lastRay; r++) rays[r].swap();
  };
  p2.extreme(p2x, p2y, lastPoint, points, lastRay, rays);

#ifdef DEBUG_HULL
  cerr << "points are: ";
  for(size_t current=0; current<lastPoint; current++)
    cerr << points[current] << ", ";
  cerr << endl << "rays are: ";
  for(size_t current=0; current<lastRay; current++)
    cerr << rays[current] << ", ";
  cerr << endl;
#endif // DEBUG_HULL

  // There must be at least one point from each polyhedron.
  assert(lastPoint>=2);

  // If this fails its probably too late.
  assert(lastPoint<=p1.facet.size()+p2.facet.size()+8);

  // Calculate the minimum side length of a square that contains all
  // absolute coordinates of the points.
  mpz_class size;
  for(size_t current=0; current<lastPoint; current++)
    points[current].updateLength(size);

  // Create the square's side length.
  mpq_class s(size, mpz_class(1));

  // Ensure that all rays are longer than 3*s.
  for(size_t current=0; current<lastRay; current++)
    rays[current].ensureLength(size*3+1);


#ifdef DEBUG_HULL
  cerr << "square size is: " << s << endl << "new rays are: ";
  for(size_t current=0; current<lastRay; current++)
    cerr << rays[current] << ", ";
  cerr << endl;
#endif // DEBUG_HULL

  // Translate the points by the rays. Use a new scope to allocate the
  // array on the stack.
  {
    // the point set Q
    Point qpoints[lastPoint*(lastRay+1)];
    // the indices of the vertices p_k_0,...,p_k_n
    size_t vertices[lastPoint*(lastRay+1)];

    // copy the existing points
    size_t lastQPoint;
    for(lastQPoint=0; lastQPoint<lastPoint; lastQPoint++)
      qpoints[lastQPoint]=points[lastQPoint];
    
    // translate each point by each ray
    for(size_t curPoint=0; curPoint<lastPoint; curPoint++) {
      for(size_t curRay=0; curRay<lastRay; curRay++) {
	qpoints[lastQPoint++]=Point(rays[curRay], points[curPoint]);
      };
    };

#ifdef DEBUG_HULL
    cerr << "qpoints are: ";
    for(size_t current=0; current<lastQPoint; current++)
      cerr << qpoints[current] << ", ";
    cerr << endl;
#endif // DEBUG_HULL

    // As opposed to the description in the convex hull paper, we do
    // look for an extreme point which we then use as a pivot
    // point. The reason is that we want the resulting inequalities to
    // be sorted and not only merly consecutive.

    // The point with the highest x and lowest y coordinate
    // (lexicographically) is the pivot point. Connecting the pivot
    // point to the next in counter-clockwise order, gives us an
    // inequality closest in angle to |a| x <= c. Find the pivot point
    // and put it into qpoints[0].
    assert(lastQPoint>0);
    for (size_t i=1; i<lastQPoint; i++)
      if (qpoints[i]<qpoints[0]) qpoints[0].swapWith(qpoints[i]);

#ifdef DEBUG_HULL
    cerr << "the pivot point is " << qpoints[0] << endl;
    cerr << "qpoints before shifting: ";
    for(size_t current=0; current<lastQPoint; current++) {
      if (current>0) cerr << ", ";
      cerr << qpoints[current];
    };
    cerr << endl;
#endif // DEBUG_HULL

    // Run quicksort on the set of points. The standard qsort function
    // that comes with every C implementation does not take additional
    // data, so there is no way to pass the pivot point to the
    // comparison function. The crucial observation here is that the
    // comparison calculates the determinant of the coordinates of
    // (p1-pp | p2-pp) where pp is the pivot point. Hence we subtract
    // the pivot point from each point before we start the comparison
    // and use a compare function that just uses the origin as pivot
    // point.

    // Set pivotPoint to the pivot point and set qpoints[0] to <0,0>.
    Point pivotPoint = Point();
    pivotPoint.swapWith(qpoints[0]);
    for (size_t current=1; current<lastQPoint; current++) {
      qpoints[current]-=pivotPoint;
      // Check if this point coincides with the pivot point. If so,
      // remove it from the list.
      if (sgn(qpoints[current].getX())==0 &&
	  sgn(qpoints[current].getY())==0) {
	lastQPoint--;
	if (lastQPoint!=current)
	  qpoints[lastQPoint].swapWith(qpoints[current--]); // No increment.
      }
    };

#ifdef DEBUG_HULL
    cerr << "qpoints before sorting: ";
    for(size_t current=0; current<lastQPoint; current++) {
      if (current>0) cerr << ", ";
      cerr << qpoints[current];
    };
    cerr << endl;
#endif // DEBUG_HULL

    qsort(&qpoints[1], lastQPoint-1, sizeof(Point),
	  (compareFunc) &Point::ordering);

#ifdef DEBUG_HULL
    cerr << "qpoints after sorting: ";
    for(size_t current=0; current<lastQPoint; current++) {
      if (current>0) cerr << ", ";
      cerr << qpoints[current];
    };
    cerr << endl;
#endif // DEBUG_HULL

    // Check which of the first points are collinear with the pivot
    // point and sort them such that their x coordinate is
    // decreasing. Note that we need to sort the last colinear points
    // before the first colinear point since they all might be
    // colinear in which case this order prevails, leaving the correct
    // origin.
    {
      size_t lastColinear=2;
      while (lastColinear<lastQPoint)
	if (Point::calcDeterminant(&qpoints[lastColinear-1],
				   &qpoints[lastColinear])==0)
	  lastColinear++; else break;
#ifdef DEBUG_HULL
      cerr << "sorting starting colinear points: ";
      for(size_t current=0; current<lastColinear; current++) {
	if (current>0) cerr << ", ";
	cerr << qpoints[current];
      };
      cerr << endl;
#endif // DEBUG_HULL
      qsort(&qpoints[0], lastColinear, sizeof(Point),
	    (compareFunc) &Point::largerX);
    };

#ifdef DEBUG_HULL
    cerr << "qpoints before moving back: ";
    for(size_t current=0; current<lastQPoint; current++) {
      if (current>0) cerr << ", ";
      cerr << qpoints[current];
    };
    cerr << endl;
#endif // DEBUG_HULL

    // Move the polyhedron back to where it was.
    for(size_t current=0; current<lastQPoint; current++)
      qpoints[current]+=pivotPoint;

#ifdef DEBUG_HULL
    cerr << "qpoints after moving back: ";
    for(size_t current=0; current<lastQPoint; current++) {
      if (current>0) cerr << ", ";
      cerr << qpoints[current];
    };
    cerr << endl;
#endif // DEBUG_HULL

    // Graham scan.
    assert(lastQPoint>=1);

    // The first point is a definitive vertex. This simplifies the
    // scan as we will never wrap around to the end of the array.
    vertices[0]=0;
    size_t verticesSize=1;
    size_t last=0;
    size_t current=1 % lastQPoint;
    size_t next=2 % lastQPoint;
    while (current) {
      int res=qpoints[current].calcBend(qpoints[last],qpoints[next]);
      bool keep=res>0;
      if (res==0) {
	// If the three points are on a line, keep is false and we
	// normally discard the current point. However, we will keep
	// the point if the points are trivially colinear since there
	// are fewer than two points that take part in the comparison.
	if (last==next || last==current) keep=true; else
	  // Discard one of the colinear points. We discard the one
	  // that is closer to the last point.
	  if (qpoints[last].calcDistance(qpoints[current])>
	      qpoints[last].calcDistance(qpoints[next])) {
#ifdef DEBUG_HULL
	    cerr << "swapping current=" << current << " and next= " << next
		 << " point." << endl;
#endif // DEBUG_HULL
	    qpoints[current].swapWith(qpoints[next]);
	    assert(qpoints[current].calcBend(qpoints[last],qpoints[next])==0);
	  }
      };
#ifdef DEBUG_HULL
      cerr << (keep ? "keep" : "discard") << " point " << current << 
	", vertices: ";
	for(size_t index=0; index<verticesSize; index++)
	  cerr << vertices[index] << ", ";
	cerr << endl;
#endif // DEBUG_HULL
      if (keep) {
	// the points [last], [current], [next] are counter clockwise
	vertices[verticesSize++]=current;
	last=current;
	current=next;
	next=(next+1) % lastQPoint;
      } else {
	// the points [last], [current], [next] are clockwise or
	// collinear (where [current] is a convex combination of
	// [last] and [next]), remove [current]
	if (last) {
	  // if last is positive, vertices contains at least two
	  // potential vertices
	  current=last; // equals vertices[verticesSize-1]
	  verticesSize--;
	  last=vertices[verticesSize-1];
	} else {
	  // don't backtrack to the first point since we know it's a vertex
	  current=next;
	  next=(next+1) % lastQPoint;
	};
      };
    };

#ifdef DEBUG_HULL
    cerr << "Done, vertices: ";
    for(size_t index=0; index<verticesSize; index++) {
      if (index>0) cerr << ", ";
      cerr << "#" << vertices[index] << qpoints[vertices[index]];
    };
    cerr << endl;
#endif // DEBUG_HULL

    // Reconstitute the polyhedron. Again, a new scope to cater for
    // the temporary array of inequalities.
    {
      // Each vertex can be the starting point of at most one inequality.
      Inequality* newFacetsArray[verticesSize];
      Inequality** newFacets=&newFacetsArray[0];
      size_t lastNewFacet=0;
      
      for(size_t current=0; current<verticesSize; current++) {
	size_t next=(current+1) % verticesSize;
#ifdef DEBUG_HULL
	cerr << qpoints[vertices[current]] << "--"
	     << qpoints[vertices[next]] << ": ";
#endif // DEBUG_HULL
	// An Inequality between the current and the next vertex has
	// to be in the output if there is a point in the square that
	// saturates this inequality. This is trivially fullfilled if
	// one of the vertices that created the inequality is in the
	// square. In the one dimensional case, the two inequalities
	// between the two vertices must always be created, even if
	// none of the vertices are in the square. Dealing with this
	// case explitily allows us to make the loop only inspect all
	// points between two adjacent vertices in the normal two
	// dimensional case.
	bool add=qpoints[vertices[current]].inBox(s) || 
	  qpoints[vertices[next]].inBox(s) ||
	  verticesSize==2;
	// In case either the x or the y coordinate is the same in
	// this pair of points, the inequality is already
	// represented as a bound.
	bool isBound =
	  (qpoints[vertices[current]].getX()==
	   qpoints[vertices[next]].getX()) ||
	  (qpoints[vertices[current]].getY()==
	   qpoints[vertices[next]].getY());

#ifdef DEBUG_HULL
      if (add)
	if (isBound) cerr << "will create bound instead of inequality" << endl;
	else cerr << "added due to vertex in square" << endl;
#endif // DEBUG_HULL

	// We need to test if any points with indices between those of
	// the vertices that constitute the inequality lie in the
	// square. In this case we need to create an inequality
	// nevertheless (stripe). However, if the line does not
	// intersect with the square, then this is never going to
	// happen and we can safely skip the loop. We also skip the
	// loop if we already know that we have a valid inequality
	// (add or isBound is true).
	Inequality* e = NULL;
	if (!(add || isBound || outsideSquare(s, qpoints[vertices[current]],
					      qpoints[vertices[next]]))) {

	  size_t check=vertices[current];
	  if (++check==lastQPoint) check=0;

	  while (!add && check!=vertices[next]) {
	    // Examine the points that lie between the two
	    // vertices. If any of them lie within the square,
	    // create the inequality and check if the point is
	    // actually on the boundary of this inequality.
	    if (qpoints[check].inBox(s)) {
	      // Create the inequality lazily. If there is no point
	      // in the square, then we don't generate an
	      // inequality. If we generate the inequality and no
	      // point staturates it, the inequality will be deleted
	      // further down (!add holds).
	      if (!e) e=new Inequality(qpoints[vertices[current]],
				       qpoints[vertices[next]], isZ);
	      add=e->isSaturated(qpoints[check]);
#ifdef DEBUG_HULL
	      cerr << endl 
		   << "candidate " << *e << " is " << (add ? "" : "not ")
		   << "saturated by " << qpoints[check] << " ";
#endif // DEBUG_HULL
	    };
	    if (++check==lastQPoint) check=0;
	  };
	};
	if (add && !isBound) {
	  if (!e) e=new Inequality(qpoints[vertices[current]],
				   qpoints[vertices[next]], isZ);
	  newFacets[lastNewFacet++]=e;
#ifdef DEBUG_HULL
	  cerr << "added" << endl;
#endif // DEBUG_HULL
	} else {
	  if (e) delete e;
#ifdef DEBUG_HULL
	  cerr << "not added" << endl;
#endif // DEBUG_HULL
	};
      };

#ifdef DEBUG_HULL
      cerr << lastNewFacet << " inequalities: " << endl;
      for(size_t current=0; current<lastNewFacet; current++)
	cerr << *newFacets[current] << endl;
#endif // DEBUG_HULL

      // Copy the created inequalities to their final destination.
      facet.clear();
      facet.reserve(lastNewFacet);
      for (size_t current=0; current<lastNewFacet; current++)
	facet.push_back(newFacets[current]);
      
      // Update the partitioning information.
      part=Partitioning(facet);

      // We just created a new polyhedron, we can't have a twin.
      twin=0;

      assert(inequalitiesSorted(facet.size(), &facet.front()));
    };
  };
}

#define removeNullEntries(label) \
  if (tightenQuadrant(tempFacet, false, xBound, yBound, \
		      tempFacetStart, tempFacetSize)) { \
    while (tempFacet[tempFacetStart]) tempFacetStart++; \
    size_t upper=tempFacetStart; \
    do { \
      while (tempFacet[upper]==NULL) if (++upper==tempFacetSize) goto label; \
      tempFacet[tempFacetStart++]=tempFacet[upper++]; \
    } while (upper!=tempFacetSize); \
  label: \
    tempFacetSize=tempFacetStart; \
  };

// Calculate the points and rays of a polyhedron.
template<bool isZ>
void PolyhedronImpl<isZ>::extreme(const Interval<isZ>& xBound,
			     const Interval<isZ>& yBound,
			     size_t& lastPoint, Point points[],
			     size_t& lastRay, Ray rays[]) const {

  // Create an array of inequalities large enough to contain all
  // facets of this polyhedron and the four bounds as inequalities.
  size_t tempFacetSize=facet.size()+4;
  Inequality* tempFacet[tempFacetSize];

  // Fill the array with the bounds and the inequalities in each
  // partition. To our own bewilderment, we fill this array with
  // pointers to stack-allocated inequalities and heap-allocated
  // inequalities.
  Inequality boundEast, boundNorth, boundWest, boundSouth;
  tempFacetSize=0;
  size_t facetIdx=0; // =part[east]; -- same thing
  if (xBound.upperIsFinite()) {
    boundEast=Inequality(xBound.getUpper().get_den(), mpz_class(0), xBound.getUpper().get_num());
    tempFacet[tempFacetSize++]=&boundEast;
  };
  size_t tempFacetStart=tempFacetSize;
  while (facetIdx<part[north]) tempFacet[tempFacetSize++]=facet[facetIdx++];
  removeNullEntries(end_reached1);
  if (yBound.upperIsFinite()) {
    boundNorth=Inequality(mpz_class(0), yBound.getUpper().get_den(), yBound.getUpper().get_num());
    tempFacet[tempFacetSize++]=&boundNorth;
  };
  tempFacetStart=tempFacetSize;
  while (facetIdx<part[west]) tempFacet[tempFacetSize++]=facet[facetIdx++];
  removeNullEntries(end_reached2);
  if (xBound.lowerIsFinite()) {
    boundWest=Inequality(-xBound.getLower().get_den(), mpz_class(0), -xBound.getLower().get_num());
    tempFacet[tempFacetSize++]=&boundWest;
  };
  tempFacetStart=tempFacetSize;
  while (facetIdx<part[south]) tempFacet[tempFacetSize++]=facet[facetIdx++];
  removeNullEntries(end_reached3);
  if (yBound.lowerIsFinite()) {
    boundSouth=Inequality(mpz_class(0), -yBound.getLower().get_den(), -yBound.getLower().get_num());
    tempFacet[tempFacetSize++]=&boundSouth;
  };
  tempFacetStart=tempFacetSize;
  while (facetIdx<part[total]) tempFacet[tempFacetSize++]=facet[facetIdx++];
  removeNullEntries(end_reached4);

  assert(facetIdx==facet.size());
  assert(tempFacetSize<=facet.size()+4);
  assert(tempFacetSize>0);

  if (tempFacetSize==1) 
    // This polyhedron consists of a single halfspace, create a
    // ray that points into the feasible space.
    rays[lastRay++]=Ray(*tempFacet[0], Ray::dirToFeasible);

  size_t first=tempFacetSize-1;
  size_t second=0;
  while (second<tempFacetSize) {
    size_t third=second+1;
    if (third==tempFacetSize) third=0;

#ifdef DEBUG_HULL
    cerr << "first: (" << first << ") " << *tempFacet[first]
	 << " second: (" << second << ") " << *tempFacet[second]
	 << " third: (" << third << ")" << *tempFacet[third] << endl;
#endif // DEBUG_HULL

    // Check if the angle between the adjacent inequalities is greater
    // or equal to pi.
    bool dPre = (tempFacetSize==1) || 
      !tempFacet[first]->lessThanPi(*tempFacet[second]);
    bool dPost = (tempFacetSize==1) ||
      !tempFacet[second]->lessThanPi(*tempFacet[third]);

#ifdef DEBUG_HULL
    cerr << "dPre is " << (dPre ? "true" : "false")
	 << ", dPost is " << (dPost ? "true" : "false") << endl;
#endif // DEBUG_HULL

    if (dPre) rays[lastRay++]=Ray(*tempFacet[second], Ray::dirAgainstArrow);
    if (dPost) rays[lastRay++]=Ray(*tempFacet[second], Ray::dirWithArrow);
    else points[lastPoint++]=Point(*tempFacet[second], *tempFacet[third]);
    if (dPre && dPost) points[lastPoint++]=Point(*tempFacet[0]);
    first=second;
    second++;
  };
};

template<bool isZ>
PolyhedronImpl<isZ>::~PolyhedronImpl() {
  assert(refCount==0);
  assert((twin ? twin->refCount==0 : true));
  assert(facet.size()==part[total]);
  for(size_t current=0; current<facet.size(); current++)
    delete facet[current];
}


// Remove inequalities that are redundant with respect to the given
// bounds.
template<bool isZ>
bool PolyhedronImpl<isZ>::tightenQuadrant(Inequality* curFacet[],
				     bool deleteInequalities,
				     const Interval<isZ>& xBound,
				     const Interval<isZ>& yBound,
				     size_t start,
				     size_t stop) {
  // Ensure we can count upwards.
  assert(stop>=start);

  bool changed = false;

  // Iterate from start to stop. Firstly delete inequalities as long
  // as they include the box.
  size_t middle = start;

  // Stop here if there are no inequalities in this quadrant.
  if (middle==stop) return changed;
#ifdef DEBUG_ENFORCE
  cerr << "sweeping upwards from " << middle << " to " << stop 
       << ", testing box" << endl;
#endif // DEBUG_ENFORCE
  while (middle!=stop && curFacet[middle]->includes(xBound, yBound)) {
#ifdef DEBUG_ENFORCE
    cerr << "deleting (" << middle << ") " << *curFacet[middle] << endl;
#endif // DEBUG_ENFORCE
    if (deleteInequalities) delete curFacet[middle];
    curFacet[middle]=0;
    changed=true;
    middle++;
  };

  // Stop here if there are no more inequalities in this quadrant.
  if (middle==stop) return changed;

  // Secondly, delete inequalities if the box and the next inequality
  // define a smaller space.
#ifdef DEBUG_ENFORCE
  cerr << "sweeping upwards from " << middle << " to " << stop 
       << ", testing inequality" << endl;
#endif // DEBUG_ENFORCE
  while (middle+1!=stop && 
	 curFacet[middle]->includes(*curFacet[middle+1], xBound, yBound)) {
#ifdef DEBUG_ENFORCE
    cerr << "deleting (" << middle << ") " << *curFacet[middle] << endl;
#endif // DEBUG_ENFORCE
    if (deleteInequalities) delete curFacet[middle];
    curFacet[middle]=0;
    changed=true;
    middle++;
  };

  // The stop index has to be one beyond the last valid one (and
  // middle indices to the last valid one).
  middle--;

  // Thirdly, iterate from stop to middle. Check if any of the
  // inequalities does not include the box.
  size_t upper=stop-1;

  // Stop here if there are no more inequalities in this quadrant.
  if (upper==middle) return changed;

#ifdef DEBUG_ENFORCE
  cerr << "sweeping downwards from " << upper << " to " << middle 
       << ", testing box" << endl;
#endif // DEBUG_ENFORCE
  while (upper!=middle && curFacet[upper]->includes(xBound, yBound)) {
#ifdef DEBUG_ENFORCE
    cerr << "deleting (" << upper << ") " << *curFacet[upper] << endl;
#endif // DEBUG_ENFORCE
    if (deleteInequalities) delete curFacet[upper];
    curFacet[upper]=0;
    changed=true;
    upper--;
  };

  // Stop here if there are no more inequalities in this quadrant.
  if (upper==middle) return changed;

#ifdef DEBUG_ENFORCE
  cerr << "sweeping downwards from " << upper << " to " << middle 
       << ", testing inequality" << endl;
#endif // DEBUG_ENFORCE
  while (upper-1!=middle &&
	 curFacet[upper]->includes(*curFacet[upper-1],xBound, yBound)) {
#ifdef DEBUG_ENFORCE
    cerr << "deleting (" << upper << ") " << *curFacet[upper] << endl;
#endif // DEBUG_ENFORCE
    if (deleteInequalities) delete curFacet[upper];
    curFacet[upper]=0;
    changed=true;
    upper--;
  };
  return changed;
}

// Remove all inequalities that are redundant with repect to the given
// bounding box. One or both intervals are ignored if the `changed'
// flag in the interval is not set, hence if both intervals are marked
// not being modified, the function is a no-op. The given intervals
// must be as tight as possible, i.e. a combination of an inequality
// from the polyhedron and one bound may not imply a tighter value on
// the other bound. This condition can be enforced by calling
// propagateBounds on this projection.
template<bool isZ>
void PolyhedronImpl<isZ>::tightenAll(const Interval<isZ>& xBound, 
				const Interval<isZ>& yBound,
				vector<int>& isNewFlags) {
  if (facet.empty()) return;
  Inequality** curFacet= &facet.front();
  size_t curSize=facet.size();

  bool noChange=true;
  // Remove inequalities that are no longer feasible with respect to
  // the upper and lower bound on x and y.
  if (tightenQuadrant(curFacet, true, xBound, yBound, 
		      part[east], part[north])) noChange=false;
  if (tightenQuadrant(curFacet, true, xBound, yBound,
		      part[north], part[west])) noChange=false;
  if (tightenQuadrant(curFacet, true, xBound, yBound,
		      part[west], part[south])) noChange=false;
  if (tightenQuadrant(curFacet, true, xBound, yBound,
		      part[south], part[total])) noChange=false;

  // Remove the holes (NULL entries) in the facet array.
  if (noChange) return;
  // There is at least one NULL entry. Find the first one.
  size_t lower=0;
  while (curFacet[lower]) lower++;
  // From here on, copy all pointers down.
  size_t upper=lower;
  // The first NULL entry is in upper.
  do {
    while (curFacet[upper]==NULL) if (++upper==curSize) goto end_reached;
    if (!isNewFlags.empty()) isNewFlags[lower]=isNewFlags[upper];
    curFacet[lower++]=curFacet[upper++];
  } while (upper!=curSize);
 end_reached:
  // lower contains the new number of valid inequalities.
  curSize=lower;

  facet.resize(curSize);
  // Update the partitioning information.
  part=Partitioning(curSize, curFacet);
}

// Update the x bound with respect to the lower y bound.
template<bool isZ>
bool PolyhedronImpl<isZ>::propagateLowerYBound(Interval<isZ>& xBound,
					  const Interval<isZ>& yBound) {
  assert(yBound.lowerIsFinite());

  bool changed=false;

  // Sweep from east to north, updating the upper x bound if necessary.
  size_t current=part[east];
  size_t stop=part[north];
  while (current!=stop) {
    Inequality* e = facet[current++];
    // If the y value of the intersection point does not update the y
    // bound then no other inequalities later in this loop will. Bail
    // out in this case.
    if (!xBound.updateUpper((e->getC() - e->getB() * yBound.getLower())/
			    e->getA())) break;
    changed=true;
  };
#ifdef DEBUG_ENFORCE
    cerr << "x Bound after east - north sweep: " << xBound << endl;
#endif // DEBUG_ENFORCE
  // Sweep from west to north, updating the lower x bound if necessary.
  current=part[west];
  stop=part[north];
  while (current!=stop) {
    Inequality* e = facet[--current];
    // If the y value of the intersection point does not update the y
    // bound then no other inequalities later in this loop will. Bail
    // out in this case.
    if (!xBound.updateLower((e->getC() - e->getB() * yBound.getLower())/
			    e->getA())) break;
    changed=true;
  };
#ifdef DEBUG_ENFORCE
    cerr << "x Bound after west - north sweep: " << xBound << endl;
#endif // DEBUG_ENFORCE
  return changed;
}

// Update the x bound with respect to the upper y bound.
template<bool isZ>
bool PolyhedronImpl<isZ>::propagateUpperYBound(Interval<isZ>& xBound,
					  const Interval<isZ>& yBound) {
  assert(yBound.upperIsFinite());

  bool changed=false;

  // Sweep from west to south, updating the lower x bound if necessary.
  size_t current=part[west];
  size_t stop=part[south];
  while (current!=stop) {
    Inequality* e = facet[current++];
    if (!xBound.updateLower((e->getC() - e->getB() * yBound.getUpper())/
			    e->getA())) break;
    changed=true;
  };
  // Sweep from east to south, updating the upper x bound if necessary.
  current=part[total];
  stop=part[south];
  while (current!=stop) {
    Inequality* e = facet[--current];
    if (!xBound.updateUpper((e->getC() - e->getB() * yBound.getUpper())/
			    e->getA())) break;
    changed=true;
  };
  return changed;
}

// Update the y bound with respect to the lower x bound.
template<bool isZ>
bool PolyhedronImpl<isZ>::propagateLowerXBound(const Interval<isZ>& xBound,
					       Interval<isZ>& yBound) {
  assert(xBound.lowerIsFinite());

  bool changed=false;

  // Sweep from south to east, updating the lower y bound if necessary.
  size_t current=part[south];
  size_t stop=part[total];
  while (current!=stop) {
    Inequality* e = facet[current++];
    if (!yBound.updateLower((e->getC() - e->getA() * xBound.getLower())/
			    e->getB())) break;
    changed=true;
  };
  // Sweep from north to east, updating the upper y bound if necessary.
  current=part[north];
  stop=part[east];
  while (current!=stop) {
    Inequality* e = facet[--current];
    if (!yBound.updateUpper((e->getC() - e->getA() * xBound.getLower())/
			    e->getB())) break;
    changed=true;
  };
  return changed;
}

// Update the y bound with respect to the upper x bound.
template<bool isZ>
bool PolyhedronImpl<isZ>::propagateUpperXBound(const Interval<isZ>& xBound,
					       Interval<isZ>& yBound) {
  assert(xBound.upperIsFinite());

  bool changed=false;

  // Sweep from north to west, updating the upper y bound if necessary.
  size_t current=part[north];
  size_t stop=part[west];
  while (current!=stop) {
    Inequality* e = facet[current++];
    if (!yBound.updateUpper((e->getC() - e->getA() * xBound.getUpper())/
			    e->getB())) break;
    changed=true;
  };
  // Sweep from south to west, updating the lower y bound if necessary.
  current=part[south];
  stop=part[west];
  while (current!=stop) {
    Inequality* e = facet[--current];
    if (!yBound.updateLower((e->getC() - e->getA() * xBound.getUpper())/
			    e->getB())) break;
    changed=true;
  };
  return changed;
}

// Propagate upper and lower bounds implied by the inequalities in
// this polyhedron to the given bounds.
template<bool isZ>
void PolyhedronImpl<isZ>::tightenBounds(Interval<isZ>& xBound,
					Interval<isZ>& yBound) {
  size_t quad1=part[north]-part[east];
  size_t quad2=part[west]-part[north];
  size_t quad3=part[south]-part[west];
  size_t quad4=part[total]-part[south];
  if (quad4>0 && quad1>0)
    xBound.updateUpper(extremeX<isZ>(true, *facet[part[total]-1], *facet[part[east]]));
  if (quad1>0 && quad2>0)
    yBound.updateUpper(extremeY<isZ>(true, *facet[part[north]-1], *facet[part[north]]));
  if (quad2>0 && quad3>0)
    xBound.updateLower(extremeX<isZ>(false, *facet[part[west]-1], *facet[part[west]]));
  if (quad3>0 && quad4>0)
    yBound.updateLower(extremeY<isZ>(false, *facet[part[south]-1], *facet[part[south]]));
}

template<bool isZ> vector<int> PolyhedronImpl<isZ>::isNew;

// Add a set of constraints to this polyhedron.
template<>
void PolyhedronImpl<false>::insertInequalities(Interval<false>& xBound,
					       Interval<false>& yBound,
					       vector<Inequality*>& ineq,
					       const size_t ineqStart) {
  // The facet vectors do not change, hence use a quicker way to access them.
  size_t sizeNew = ineq.size()-ineqStart;
  Inequality** facetNew = (sizeNew ? &ineq[ineqStart] : 0);
  size_t sizeCur = facet.size();
  Inequality** facetCur = (sizeCur ? &facet.front() : 0);

  // Create a temporary array that can hold the inequalities in this
  // polyhedron and the new ones.
  Inequality* temporaryFacet[sizeCur+sizeNew];
  Inequality** tempFacet = temporaryFacet;

#ifdef DEBUG_MEMORY
    bzero(tempFacet, (sizeCur+sizeNew)*sizeof(Inequality*));
#endif // DEBUG_MEMORY

  // Make sure the inequalities are really sorted.
  assert(inequalitiesSorted(sizeNew, facetNew));

  // Setup an array to keep track of which entries in tempFacet are
  // from this polyhedron (entry: false) and which entries are from
  // the passed-in set of inequalities (entry: true). Each operation
  // on tempFacet is mirrored in tempFacetIsNew.
  int temporaryFacetIsNew[sizeCur+sizeNew];
  int* tempFacetIsNew = temporaryFacetIsNew;

  //  Do a merge of the two sorted sequences into the array tempFacet.
  size_t idxCur=0;
  size_t idxNew=0;
  size_t tempLast=0;
  while (idxCur<sizeCur && idxNew<sizeNew) {
    int angle = facetCur[idxCur]->compareWith(*facetNew[idxNew]);
    if (angle==0) {
      // The next two inequalities have the same angle, choose the one
      // which defines a subspace of the other. In case both define the
      // the same halfspace, keep the original in the hope the marking the
      // new inequality as redundant means less work later on.
      if (facetNew[idxNew]->includes(*facetCur[idxCur])) {
	tempFacetIsNew[tempLast] = false;
	tempFacet[tempLast++] = facetCur[idxCur++];
	delete facetNew[idxNew++];
      } else {
	tempFacetIsNew[tempLast] = true;
	tempFacet[tempLast++] = facetNew[idxNew++];
	delete facetCur[idxCur++];
      };
    } else if (angle>0) {
      tempFacetIsNew[tempLast] = true;
      tempFacet[tempLast++] = facetNew[idxNew++];
    } else {
      tempFacetIsNew[tempLast] = false;
      tempFacet[tempLast++] = facetCur[idxCur++];
    }
  };
  // Copy the remaining inequalities.
  while (idxNew<sizeNew) {
    tempFacetIsNew[tempLast] = true;
    tempFacet[tempLast++] = facetNew[idxNew++];
  };
  while (idxCur<sizeCur) {
    tempFacetIsNew[tempLast] = false;
    tempFacet[tempLast++] = facetCur[idxCur++];
  };

  // Reset the size of the input vector to zero.
  ineq.resize(ineqStart);

#ifdef DEBUG_REDUNDANT
  cerr << "bounds and inequalities after merge:\n";
  cerr << "x in " << xBound << ", y in " << yBound << endl;
  for (size_t i=0; i<tempLast; i++)
    cerr << "(" << i << "): " << *tempFacet[i] 
	 << " direction " << tempFacet[i]->calcDirection()
	 << " is " << (tempFacetIsNew[i] ? "new" : "old") << endl;
#endif // DEBUG_REDUNDANT

  // Remove redundant inequalities.
  removeRedundant(xBound, yBound, tempLast,
		  tempFacet, tempFacetIsNew);
#ifdef DEBUG_REDUNDANT
  cerr << "bounds and inequalities after redundancy removal:\n";
  cerr << "x in " << xBound << ", y in " << yBound << endl;
  if (tempLast==0) cerr << "<none>" << endl; else 
    for (size_t i=0; i<tempLast; i++)
      cerr << "(" << i << "): " << *tempFacet[i] 
	   << " direction " << tempFacet[i]->calcDirection()
	   << " is " << (tempFacetIsNew[i] ? "new" : "old") << endl;
#endif // DEBUG_REDUNDANT

  facet.clear();
  facet.reserve(tempLast);
  for(size_t current=0; current<tempLast; current++) {
    facet.push_back(tempFacet[current]);
    if (tempFacetIsNew[current]) ineq.push_back(tempFacet[current]);
  };    

  // Update the partitioning information.
  part=Partitioning(facet);
}

template<>
void PolyhedronImpl<true>::insertInequalities(Interval<true>& xBound,
					      Interval<true>& yBound,
					      vector<Inequality*>& ineq,
					      size_t ineqStart) {
  // The facet vectors do not change, hence use a quicker way to access them.
  size_t sizeNew = ineq.size()-ineqStart;
  Inequality** facetNew = (sizeNew ? &ineq[ineqStart] : 0);
  size_t sizeCur = facet.size();

  // Create a temporary array that can hold the inequalities in this
  // polyhedron and the new ones.
  Inequality* temporaryFacet[sizeCur+sizeNew];
  Inequality** tempFacet = temporaryFacet;

#ifdef DEBUG_MEMORY
  bzero(tempFacet, (sizeCur+sizeNew)*sizeof(Inequality*));
#endif // DEBUG_MEMORY

  // Make sure the inequalities are really sorted.
  assert(inequalitiesSorted(sizeNew, facetNew));

  // Setup an array to keep track of which entries in tempFacet are
  // from this polyhedron (entry: false) and which entries are from
  // the passed-in set of inequalities (entry: true). Each operation
  // on tempFacet is mirrored in tempFacetIsNew.
  int temporaryFacetIsNew[sizeCur+sizeNew];
  int* tempFacetIsNew = temporaryFacetIsNew;

  //  Do a merge of the two sorted sequences into the array tempFacet.
  size_t idxCur=0;
  size_t idxNew=0;
  size_t tempLast=0;
  while (idxCur<sizeCur && idxNew<sizeNew) {
    int angle = facet[idxCur]->compareWith(*facetNew[idxNew]);
    if (angle==0) {
      // The next two inequalities have the same angle, choose the one
      // which defines a subspace of the other. In case both define the
      // the same halfspace, keep the original in the hope marking the
      // new inequality as redundant means less work later on.
      if (facetNew[idxNew]->includes(*facet[idxCur])) {
	tempFacetIsNew[tempLast] = false;
	tempFacet[tempLast++] = facet[idxCur++];
	delete facetNew[idxNew++];
      } else {
	tempFacetIsNew[tempLast] = true;
	tempFacet[tempLast++] = facetNew[idxNew++];
	delete facet[idxCur++];
      };
    } else if (angle>0) {
      tempFacetIsNew[tempLast] = true;
      tempFacet[tempLast++] = facetNew[idxNew++];
    } else {
      tempFacetIsNew[tempLast] = false;
      tempFacet[tempLast++] = facet[idxCur++];
    }
  };
  // Copy the remaining inequalities.
  while (idxNew<sizeNew) {
    tempFacetIsNew[tempLast] = true;
    tempFacet[tempLast++] = facetNew[idxNew++];
  };
  while (idxCur<sizeCur) {
    tempFacetIsNew[tempLast] = false;
    tempFacet[tempLast++] = facet[idxCur++];
  };

  // Reset the size of the input vector to zero.
  ineq.resize(ineqStart);

#ifdef DEBUG_REDUNDANT
  cerr << "bounds and inequalities after merge:\n";
  cerr << "x in " << xBound << ", y in " << yBound << endl;
  for (size_t i=0; i<tempLast; i++)
    cerr << "(" << i << "): " << *tempFacet[i] 
	 << " direction " << tempFacet[i]->calcDirection()
	 << " is " << (tempFacetIsNew[i] ? "new" : "old") << endl;
#endif // DEBUG_REDUNDANT

  // Remove redundant inequalities.
  bool hasChanged=removeRedundant(xBound, yBound, tempLast,
				  tempFacet, tempFacetIsNew);

#ifdef DEBUG_INTEGRALHULL
  cerr << "bounds and inequalities before calculating cuts:\n";
  cerr << "x in " << xBound << ", y in " << yBound << endl;
  if (tempLast==0) cerr << "<none>" << endl; else 
    for (size_t i=0; i<tempLast; i++)
      cerr << "(" << i << "): " << *tempFacet[i] 
	   << " direction " << tempFacet[i]->calcDirection()
	   << " is " << (tempFacetIsNew[i] ? "new" : "old") << endl;
#endif // DEBUG_INTEGRALHULL

  // There is no need to generate cuts if all new inequalities were
  // redundant.
  if (!hasChanged) {
    // However, we need to copy the inequalities from tempFacet to
    // facet since some of the original inequalities may have been
    // redundant. If they happened to have the same angle as some of
    // the new inequalities, then they were not removed by
    // 'removeRedundant' (but during the merge earlier) and
    // 'hasChanged' would be false even though facet contains stale
    // pointers.
    facet.resize(tempLast);
    for (size_t i=0; i<tempLast; i++) facet[i]=tempFacet[i];
    part=Partitioning(facet);
    return;
  };

  // Shuffle the inequalities form the tempFacet array into the facet
  // list of this polyhedron. Within each quadrant insert integral
  // cuts that possibly replace the inequalities in tempFacet.
  isNew.clear();
  facet.clear();

  // Create cuts in each quadrant.
  Partitioning tempPart=Partitioning(tempLast, tempFacet);
  for (Direction curQuad=east; curQuad<total; curQuad=Direction(curQuad+1)) {
    Direction nextQuad=Direction(curQuad+1);
    size_t qStart=tempPart[curQuad];
    size_t qEnd=tempPart[nextQuad];
    if (nextQuad==total) nextQuad=east;

    // Ensure there is at least one inequality in this quadrant to tighten.
    if (qStart==qEnd) continue;
#ifdef DEBUG_INTEGRALHULL
    cerr << "tightening quadrant " << curQuad << ", indices "
	 << qStart << " to " << qEnd << endl;
#endif // DEBUG_INTEGRALHULL

    // Track the last emitted inequality for comparison purposes. Key
    // property of this inequality is that it is never redundant. It
    // points to the second last inequality in the output queue except
    // in the beginning where it might point to eBound.
    Inequality* prevCut=NULL;

    // Create a first cut at the bound by turning the bound into an
    // inequality.
    Inequality eBound;
    switch (curQuad) {
    case east: if (xBound.upperIsFinite()) {
      eBound=Inequality(1, 0, xBound.getUpper().get_num());
      prevCut=&eBound;
    }; break;
    case north: if (yBound.upperIsFinite()) {
      eBound=Inequality(0, 1, yBound.getUpper().get_num());
      prevCut=&eBound;
    }; break;
    case west: if (xBound.lowerIsFinite()) {
      eBound=Inequality(-1,0,-xBound.getLower().get_num());
      prevCut=&eBound;
    }; break;
    case south: if (yBound.lowerIsFinite()) {
      eBound=Inequality(0,-1,-yBound.getLower().get_num());
      prevCut=&eBound;
    }; break;
    default: {
      assert(false);
    }; break;
    }

    // Calculate the first definitive cut stemming from the bound.
    Inequality* cut=NULL;
    if (prevCut) cut=prevCut->calculateCut(*tempFacet[qStart]);

    size_t current=qStart;

    // Move the first integral inequality to the output list.
    if (cut) {
#ifdef DEBUG_INTEGRALHULL
      cerr << "created first cut " << *cut << endl;
#endif // DEBUG_INTEGRALHULL

      // A new cut should never imply a tighter bound since the bounds
      // inferred in the redundancy removal function should be as tight as
      // possible.
      assert(!cut->tightenBounds(xBound, yBound));

      facet.push_back(cut);
      isNew.push_back(true);
    } else { 
#ifdef DEBUG_INTEGRALHULL
      cerr << "reused inequality " << *tempFacet[qStart] << " as first cut" 
	   << endl;
#endif // DEBUG_INTEGRALHULL
      cut=tempFacet[current];
      facet.push_back(cut);
      isNew.push_back(tempFacetIsNew[current]);
      current++;
    };

    while (current<qEnd) {
      assert(facet.size()==isNew.size());

      size_t next=current+1;
      // If possible, look at three inequalities (cut,
      // tempFacet[current] and tempFacet[next]) in sequence. In case
      // there is no third inequality, use the next bound.
      if (next<qEnd ? 
	  tempFacet[current]->includes(*cut, *tempFacet[next]) :
	  tempFacet[current]->includes(*cut, xBound, yBound)) {
	// Replace the redundant current inequality by a tighter cut.
#ifdef DEBUG_INTEGRALHULL
	cerr << "delete " << *tempFacet[current];
#endif // DEBUG_INTEGRALHULL

	delete (tempFacet[current]);
	tempFacet[current]=(next<qEnd ? 
			    cut->calculateCut(*tempFacet[next]) :
			    cut->calculateCut(xBound, yBound, nextQuad));
	tempFacetIsNew[current]=true;
#ifdef DEBUG_INTEGRALHULL
	if (tempFacet[current]) 
	  cerr << " and replace by " << *tempFacet[current] << endl;
	else 
	  cerr << " without replacement" << endl;
#endif // DEBUG_INTEGRALHULL
	if (!tempFacet[current]) current++;
      } else {
	// Calculate more cuts in front of the current inequality.
	cut=cut->calculateCut(*tempFacet[current]);
#ifdef DEBUG_INTEGRALHULL
	if (!cut) cerr << "no ";
	cerr << "new cut in front of " 
	     << (tempFacetIsNew[current] ? "new" : "old")
	     << *tempFacet[current];
	if (cut) cerr << ":" << *cut;
#endif // DEBUG_INTEGRALHULL
	bool cutIsNew=true;
	if (!cut) {
	  cut=tempFacet[current];
	  cutIsNew=tempFacetIsNew[current];
	  current++;
	} else {
          // A new cut should never imply a tighter bound since the bounds
          // inferred in the redundancy removal function should be as tight as
          // possible.
          assert(!cut->tightenBounds(xBound, yBound));
	};
	// The new cut replaces the top inequality in the output if
	// the top inequality is redundant.
	if (prevCut && facet.back()->includes(*prevCut, *cut)) {
#ifdef DEBUG_INTEGRALHULL
	  cerr << " repl. " << (isNew.back() ? "new " : "old ")
	       << *facet.back() << endl;
#endif // DEBUG_INTEGRALHULL
	  delete(facet.back());
	  facet.back()=cut;
	  isNew.back()=cutIsNew;
	} else {
#ifdef DEBUG_INTEGRALHULL
	  cerr << " conf. " << (isNew.back() ? "new " : "old ")
	       << *facet.back() << endl;
#endif // DEBUG_INTEGRALHULL
	  prevCut=facet.back();
	  facet.push_back(cut);
	  isNew.push_back(cutIsNew);
	}
      }
    }
#ifdef DEBUG_INTEGRALHULL
    cerr << "finishing quadrant: bounds x : " << xBound
	 << " y : " << yBound << endl;
#endif // DEBUG_INTEGRALHULL

    while(1) {
      Inequality* newCut=cut->calculateCut(xBound, yBound, nextQuad);
#ifdef DEBUG_INTEGRALHULL
      if (newCut) cerr << "finishing quadrant with new cut " 
		       << *newCut << endl;
#endif // DEBUG_INTEGRALHULL
      if (!newCut) break;
      if (newCut->includes(*cut, xBound, yBound)) {
	delete(newCut);
	break;
      };

      // A new cut should never imply a tighter bound since the bounds
      // inferred in the redundancy removal function should be as tight as
      // possible.
      assert(!newCut->tightenBounds(xBound, yBound));

#ifdef DEBUG_INTEGRALHULL
      cerr << "adding cut " << *newCut;
#endif // DEBUG_INTEGRALHULL
      if (prevCut && facet.back()->includes(*prevCut, *newCut)) {
#ifdef DEBUG_INTEGRALHULL
	cerr << " repl. " << *facet.back() << endl;
#endif // DEBUG_INTEGRALHULL
	delete(facet.back());
	facet.back()=newCut;
	isNew.back()=true;
      } else {                                         
#ifdef DEBUG_INTEGRALHULL
	cerr << " conf. " << *facet.back() << endl;
#endif // DEBUG_INTEGRALHULL
	prevCut=facet.back();
	facet.push_back(newCut);
	isNew.push_back(true);
      }
      cut=newCut;
    }
  };
  
  assert(facet.size()==isNew.size());

  // Update the partitioning information.
  size_t facetSize=facet.size();
  part=Partitioning(facet);
#ifdef DEBUG_INTEGRALHULL
  cerr << "tightening resulted in " << facetSize
       << " inequalities, partitioned " << part << endl;
#endif // DEBUG_INTEGRALHULL

  // Remove inequalities that became redundant with respect to the bounds.
  tightenAll(xBound, yBound, isNew);
#ifdef DEBUG_INTEGRALHULL
  if (facet.size()!=facetSize)
    cerr << "removed " << facetSize-facet.size()
	 << " inequalities due to stricter bound." << endl;
#endif // DEBUG_INTEGRALHULL
  facetSize=facet.size();
  isNew.resize(facetSize);

#ifdef DEBUG_INTEGRALHULL
  cerr << "bounds and final inequalities:\n";
  cerr << "x in " << xBound << ", y in " << yBound 
       << ", partitioning " << part << endl;
#endif // DEBUG_INTEGRALHULL
  // Copy all inequalities that are new to the vector that was passed in.
  ineq.reserve(ineqStart+facetSize);
  for(size_t current=0; current<facetSize; current++) {
    if (isNew[current]) ineq.push_back(facet[current]);
#ifdef DEBUG_INTEGRALHULL
    cerr << "(" << current << "): " << *facet[current] 
	 << " direction " << facet[current]->calcDirection() 
	 << " is " << (isNew[current] ? "new" : "old") << endl;
#endif // DEBUG_INTEGRALHULL
  };    

}

// Remove redundant inequalities. The filter function in our original
// TVPI paper a) uses a doubly linked list and b) is broken. An
// inequality is non-redundant if it is not entailed by its two
// non-redundant neighbours. The test is easy except that it is not
// known that the neigbours we test with are non-redundant. Assume
// that the chain e_r_0..e_r_m is a maximal redundant chain of
// inequalities. It can be shown that at least one of the enailment
// checks {e_(r_0-1), e_r_1} |= e_r_0 or {e_r_(m-1), e_(r_m+1)} |=
// e_r_m will hold which makes it clear that e_r_0 or e_r_(m-1),
// respectively, is redundant. By induction, the whole chain can be
// eliminated by testing the start and the end of the chain
// repeatedly. Observe, however, that in general it is not possible to
// eliminate the whole chain by traversing the set of inequalities
// once: In case the loop starts in the middle of the chain and the
// last element of the chain does not seem to be redundant, the loop
// will only start to remove inequalities when it encounters the
// beginning of the chain. To obtain a linear process, we choose the
// following strategy: Given n inequalities, iterate n steps in search
// of the first redundant inequality in the array. This could be the
// end of a chain, so backtrack and test the previous inequality. This
// backtracking step might take us beyond the beginning of the array. 
// Stop at the first inequality that we cannot prove redundant. Now iterate
// forwards again, applying the same rules until we reach the last
// inequality which we could not prove non-redundant. Thus after a
// maximum of two round trips, we will manifest that the set is
// non-redundant. The return value is true if the original system did
// not change, that is, all new inequalities were redundant and all
// old inequalities were non-redundant.
template<bool isZ>
bool PolyhedronImpl<isZ>::removeRedundant(Interval<isZ>& xBound, 
					  Interval<isZ>& yBound,
					  size_t& tempLast,
					  Inequality**& tempFacet,
					  int*& tempFacetIsNew) {
  bool changed=false;
  // Iterate upwards through the array, move entries down if we remove
  // inequalities.
  size_t prev=tempLast-1;
  size_t current=0;
  // index keeps track of how many inequalities we still need to examine
  for(size_t index=tempLast; index>0; index--) {
    // Stop right here if the facet list is empty.
    if (tempLast==0) break;

    // let "next" point to the following inequality
    size_t next = (current + 1) % tempLast;
    // For easier access, put the relevant inequalities in new variables.
    Inequality* prevF=tempFacet[prev];
    Inequality* currentF=tempFacet[current];
    Inequality* nextF=tempFacet[next];
    assert(prevF);
    assert(currentF);
    assert(nextF);

#ifdef DEBUG_REDUNDANT
    cerr << "#" << index 
	 << " prev(" << prev << "):" << *prevF
	 << " currrent(" << current << "):" << *currentF
	 << " next(" << next << "):" << *nextF << endl;
#endif // DEBUG_REDUNDANT

    // Calculate the directions of the underlying inequalities while
    // considering wrapping of the indices, i.e. only if the indices
    // prev, current, next are increasing they constitue non-wrapped
    // directions.
    int pDir=prevF->calcDirection();
    if (prev>=current) pDir-=4;
    Direction cDir=currentF->calcDirection();
    int nDir=nextF->calcDirection();
    if (current>=next) nDir+=4;
    int prevCurrent = int(cDir)-pDir;
    int currentNext = nDir-int(cDir);

#ifdef DEBUG_REDUNDANT
    cerr << "prev -> current: " << prevCurrent 
	 << ", current -> next: " << currentNext << endl;
#endif // DEBUG_REDUNDANT

    if (currentNext==0) {
      if (prevCurrent==0) {
	// Keep the inequality if the current facet cuts off the
	// intersection point of the other two.
	if (currentF->includes(*prevF,*nextF))
	  goto redundant;
      } else {
	assert(prevCurrent>0 && currentNext==0);
	// The previous inequality is in one of the previous
	// quadrants. The current and next inequality are in the same
	// quadrant. The current inequality might be redundant with
	// regards to the box with the edge chopped off by the next
	// inequality.
	if (currentF->includes(*nextF, xBound, yBound))
	  goto redundant;
      }
    } else {
      // The next inequality is in a different quadrant.
      assert(currentNext>0);

      // There is potential that the current and the next inequality
      // intersect in a point that can reduce a bound of the
      // intervals. Only the (current, next) and not on the (prev,
      // current) pair is examined for an effect on the bounds. This
      // is sufficient since every inequality will be the current one
      // once.
      if (currentNext==1) {
	// The current and the next inequality are in adjacent
	// quadrants. This implies that their intersection point can
	// reduce the bound that lies anglewise between them.
	switch (cDir) {
	case east: 
	  yBound.updateUpper(extremeY<isZ>(true, *currentF, *nextF)); break;
	case north:
	  xBound.updateLower(extremeX<isZ>(false, *currentF, *nextF)); break;
	case west: 
	  yBound.updateLower(extremeY<isZ>(false, *currentF, *nextF)); break;
	case south:
	  xBound.updateUpper(extremeX<isZ>(true, *currentF, *nextF)); break;
	default: break;
	}
      } else if (currentF->lessThanPi(*nextF)) {
	// Make sure the loop will go round at least once more since
	// tightening the boundary might make the next inequality
	// redundant.
	if (index==1) index=2;
	// The current and the next inequality intersect in an extreme
	// point. Shrink the appropriate bound of the box to coincide
	// with this extreme point.
	switch (cDir) {
	case east: {
	  bool xChanged = 
	    xBound.updateLower(extremeX<isZ>(false, *currentF, *nextF));
	  bool yChanged =
	    yBound.updateUpper(extremeY<isZ>(true, *currentF, *nextF));
	  // The bound values are at most those of the point. In case
	  // one of the bounds was tighter than what the intersection
	  // point suggests, the other bounds needs to be tightened
	  // with the intersection of the appropriate inequality and
	  // the tighter bound. Nothing needs to be done if
	  // the point has reduced both bounds, except in the integral
	  // case where the bounds are rounded and thus may require
	  // a tightening with the bounds intersection. See
	  // tests/testPolyhedronWiden1.cpp.
	  if (isZ || !yChanged || !xChanged) {
      	    xBound.updateLower(nextF->calculateX(yBound.getUpper()));
	    yBound.updateUpper(currentF->calculateY(xBound.getLower()));
          }
	}; break;
	case north: {
	  bool xChanged =
	    xBound.updateLower(extremeX<isZ>(false, *currentF, *nextF));
	  bool yChanged =
	    yBound.updateLower(extremeY<isZ>(false, *currentF, *nextF));
	  if (isZ || !yChanged || !xChanged) {
	    xBound.updateLower(currentF->calculateX(yBound.getLower()));
	    yBound.updateLower(nextF->calculateY(xBound.getLower()));
          }
	}; break;
	case west: {
	  bool xChanged =
	    xBound.updateUpper(extremeX<isZ>(true, *currentF, *nextF));
	  bool yChanged =
	    yBound.updateLower(extremeY<isZ>(false, *currentF, *nextF));
	  if (isZ || !yChanged || !xChanged) {
	    yBound.updateLower(currentF->calculateY(xBound.getUpper()));
	    xBound.updateUpper(nextF->calculateX(yBound.getLower()));
          }
	}; break;
	case south: {
	  bool xChanged =
	    xBound.updateUpper(extremeX<isZ>(true, *currentF, *nextF));
	  bool yChanged =
	    yBound.updateUpper(extremeY<isZ>(true, *currentF, *nextF));
	  if (isZ || !yChanged || !xChanged) {
	    yBound.updateUpper(nextF->calculateY(xBound.getUpper()));
	    xBound.updateUpper(currentF->calculateX(yBound.getUpper()));
          }
	}; break;
	default: assert(false);
	};
	// Calculates cuts from currentF to nextF, possibly using one of the bounds.
	// The resulting bounds will be the correct integral bounds.
	if (isZ) refineBoundsToZ(xBound, yBound, currentF, nextF);

#ifdef DEBUG_REDUNDANT
	cerr << "updating direction " << cDir
	     << " bounds: x: " << xBound << " y: " << yBound << endl;
#endif // DEBUG_REDUNDANT
	
      } else {
	    
	// The angle between the current and the next inequality is
	// greater or equal to pi. Update the next bound with respect
	// to the next but one bound.
	currentF->tightenBounds(xBound, yBound);

	// Deal with the special case where the
	// current inequality intersects with the previous but one
	// bound and thus provides an new bound on the previous bound.
	switch (cDir) {
	case east: {
	  if (yBound.lowerIsFinite())
	    xBound.updateUpper(currentF->calculateX(yBound.getLower()));
	}; break;
	case north: {
	  if (xBound.upperIsFinite())
	    yBound.updateUpper(currentF->calculateY(xBound.getUpper()));
	}; break;
	case west: {
	  if (yBound.upperIsFinite())
	    xBound.updateLower(currentF->calculateX(yBound.getUpper()));
	}; break;
	case south: {
	  if (xBound.lowerIsFinite())
	    yBound.updateLower(currentF->calculateY(xBound.getLower()));
	}; break;
	default: assert(false);
	};
	// The symmetric case: with respect to the next inequality,
	// update the previous bound with the intersection point of
	// the inequality and the previous but one bound.
	switch (nextF->calcDirection()) {
	case east: {
	  if (yBound.lowerIsFinite())
	    xBound.updateUpper(nextF->calculateX(yBound.getLower()));
	}; break;
	case north: {
	  if (xBound.upperIsFinite())
	    yBound.updateUpper(nextF->calculateY(xBound.getUpper()));
	}; break;
	case west: {
	  if (yBound.upperIsFinite())
	    xBound.updateLower(nextF->calculateX(yBound.getUpper()));
	}; break;
	case south: {
	  if (xBound.lowerIsFinite())
	    yBound.updateLower(nextF->calculateY(xBound.getLower()));
	}; break;
	default: assert(false);
	};

	currentF->tightenBounds(xBound, yBound);

#ifdef DEBUG_REDUNDANT
	cerr << "updating bound due to empty prev quadrant of " << cDir
	     << ", bounds: x: " << xBound << " y: " << yBound << endl;
#endif // DEBUG_REDUNDANT

	// Ensure the bounds are integral bounds.
	if (isZ && currentF->areEquality<isZ>(*nextF))
	  refineBoundsToZ(xBound, yBound, currentF, nextF);

      };

      if (prevCurrent==0) {
	// The current inequality is the last in this quadrant.
	assert(prevCurrent==0 && currentNext>0);
	
	// Since the previous inequality is in this quadrant, too,
	// the current might be redundant with respect to the
	// bounding box and the previous inequality.
	if (currentF->includes(*prevF, xBound, yBound))
	  goto redundant;
      } else {
	// The previous and next inequalities are more than pi/2 apart.
	assert(currentNext>0 && prevCurrent>0);

	// The current inequality is the only one in this
	// quadrant. One of the corners of the bounding box is the
	// feasible point closest to this inequality. The inequality
	// has to chop off this point in order to be non-redundant.
	if (currentF->includes(xBound, yBound)) goto redundant;
      }
    };

    // The current inequality is non-redundant. If the inequality came
    // from the new set, remember this.
    if (tempFacetIsNew[current]) changed=true;

#ifdef DEBUG_REDUNDANT
    cerr << " currrent(" << current << "):" << *currentF
	 << " non redundant" << endl;
#endif // DEBUG_REDUNDANT
    // shuffle the non-redundant inequality down
    prev = (prev + 1) % tempLast;
    tempFacetIsNew[prev]=tempFacetIsNew[current];
    tempFacet[prev]=tempFacet[current];
    // debugging aid: set the free entry to NULL
    assert((current!=prev ? tempFacet[current]=NULL, true : true));
    // In case current is about to wrap we have to restrict the size
    // of the array to all non-zero entries.
    if (++current==tempLast) {
      current=0;
      tempLast=prev+1;
    };
    continue;

 redundant:
#ifdef DEBUG_REDUNDANT
    cerr << " current(" << current << "):" << *currentF
	 << " is redundant" << endl;
#endif // DEBUG_REDUNDANT
    // The current inequality is redundant, remove it. 
    delete currentF;

    // If the inequality came from the old set, remember this.
    if (!tempFacetIsNew[current]) changed=true;

    // Since this inequality was redundant the previously tested might
    // now be redundant, too. Set the index to the amount of valid
    // inequalities in the array so that the loop iterates over every
    // other inequality up to and including the previous
    // inequality. Note that (current-prev % tempLast) is the number
    // of inequalities that have already been removed (including this
    // one). Add one to cater for the decrement of the index at the
    // end of the loop.
    index=1+tempLast-((current-prev+tempLast) % tempLast);
#ifdef DEBUG_REDUNDANT
    cerr << "setting index to " << index << endl;
#endif // DEBUG_REDUNDANT

    // Backtracking: We differ between two cases: the beginning of
    // the array is filled with inequalities (prev<current) or the
    // array starts with NULL entries.
    if (prev<current) {
      // Instead of setting "current" to the previous value "prev",
      // we move the last entry up to where "current" points to. The
      // only thing left to do is to decrement prev.
      tempFacetIsNew[current]=tempFacetIsNew[prev];
      tempFacet[current]=tempFacet[prev];
      assert((tempFacet[prev]=NULL, true)); // debugging aid: insert NULL
    } else {
      // Whenever prev wraps to the end of the array (due to
      // backtracking) the NULL prefix of the array is chopped off,
      // thus current=0. Thus:
      assert(current==0);
      assert(prev==tempLast-1);
      // There is nothing we can move up if we delete the head of
      // the array. So chop off the leading NULL entry.
      tempFacetIsNew++;
      tempFacet++;
      tempLast--;
    };
    
    // Decrement prev. In case prev is wrapping to the end of the
    // array we need to shorten the array in the beginning by moving
    // the array pointer up.
    if (prev) {
      prev--;
    } else { 
      // chop off all leading NULL entries in the array
      tempFacetIsNew+=current;
      tempFacet+=current;
      // update indices
      tempLast-=current;
      prev=tempLast-1;
      current=0;
    };
  };
  return changed;
}

// Calculate the resultants of this polyhedron and a set of
// inequalities.
template<bool isZ>
void PolyhedronImpl<isZ>::calcResultants(
    bool commonThis,
    const Interval<isZ>& boundThis,
    const Inequality* const ineqOther[],
    const Partitioning& partOther,
    bool commonOther,
    const Interval<isZ>& boundOther,
    vector<Inequality*>& ineqOut,
    bool firstOut) const {
  if (facet.empty()) return;
  size_t sizeOut=ineqOut.size();
  if (firstOut)
    calcResultantsOfSets(ineqOther,
			 partOther,
			 commonOther,
			 boundOther,
			 &facet.front(),
			 part,
			 commonThis,
			 boundThis,
			 ineqOut);
  else
    calcResultantsOfSets(&facet.front(),
			 part,
			 commonThis,
			 boundThis,
			 ineqOther,
			 partOther,
			 commonOther,
			 boundOther,
			 ineqOut);
  if (sizeOut==ineqOut.size()) return;
  sortAndRemoveQuasiSyntacticRed(sizeOut, ineqOut);
}

// Swap the variables of the polyhedron around. 
template<bool isZ>
void PolyhedronImpl<isZ>::swapVars() {
  assert(part[total]==facet.size());
  assert(part[total]>0);
  assert(inequalitiesSorted(part[total], &facet.front()));
  
  for (size_t i=0; i<facet.size(); i++) {
    facet[i]->swapVars();
  };

  Inequality** facetP = &facet.front();

  // Swapping the variables around destroys the order on the angles of
  // the inequalities. Reversing inequalities in the range from east
  // to north and from north to the end resurrects the order.
  Inequality** beg=&facetP[part[east]];
  Inequality** end=&facetP[part[north]-1];

  // Flip [east..northWest)].
  while (beg<end) {
    Inequality* tmp=*beg;
    *beg=*end;
    *end=tmp;
    beg++;
    end--;
  };
    
  beg=&facetP[part[north]];
  end=&facetP[part[total]-1];

  // Flip [northWest..total).
  while (beg<end) {
    Inequality* tmp=*beg;
    *beg=*end;
    *end=tmp;
    beg++;
    end--;
  };

  // Check that we reshuffled correctly.
  assert(inequalitiesSorted(part[total], &facet.front()));

  // Update the partitioning info.
  part=Partitioning(facet);
}

// Check if all inequalities in this polyhedron include the other
// polyhedron which is cut off at the given bounds.
template<bool isZ>
bool PolyhedronImpl<isZ>::includes(const Interval<isZ>& xBound,
				   const Interval<isZ>& yBound,
				   const PolyhedronImpl& other) const {
  // Speed up access to the vector.
  size_t tSize=facet.size();
  Inequality* const* tfacet=(tSize ? &facet.front() : 0);

  // The way the reduced product between intervals and inequalities is
  // implemented allows for redundant inequalities with respect to the
  // bounds. Create an array to hold all inequalities in the other
  // polyhedron that are non-redundant with respect to the given
  // bounds.
  Inequality* tempFacet[other.facet.size()];
  
  size_t facetIdx=0; // =other.part[east]; -- same thing
  size_t tempFacetSize=0;
  size_t tempFacetStart=0;
  while (facetIdx<other.part[north])
    tempFacet[tempFacetSize++]=other.facet[facetIdx++];
  removeNullEntries(end_reached1);
  
  tempFacetStart=tempFacetSize;
  while (facetIdx<other.part[west])
    tempFacet[tempFacetSize++]=other.facet[facetIdx++];
  removeNullEntries(end_reached2);
  
  tempFacetStart=tempFacetSize;
  while (facetIdx<other.part[south])
    tempFacet[tempFacetSize++]=other.facet[facetIdx++];
  removeNullEntries(end_reached3);
  
  tempFacetStart=tempFacetSize;
  while (facetIdx<other.part[total])
    tempFacet[tempFacetSize++]=other.facet[facetIdx++];
  removeNullEntries(end_reached4);
  
  assert(facetIdx==other.facet.size());

#ifdef DEBUG_ENTAILMENT
  if (tempFacetSize!=other.facet.size()) {
    cerr << "removed " << (other.facet.size()-tempFacetSize)
	 << " inequalities before entailment check, leaving" << endl;
    for (size_t i=0; i<tempFacetSize; i++) cerr << tempFacet[i] << endl;
  };
#endif // DEBUG_ENTAILMENT    
  
  // In case there are no inequalities in other polyhedron, inclusion
  // can be checked by verifying that all inequalities in this
  // polyhedron include the bounding box.
  if (tempFacetSize==0) {
    for(size_t i=0; i<tSize; i++) 
      if (!tfacet[i]->includes(xBound, yBound)) return false;
    return true;
  };
  
  // For each inequality in the this polyhedron, look for a pair of
  // adjacent inequalities in the other polyhedron such that the first
  // has a lower angle and the second one an equal or higher
  // angle.
  size_t u=0;
  assert(tempFacetSize>0);
  size_t l=tempFacetSize-1;
  // Assume that l has a direction which is lagging one rotation.
  int lDir=int(tempFacet[l]->calcDirection())-4;
  for(size_t t=0; t<tSize; t++) {
    while (u<tempFacetSize && *tempFacet[u]<*tfacet[t]) {
      l=u++;
      lDir=int(tempFacet[l]->calcDirection());
    };

    // Keep a wrapped version of u.
    size_t uW=u;
    int uDir;
    if (u==tempFacetSize) {
      uW=0;
      uDir=4;
    } else uDir = int(tempFacet[u]->calcDirection());

    // Calculate the quadrant of the this facets to decide how
    // to check for entailment. 
    int tDir=int(tfacet[t]->calcDirection());

    // Based on this difference, distinguish four entailment checks.
    size_t tu = (uDir-tDir);
    size_t lt = (tDir-lDir);

#ifdef DEBUG_ENTAILMENT
    cerr << lt << ":" << tu << "(" << tDir << ")" << *tfacet[t] << " incl. ("
	 << lDir << ")" << *tempFacet[l] << " & ("
	 << uDir << ")" << *tempFacet[uW] << endl;
#endif // DEBUG_ENTAILMENT    

    if (tu==0) if (lt==0) {
      assert(tu==0 && lt==0);
      // All three inequalities are in the same quadrant. The ith
      // inequality of this polyhedron includes the other
      // polyhedron iff the intersection point of the following two
      // inequalities is in its feasible region.
      if (!tfacet[t]->includes(*tempFacet[l],
			       *tempFacet[uW])) return false;
    } else {
      assert(tu==0 && lt>0);
      // The ith inequality lies in the same quadrant as the upper
      // inequality. Intersect the upper inequality with the bounding
      // box and check if that geometric object is entailed.
      if (!tfacet[t]->includes(*tempFacet[uW], 
			       xBound, yBound)) return false;
    } else if (lt==0) {
      assert(tu>0 && lt==0);
      // As above, now for the lower inequality of the pair.
      if (!tfacet[t]->includes(*tempFacet[l], 
			       xBound, yBound)) return false;
    } else {
      assert(tu>0 && lt>0);
      // The ith inequality of this polyhedron is in a different
      // quadrant than the closest inequalities in the other
      // polyhedron. The point that comes closest to saturating the
      // ith inequality is thus a corner of the box.
      if (!tfacet[t]->includes(xBound, yBound)) return false;
    }
  };
  return true;
}

  // Return the number of dimensions of the space that this polyhedron
  // describes. Returns 0, 1 or 2.
  template<bool isZ>
  int PolyhedronImpl<isZ>::dimensionality(const Interval<isZ>& x,
				     const Interval<isZ>& y) const {
    if (x.isSingleton()) return (y.isSingleton() ? 0 : 1);
    if (y.isSingleton()) return 1;
    if (facet.size()==2 && facet[0]->areEquality<isZ>(*facet[1])) return 1;
    return 2;
  }


  // Extrapolate changes in this polyhedron with respect to the given
  // other polyhedron (and its bounds). The extrapolate parameter
  // detemines how many times the change should be extrapolated. The
  // resulting set of inequalities may have redundancies and imply
  // tighter bounds. The partitioning information is invalid on return
  // of this function.
  template<bool isZ>
  void PolyhedronImpl<isZ>::widen(const PolyhedronImpl& other,
				  const mpz_class& extrapolate,
				  bool keepDifferent) {
    // The representation of polyhedra in the reduced product of TVPI
    // and intervals is unequivocal. Hence the revised Halbwachs and
    // the original Halbwachs widening have the same effect so that the
    // latter suffices. Thus widening can be implemented by simply
    // removing all inequalities which are not also in the other
    // polyhedron. Proceed by finding a matching inequality in the other
    // polyhedron for each inequality in this polyhedron.
#ifdef DEBUG_WIDEN
    if (keepDifferent) cerr << "(keepDiff) ";
    cerr << "widening " << *this << "with respect to " << other << endl;
#endif // DEBUG_WIDEN

    // Create two indices for this polyhedron, one to store inequalities
    // and one to read them, writeIdx <= readIdx holds all the time.
    size_t readIdx=0;
    size_t writeIdx=0;
    // The index into the other polyhedron which is entailed by this
    // polyhedron.
    size_t otherIdx=0;
    while (readIdx<facet.size()) {
      // Find a pair of inequalities in this and the other polyhedron
      // with the same angle.
      do {
	// For the given inequality in this polyehdron find the
	// inequality in the other (entailed) polyhedron with at
	// least the same slope.
	while (otherIdx<other.facet.size() &&
	       *other.facet[otherIdx]<*facet[readIdx]) otherIdx++;

#ifdef DEBUG_WIDEN
	cerr << "indices: readIdx=" << readIdx << ", writeIdx=" << writeIdx
	     << ", otherIdx=" << otherIdx << endl;
	if (otherIdx!=other.facet.size())
	  cerr << "picking (" << otherIdx << ") " << *other.facet[otherIdx]
	       << " in previous polyhedron" << endl;
#endif // DEBUG_WIDEN
	// For the given inequality in the other polyhedron, find an
	// inequaility with the same angle in this polyhedron. All
	// inequalities that are skipped in this loop do not pair up
	// with any inequalities in the other polyhedron, hence they are
	// different and we discard them.
	while (otherIdx==other.facet.size() ||
	       *facet[readIdx]<*other.facet[otherIdx]) {
	  if (keepDifferent) {
#ifdef DEBUG_WIDEN
	    cerr << "  smaller angle (" << readIdx << ") " << *facet[readIdx]
		 << " in new polyhedron kept" << endl;
#endif // DEBUG_WIDEN

	    facet[writeIdx++]=facet[readIdx++];
	  } else {
#ifdef DEBUG_WIDEN
	    cerr << "  smaller angle (" << readIdx << ") " << *facet[readIdx]
		 << " in new polyhedron removed" << endl;
#endif // DEBUG_WIDEN
	    delete facet[readIdx++];
	  }
	  // Stop here if there are no more constraints to compare.
	  if (readIdx==facet.size()) goto stop;
	};
	
	// We should have bailed out if we exhausted the inequalities in
	// the other polyhedron.
	assert(otherIdx!=other.facet.size());
	
	// At this point, the angle of the inequality at readIdx
	// might have become larger than that of the other inequality
	// under otherIdx. If this is the case, we need to find a new
	// otherIdx to compare against.
      } while (*facet[readIdx]>*other.facet[otherIdx]);
      
#ifdef DEBUG_WIDEN
      cerr << "  (" << readIdx << ") " << *facet[readIdx]
	   << " parallel (" << otherIdx << ") " << *other.facet[otherIdx]
	   << " ";
#endif // DEBUG_WIDEN
      // We have found two inequalities with the same angle.
      assert(*facet[readIdx]==*other.facet[otherIdx]);
      
      // The entailment realtionship between two inequalities
      // with the same angle can easily be determined by comparing the
      // constants.
      if (other.facet[otherIdx]->getC()!=facet[readIdx]->getC()) {
	assert(other.facet[otherIdx]->getC()<facet[readIdx]->getC());
	// The inequality currently pointed to by readIdx entails less
	// space than the corresponding one of the other
	// polyhedron.

	if (sgn(extrapolate)<=0) {
	  // In case we should extrapolate an arbitrary amount, this
	  // inequality is simply widened away.
#ifdef DEBUG_WIDEN
	  cerr << "discarded" << endl;
#endif // DEBUG_WIDEN
	  delete facet[readIdx++];
	} else {
	  // Extrapolate the change of this inequality.
	  mpz_class diff = (facet[readIdx]->getC()-
			    other.facet[otherIdx]->getC());
	  facet[readIdx]->getC()+=diff*extrapolate;
#ifdef DEBUG_WIDEN
	  cerr << "changed inequality: " << *facet[readIdx] << endl;
#endif // DEBUG_WIDEN
	  facet[writeIdx++]=facet[readIdx++];
	}
      } else {
#ifdef DEBUG_WIDEN
	cerr << "kept" << endl;
#endif // DEBUG_WIDEN
	facet[writeIdx++]=facet[readIdx++];
      }
    };
  stop:
    facet.resize(writeIdx);
  }
  
  // Return the maximum integral value in this polyhedron that lies in
  // the direction given by the coefficients of the inequality. If the
  // found value is finite, the function returns true and set the
  // constant of the inequality result. In case the polyhedron extends
  // to infinity towards the given direction, the function returns
  // false.
  template<bool isZ>
  bool PolyhedronImpl<isZ>::linOpt(const Interval<isZ>& xBound,
				   const Interval<isZ>& yBound,
				   Inequality& goal) const {

    // Copy'n'Paste from extreme().
    size_t tempFacetSize=facet.size()+4;
    Inequality* tempFacet[tempFacetSize];

#ifdef DEBUG_LINOPT
    cerr << "linOpt: called on " << facet.size() << " equations and x="
	 << xBound << ", y=" << yBound << endl;
#endif // DEBUG_LINOPT
    
    // Fill the array with the bounds and the inequalities in each
    // partition.
    Inequality boundEast, boundNorth, boundWest, boundSouth;
    tempFacetSize=0;
    size_t facetIdx=0; // =part[east]; -- same thing
    if (xBound.upperIsFinite()) {
      boundEast=Inequality(xBound.getUpper().get_den(), mpz_class(0), xBound.getUpper().get_num());
      tempFacet[tempFacetSize++]=&boundEast;
    };
    size_t tempFacetStart=tempFacetSize;
    while (facetIdx<part[north]) tempFacet[tempFacetSize++]=facet[facetIdx++];
    removeNullEntries(end_reached1);
    if (yBound.upperIsFinite()) {
      boundNorth=Inequality(mpz_class(0), yBound.getUpper().get_den(), yBound.getUpper().get_num());
      tempFacet[tempFacetSize++]=&boundNorth;
    };
    tempFacetStart=tempFacetSize;
    while (facetIdx<part[west]) tempFacet[tempFacetSize++]=facet[facetIdx++];
    removeNullEntries(end_reached2);
    if (xBound.lowerIsFinite()) {
      boundWest=Inequality(-xBound.getLower().get_den(), mpz_class(0), -xBound.getLower().get_num());
      tempFacet[tempFacetSize++]=&boundWest;
    };
    tempFacetStart=tempFacetSize;
    while (facetIdx<part[south]) tempFacet[tempFacetSize++]=facet[facetIdx++];
    removeNullEntries(end_reached3);
    if (yBound.lowerIsFinite()) {
      boundSouth=Inequality(mpz_class(0), -yBound.getLower().get_den(), -yBound.getLower().get_num());
      tempFacet[tempFacetSize++]=&boundSouth;
    };
    tempFacetStart=tempFacetSize;
    while (facetIdx<part[total]) tempFacet[tempFacetSize++]=facet[facetIdx++];
    removeNullEntries(end_reached4);

    assert(facetIdx==facet.size());
    assert(tempFacetSize<=facet.size()+4);

#ifdef DEBUG_LINOPT
    cerr << "linOpt: searching " << goal << " in:" << endl;
    for (size_t i = 0; i < tempFacetSize; i++)
      cerr << "(" << i << ") : " << *tempFacet[i] << endl;
#endif // DEBUG_LINOPT
      
    // No finite value if the polyhedron is unrestricted.
    if (tempFacetSize==0) return false;

    // Do a binary search for an inequality that has the same or
    // larger angle. Result is in 'lower'.
    size_t lower = 0;
    size_t upper = tempFacetSize;
    while (lower<upper) {
      size_t mid = (lower+upper)/2;
      if (*tempFacet[mid]<goal) lower=mid+1; else upper=mid;
    };

    size_t larger = lower;
    if (larger==tempFacetSize) larger=0;

#ifdef DEBUG_LINOPT
    cerr << "linOpt: found index: " << larger;
#endif // DEBUG_LINOPT

    // Check for direct hit.
    if (goal==*tempFacet[larger]) {
#ifdef DEBUG_LINOPT
      cerr << " which is exact angle" << endl;
#endif // DEBUG_LINOPT
      goal.getC() = tempFacet[larger]->getC();
      return true;
    };

    size_t smaller = (larger ? larger-1 : tempFacetSize-1);

    // Check for intersection point.
#ifdef DEBUG_LINOPT
    cerr << " which is larger and " << smaller << " which is smaller" << endl;
#endif // DEBUG_LINOPT

    // Don't heed a single inequality.
    if (smaller==larger) return false;

    // The result is the intersection point of the two adjacent
    // inequalities.
    if (tempFacet[smaller]->lessThanPi(*tempFacet[larger])) {
      Point p = Point(*tempFacet[smaller], *tempFacet[larger]);
      goal.getC()=goal.getA()*p.getX()+goal.getB()*p.getY();
      return true;
    };
    return false;
  }

  template<bool isZ>
  std::ostream& output(std::ostream& stream,
		       const PolyhedronImpl<isZ>& p,
		       const DomVar varNameX,
		       const DomVar varNameY) {
    using namespace std;
    if (varNameX==0 && varNameY==0) {
      if (p.refCount>1) stream << p.refCount << "* shared ";
      stream << "polyhedron partitioned " << p.part
	     << " inequalities:" << endl;
    };
    for(size_t current=0; current<p.facet.size(); current++) {
      if (varNameX==0 && varNameY==0) stream << "e" << current << ": ";
      p.facet[current]->output(stream, varNameX, varNameY);
      size_t next=current+1;
      if (next==p.facet.size()) next=0;
      if (p.facet[current]->lessThanPi(*p.facet[next]))
	stream << " intersects at " << Point(*p.facet[current], *p.facet[next])
	       << (next ? " with" : " first");
      stream << endl;
    };
    return stream;
  }

// Polyhedron: a reference-counting wrapper around PolyehdronImpl

// Instead of generating tons of new empty polyhedra we create an
// empty polyhedron once and copy that each time a new, emtpy
// polyhedron is created. In case the memory leak tracking is on,
// calling the class' new operation would trigger a write operation
// which crashes on FreeBSD. Hence, use the global new operator.
template<bool isZ>
PolyhedronImpl<isZ>*
  Polyhedron<isZ>::emptyPolyhedronImpl = ::new PolyhedronImpl<isZ>();

// Create an empty polyhedron.
template<bool isZ>
Polyhedron<isZ>::Polyhedron() {
  // Make a shallow copy of the null polyhedron.
  poly=emptyPolyhedronImpl;
  poly->refCount++;
};

// Create a duplicate of this polyhedron (copy constructor).
template<bool isZ>
Polyhedron<isZ>::Polyhedron(const Polyhedron& p1) {
  poly=p1.poly;
  poly->refCount++;
};

template<bool isZ>
Polyhedron<isZ>::~Polyhedron() {
  assert(poly);
  assert(poly->refCount>0);
  if (!--poly->refCount) {
    if (poly->twin) {
      assert(poly->twin->twin==poly);
      if (poly->twin->refCount==0) {
	poly->twin->twin=0;
	delete poly->twin;
	poly->twin=0;
	delete poly;
      }
    }
  }
};

template<bool isZ>
Polyhedron<isZ>& Polyhedron<isZ>::operator=(const Polyhedron& other) {
  assert(poly);
  assert(poly->refCount>0);
  assert(other.poly);
  other.poly->refCount++;
  if (!--poly->refCount) {
    if (poly->twin) {
      assert(poly->twin->twin==poly);
      if (poly->twin->refCount==0) {
	poly->twin->twin=0;
	delete poly->twin;
	poly->twin=0;
	delete poly;
      }
    }
  }
  poly=other.poly;
  return *this;
};


// Propagate a tighter bound on the x axis to the y axis.
template<bool isZ>
void Polyhedron<isZ>::propagateXBounds(Interval<isZ>& xBound,
				  Interval<isZ>& yBound) {
  if (!xBound.hasChanged()) return;
  makeUnique();
  if (xBound.lowerIsFinite()) poly->propagateLowerXBound(xBound, yBound);
  if (xBound.upperIsFinite()) poly->propagateUpperXBound(xBound, yBound);
}

// Propagate a tighter bound on the y axis to the x axis.
template<bool isZ>
void Polyhedron<isZ>::propagateYBounds(Interval<isZ>& xBound,
				  Interval<isZ>& yBound) {
  if (!yBound.hasChanged()) return;
  makeUnique();
  if (yBound.lowerIsFinite()) poly->propagateLowerYBound(xBound, yBound);
  if (yBound.upperIsFinite()) poly->propagateUpperYBound(xBound, yBound);
}

// Add a set of constraints to this polyhedron. The constraints must
// be sorted by angle, no constrains may have the same angle and no
// constraint may have zero coefficients. The bounds are tightened if
// necessary.
template<bool isZ>
void Polyhedron<isZ>::addInequalitySet(Interval<isZ>& xBound,
				       Interval<isZ>& yBound,
				       vector<Inequality*>& ineq,
				       const size_t ineqStart) {
  assert(poly);
  assert(ineqStart<=ineq.size());
  if (ineqStart==ineq.size()) return;
  makeUnique();
  poly->insertInequalities(xBound, yBound, ineq, ineqStart);
}

// Retrieve the number of inequalities with positive and negative
// coeffiecients. The coefficient which is examined is 'a' if
// secondCoeff is true and 'b' otherwise.
template<bool isZ>
void Polyhedron<isZ>::getNumOfInequalitiesWithSign(bool secondCoeff,
					      size_t* plusCoeff,
					      size_t* minusCoeff) {
  assert(poly);
  poly->part.getNumOfInequalitiesWithSign(secondCoeff,
					  plusCoeff,
					  minusCoeff);
}

// Print a polyhedron with variable names.
template<bool isZ>
std::ostream& output(std::ostream& stream,
		     const Polyhedron<isZ>& p,
		     const DomVar varNameX,
		     const DomVar varNameY) {
  assert(p.poly);
  output(stream, *p.poly, varNameX, varNameY);
  return stream;
}


// Get a pointer to the idx'th inequality. The pointer belongs to
// the polyhedron and may only be used until the next modification
// of this polyhedron.
template<bool isZ>
const Inequality* Polyhedron<isZ>::operator[](int idx) const {
  assert(poly);
  return poly->facet[idx];
}


// Ensure that the underlying PolyhedronImpl object is referenced only by
// this class. This function is called before any destructive actions
// are taken.
template<bool isZ>
void Polyhedron<isZ>::makeUnique() {
  assert(poly);
  assert(poly->refCount>0);
  if (poly->refCount==1) {
    // Loose the twin.
    if (poly->twin) {
      assert(poly->twin->twin==poly);
      poly->twin->twin=0;
      if (poly->twin->refCount==0) delete poly->twin;
      poly->twin=0;
    };
    return;
  };
  // Somebody else is refering to this data, thus make a deep copy. 
  // This gets rid of the twin automatically.
  poly->refCount--;
  poly = new PolyhedronImpl<isZ>(*poly);
  assert(poly->twin==0);
}

defMemDbg(,PolyhedronImpl<isZ>,p,P)

template class PolyhedronImpl<false>;
template class PolyhedronImpl<true>;
template class Polyhedron<false>;
template class Polyhedron<true>;
template std::ostream& output<false>(std::ostream&, const Polyhedron<false>&,
				     const DomVar, const DomVar);
template std::ostream& output<true>(std::ostream&, const Polyhedron<true>&,
				    const DomVar, const DomVar);
template std::ostream& output<false>(std::ostream&,
				     const PolyhedronImpl<false>&,
				     const DomVar, const DomVar);
template std::ostream& output<true>(std::ostream&,
				    const PolyhedronImpl<true>&,
				    const DomVar, const DomVar);


}; // namespace Tvpi
