/* This file is part of the TVPI library, adapted for inclusion in Crab.
 * Original author: Axel Simon <A.Simon@kent.ac.uk>
 * Original source: https://github.com/axelsimon/tvpi
 *
 * This file is distributed under the GNU General Public License version 2.
 * See COPYING in this directory for the full license text.
 *
 * Modifications for Crab (2026): Yusen Su <yusen.su@uwaterloo.ca>
 *   - Adjusted include paths for bundled build
 *   - Added public accessors to Inequality (getA/getB/getC) in planar.hh
 *   - Added public accessors to Polyhedron (getNoOfInequalities,
 *     getInequalityAt) in polyhedron.hh
 */
// polyhedron.hh
// A planar polyhedron represented by a finite number of inequalities.
//
// This file defines two classes, Polyhedron and PolyhedronImpl. The
// Polyhedron class is a wrapper of around PolyhedronImpl and provides
// copy-on-write behaviour for the PolyhedronImpl class.

#ifndef __POLYHEDRON_H
#define __POLYHEDRON_H

#include "common.hh"
#include "planar.hh"
#include "interval.hh"
#include "tvpiexception.hh"
#include "memory.hh"
#include <iostream>
#include <vector>
#include <assert.h>


// TODO:

// - widening could inline linOpt such that the inequalities don't have
// to be copied repeatedly and such that an (log n) lookup can become
// constant time.

namespace Tvpi {


template<bool isZ> class Polyhedron;
template<bool isZ> class PolyhedronImpl;

template<bool isZ>
std::ostream& output(std::ostream& stream,
		     const Polyhedron<isZ>& p,
		     const DomVar varNameX,
		     const DomVar varNameY);

template<bool isZ>
std::ostream& output(std::ostream& stream,
		     const PolyhedronImpl<isZ>& p,
		     const DomVar varNameX,
		     const DomVar varNameY);

// A Partitioning contains information on which inequalities of a
// Polyhedron correspond to which quadrant. This auxilliary class is
// used in caluclating the set of resultants (the set of inequalities
// resulting in eliminating a common variable of two polyehdra).
class Partitioning {
  friend class PolyhedronImpl<true>;
  friend class PolyhedronImpl<false>;

  // The partitions are a set of indices, enumerating all sign
  // combinations of the inequality coefficients in a
  // counter-clockwise sequence. The sign combinations correspond to
  // the normal vectors of the inequalities. The sequence is: <0,->,
  // <+,->, <+,0>, <+,+>, <0,+>, <-,+>, <-,0>, <-,->, <0,0>. The last
  // element corresponds to a tautology which are not allowed. 
    
  // This is the array containing the indices where each partition
  // starts. The total number of inequalities is kept in the last
  // entry which proves to be useful when calculating the resultants.
  size_t dirStart[numDirs];

public:
  // Create an empty partitioning.
  Partitioning() {
    for(int d=east; d<numDirs; d++) dirStart[d]=0;
  };

  // Create partitioning information for a set of inequalities.
  Partitioning(size_t facetSize, Inequality** facet);

  // Create partitioning information for a vector of inequalities.
  Partitioning(const std::vector<Inequality* > facet);

private:

  // Determine the indices where the inequalities change the
  // direction.
  void binarySearch(Inequality* const facet[],
		    Direction dirLower,
		    Direction dirUpper);

public:
  // Retrieve the number of inequalities with positive and negative
  // coeffiecients. The coefficient which is examined is 'a' if
  // secondCoeff is true and 'b' otherwise.
  void getNumOfInequalitiesWithSign(bool secondCoeff,
				    size_t* plusCoeff,
				    size_t* minusCoeff);


public:
  size_t operator[](Direction dir) const {
    assert(dir<numDirs);
    return dirStart[dir];
  };


  friend std::ostream& operator<<(std::ostream& stream,
				  const Partitioning& p) {
    stream << "(";
    for(int current=int(east); current<int(total); current++)
      stream << p.dirStart[current] << ", ";
    stream << p.dirStart[total] << ")";
    return stream;
  };

};

// Calculate the resultants of two sets of inequalities. The set of
// resultants corresponds to eliminating a common variable between the
// first and the second set of inequalities. The common variable is
// the determined by the flags f1CommonCoeff and f2CommonCoeff. If
// f1CommonCoeff is false, then the common variable is x, otherwise it
// is y. Similarly for f2CommonCoeff. The partitioning contains
// information on how many inequalities there are in each set with
// positive and negative coefficients. From the two partitions and the
// knowledge which variable to elimintate, the maximum size of the
// resulatants set can be calculated. This resulting set of
// inequalities is stored in ineqOut which has to be large enough to
// hold all resultants. The variable sizeOut is used as an index into
// the ineqOut array and increases by the number of added resultants
// in ineqOut. It must be initialized to 0 if the array ineqOut is to
// be filled from the start. The inequalities of the output set will
// have the first coefficients from the first input set facet1 and the
// second coefficient will be the remaining coefficients from the
// second input set facet2.  The input set of inequalities must be
// ordered by angle but may contain redundant constraints. The output
// set is unordered and may contain redundant constraints, in
// particular it may contain tautologies and constraints with the same
// angle.
template<bool isZ>
void calcResultantsOfSets(const Inequality* const facet1[],
			  const Partitioning& f1Part,
			  bool f1CommonCoeff,
			  const Interval<isZ>& bound1,
			  const Inequality* const facet2[],
			  const Partitioning& f2Part,
			  bool f2CommonCoeff,
			  const Interval<isZ>& bound2,
			  std::vector<Inequality*>& ineqOut);

// Sort a set of constraints and remove all quasi syntactic redundant
// inequalities (so called by Jaffar). Quasi syntactic redundant
// inequalities are inequalities with the same angle where the
// constant coefficient conveys straight away which of the two
// inequalities is redundant (weaker).
void sortAndRemoveQuasiSyntacticRed(size_t ineqStart,
				    std::vector<Inequality*>& ineq);


//  This class implements the operations on satisfiable polyhedra.
template<bool isZ> 
class PolyhedronImpl {
  friend class Polyhedron<isZ>;
  friend class Partitioning;

  // The container for the inequalities.
  std::vector<Inequality*> facet;

  // The number of Polyhedron objects that use this set of inequalities.
  size_t refCount;

  // The facet array can be partitioned into sequences that have the
  // same Direcition. The Partitioning object contains information on
  // the indices of each of the eight Directions.
  Partitioning part;

  // Optionally, the polyhedron may retain a pointer to another
  // polyhedron that is it's mirror image. This pointer may be 0.
  PolyhedronImpl* twin;

public:
  // Create an empty polyhedron.
  PolyhedronImpl();

  // Copy constructor. Creates a deep copy of the facet array.  Called
  // by class Polyhedron if it wants to modify this polyhedron and
  // refCount is greater than 1.
  PolyhedronImpl(const PolyhedronImpl& p);

  PolyhedronImpl(const PolyhedronImpl& p1,
		 const Interval<isZ>& p1x,
		 const Interval<isZ>& p1y,
		 const PolyhedronImpl& p2,
		 const Interval<isZ>& p2x,
		 const Interval<isZ>& p2y,
		 bool flipFirst) {
    convexHull(p1, p1x, p1y, p2, p2x, p2y, flipFirst);
  };

  // Calculate the convex hull of two non-empty polyhedra. The result
  // is undefined if one of the polyhedra is not satisfiable. Each
  // polyhedron (set of inequalities) must be non-redundant but they
  // may contain redundant inequalities with respect to the given
  // bounds. If the flag flipFirst is set the x and y axes of the
  // first polyhedron (including the passed-in bounds) are swapped
  // around.
  void convexHull(const PolyhedronImpl& p1,
		  const Interval<isZ>& p1x,
		  const Interval<isZ>& p1y,
		  const PolyhedronImpl& p2,
		  const Interval<isZ>& p2x,
		  const Interval<isZ>& p2y,
		  bool flipFirst);

private:
  // Calculate the points and rays of a polyhedron.
  void extreme(const Interval<isZ>& xBound,
	       const Interval<isZ>& yBound,
	       size_t& lastPoint, Point points[],
	       size_t& lastRay, Ray rays[]) const;

public:
  // Destructor.
  ~PolyhedronImpl();

private:
  // Remove inequalities that are redundant with respect to the given
  // bounds. Inequalities between start and stop are examined
  // and their array entry set to NULL if redundant. The removed
  // inequalities are freed if deleteInequalities is set.
  static bool tightenQuadrant(Inequality* curFacet[],
			      bool deleteInequalities,
			      const Interval<isZ>& xBound,
			      const Interval<isZ>& yBound,
			      size_t start,
			      size_t stop);

  // Remove all inequalities that are redundant with repect to the
  // given bounding box. One or both intervals are ignored if the
  // `changed' flag in the interval is not set, hence if both
  // intervals are marked as unmodified, the function is a no-op.
  // The third element is either NULL or an array with elements for
  // each inequality. The array is reshuffled just as the vector if
  // inequalities is.
  void tightenAll(const Interval<isZ>& xBound, const Interval<isZ>& yBound,
		  std::vector<int>& isNewFlags);

  // Propagate relational information from the lower y bound to the x
  // bounds. Sweep over the inequalities in the two quadrants spanning
  // from east over north to west and determine with which (if any)
  // inequality the bound intersects. If an inequality is found, the
  // intersection point determines an upper (east - north quadrant) or
  // lower (north - west quadrant) bound for x. Update the x bound in
  // this case. This function must be called after a bound has changed
  // and before tightenAll is called which removes inequalities that
  // became redundant with respect to the new bounds. The function
  // returns true if any bound has changed.
  bool propagateLowerYBound(Interval<isZ>& xBound,
			    const Interval<isZ>& yBound);

  // Similarly for the other bounds.
  bool propagateUpperYBound(Interval<isZ>& xBound,
			    const Interval<isZ>& yBound);
  bool propagateLowerXBound(const Interval<isZ>& xBound,
			    Interval<isZ>& yBound);
  bool propagateUpperXBound(const Interval<isZ>& xBound,
			    Interval<isZ>& yBound);


public:
  // Propagate upper and lower bounds implied by the inequalities in
  // this polyhedron to the given bounds.
  void tightenBounds(Interval<isZ>& xBound, Interval<isZ>& yBound);

  // Add a set of constraints to this polyhedron. The constraints must
  // be sorted by angle, no constrains may have the same angle and no
  // constraint may have a or b as zero. The vector must contain
  // pointers to heap allocated inequalities. The ownership of the
  // inequalities is transfered to this polyhedron. The vector will
  // contain inequalities starting at ineqStart which are new in this
  // polyhedron. The inequalities in the vector are owned by this
  // polyhedron and may not be freed or modified. This operation may
  // render the polyhedron unsatisfiable. If ineqLast is zero, all
  // inequalities from ineqStart to the end of the vector are
  // inserted. In this case the vector is resized and all inequalities
  // that were new in this polyhedron are copied into the vector.
  void insertInequalities(Interval<isZ>& xBound,
			  Interval<isZ>& yBound,
			  std::vector<Inequality*>& ineq,
			  const size_t ineqStart=0);

private:
  // A dynamically growing vector of flags determining the source of
  // an inequality. The type should be bool, really, but that is
  // specialized to a compacted bit-vector which is no good to pass to
  // removeRedundant.
  static std::vector<int> isNew;

  // Shuffle the inequalities form the tempFacet array into the facet
  // list of this polyhedron. Within each quadrant insert integral cuts
  // that possibly replace the inequalities in tempFacet.
  void integralHull(Interval<isZ>& xBound, Interval<isZ>& yBound,
		    size_t tempLast, Inequality** tempFacet,
		    int* tempFacetIsNew);

  // Remove redundant inequalities. Used by insertInequalities and
  // enforceBounds.
  static bool removeRedundant(Interval<isZ>& xBound, Interval<isZ>& yBound,
			      size_t& tempSize, Inequality**& tempFacet,
			      int*& tempFacetXRef);

public:
  // Calculate the resultants of this polyhedron and a set of
  // inequalities. The set of resultants corresponds to eliminating a
  // common variable between this polyhedron and the set of
  // inequalities. The common variable is the determined by the flags
  // firstThis and firstOther. If firstThis is false, then the common
  // variable is that following a, otherwise it is that following
  // b. Similarly for firstOther. The resulting set of inequalities is
  // stored in ineqOut and sizeOut where firstOut determines whether
  // the first coefficient in the output corresponds to this
  // polyhedron (firstOut=false) or to the other set of inequalities
  // (firstOut=true). The function sets sizeOut to 0 if the set of
  // resultants is empty. In this case the output variable ineqOut is
  // not modified. The input set of inequalities must be ordered by
  // angle but may contain redundant constraints. The output set is
  // totally ordered by angle and may contain redundant constraints.
  void calcResultants(bool commonThis,
		      const Interval<isZ>& boundThis,
		      const Inequality* const ineqOther[],
		      const Partitioning& partOther,
		      bool commonOther,
		      const Interval<isZ>& boundOther,
		      std::vector<Inequality*>& ineqOut,
		      bool firstOut) const;

  // Swap the variables of the polyhedron around.
  void swapVars();

  // Check if all inequalities in this polyhedron include the other
  // polyhedron which is cut off at the given bounds.
  bool includes(const Interval<isZ>& xBound,
		const Interval<isZ>& yBound,
		const PolyhedronImpl<isZ>& other) const;

  // Return the number of dimensions of the space that this polyhedron
  // describes. Returns 0, 1 or 2.
  int dimensionality(const Interval<isZ>& x,
		     const Interval<isZ>& y) const;

  // Check if this polyhedron contains an equality.
  inline bool isEquality() const {
    return (facet.size()==2 && facet[0]->areEquality<isZ>(*facet[1]));
  };

  // Extrapolate changes in this polyhedron with respect to the given
  // other polyhedron (and its bounds). The extrapolate parameter
  // detemines how many times the change should be extrapolated. The
  // resulting set of inequalities may have redundancies and imply
  // tighter bounds. The partitioning information is invalid on return
  // of this function.
  void widen(const PolyhedronImpl<isZ>& other,
	     const mpz_class& extrapolate,
	     bool keepDifferent);

  // Stretch the polyhedron in x direction.
  inline void stretchX(Mult m) {
    for(typename std::vector<Inequality*>::iterator i=facet.begin();
	i!=facet.end(); i++) (*i)->stretchX(m);
  };

  // Stretch the polyhedron in y direction.
  inline void stretchY(Mult m) {
    for(typename std::vector<Inequality*>::iterator i=facet.begin();
	i!=facet.end(); i++) (*i)->stretchX(m);
  };
  
  // Return the maximum integral value in this polyhedron that lies in
  // the direction given by the coefficients of the inequality. If the
  // found value is finite, the function returns true and set the
  // constant of the inequality result. In case the polyhedron extends
  // to infinity towards the given direction, the function returns
  // false.
  bool linOpt(const Interval<isZ>& xBound,
	      const Interval<isZ>& yBound,
	      Inequality& goal) const;

  friend std::ostream& operator<<(std::ostream& stream,
				  const PolyhedronImpl<isZ>& p) {
    output<isZ>(stream, p, 0, 0);
    return stream;
  };

  friend std::ostream& output<isZ>(std::ostream& stream,
				   const PolyhedronImpl<isZ>& p,
				   const DomVar varNameX,
				   const DomVar varNameY);

  declMemDbg;

};


template<bool isZ>
class Polyhedron {
  friend class Partitioning;

  // The (possibly shared) data of this polyhedron.
  PolyhedronImpl<isZ>* poly;

  // A pointer to an empty polyhedron.
  static PolyhedronImpl<isZ>* emptyPolyhedronImpl;

 public:
  // Create an empty polyhedron. 
  Polyhedron();

  // Create a duplicate of this polyhedron (copy constructor).
  Polyhedron(const Polyhedron& p1);

  // Calculate the convex hull of the two planar polyhedpra p1 and p2.
  // Or rather: try to avoid it.
  Polyhedron(const Polyhedron& p1,
	     const Interval<isZ>& p1x,
	     const Interval<isZ>& p1y,
	     const Polyhedron& p2,
	     const Interval<isZ>& p2x,
	     const Interval<isZ>& p2y,
	     bool flipFirst = false) {
    assert(p1.poly);
    assert(p2.poly);
    // Whenever one of the input polyhedra is empty, the result is
    // an empty polyhedron, so just copy it.
    if (p1.poly->facet.size()==0 &&
	p1x.upperIsInfinite() &&
	p1x.lowerIsInfinite() &&
	p1y.upperIsInfinite() &&
	p1y.lowerIsInfinite()) {
      poly=p1.poly;
      poly->refCount++;
      return;
    };
    if (p2.poly->facet.size()==0 &&
	p2x.upperIsInfinite() &&
	p2x.lowerIsInfinite() &&
	p2y.upperIsInfinite() &&
	p2y.lowerIsInfinite()) {
      poly=p2.poly;
      poly->refCount++;
      return;
    };
    if (flipFirst) {
      if (p1.poly->twin && p1.poly->twin==p2.poly &&
	  p1y==p2x && p1x==p2y) {
	// Whoha! If the pointers are the same, the underlying polyhedra
	// must be the same. The convex hull of two identical polyhedra is
	// just one of them.
	poly=p1.poly->twin;
	poly->refCount++;
	// Sanity check: There are at least two Polyhedron objects pointing
	// to the same PolyhedronImpl object (the twin might have a refcount
	// of 0).
	assert(poly->refCount>=2);
	return;
      } else if (p2.poly->twin && p1.poly==p2.poly->twin &&
		 p1y==p2x && p1x==p2y) {
	poly=p1.poly;
	poly->refCount++;
	assert(poly->refCount>=2);
	return;
      } 
    } else if (p1.poly==p2.poly && p1x==p2x && p1y==p2y) {
      poly=p1.poly;
      poly->refCount++;
      assert(poly->refCount>=3);
      return;
    };
 
    // Calculate the convex hull: Pass the two non-empty input polyhedra
    // to the constructor of PolyhedronImpl.
    poly=new PolyhedronImpl<isZ>(*p1.poly, p1x, p1y,
				 *p2.poly, p2x, p2y, flipFirst);
  };

  ~Polyhedron();

  Polyhedron& operator=(const Polyhedron& other);

  // Reduce the feasible space of the current polyhedron by adding
  // more constraints.

  // Remove all inequalities that are redundant with repect to the
  // given bounding box. One or both intervals are ignored if the
  // `changed' flag in the interval is not set, hence if both
  // intervals are marked not being modified, the function is a no-op.
  void enforceBounds(const Interval<isZ>& xBound,
		     const Interval<isZ>& yBound) {
    makeUnique();
    std::vector<int> fake;
    poly->tightenAll(xBound, yBound, fake);
  };

  // Propagate a tighter bound on the x axis to the y axis.
  void propagateXBounds(Interval<isZ>& xBound, Interval<isZ>& yBound);

  // Propagate a tighter bound on the y axis to the x axis.
  void propagateYBounds(Interval<isZ>& xBound, Interval<isZ>& yBound);

  // Add a set of constraints to this polyhedron. The constraints must
  // be sorted by angle, no constrains may have the same angle and no
  // constraint may have zero coefficient. The bounds are tightened if
  // necessary. If ineqLast is zero, all inequalities from ineqStart
  // to the end of the vector are inserted. In this case the vector is
  // resized and all inequalities that were new in this polyhedron are
  // copied into the vector from ineqStart onwards.
  void addInequalitySet(Interval<isZ>& xBound,
			Interval<isZ>& yBound,
			std::vector<Inequality*>& ineq,
			const size_t ineqStart=0);

  // Retrieve the number of inequalities with positive and negative
  // coeffiecients. The coefficient which is examined is 'a' if
  // secondCoeff is true and 'b' otherwise.
  void getNumOfInequalitiesWithSign(bool secondCoeff,
				    size_t* plusCoeff,
				    size_t* minusCoeff);

  // Calculate the resultants of this polyhedron and a set of
  // inequalities. The set of resultants corresponds to eliminating a
  // common variable between this polyhedron and the set of
  // inequalities. The common variable is the determined by the flags
  // commonThis and commonOther. If commonThis is false, then the
  // common variable is x, otherwise it is y. Similarly for
  // commonOther. The resulting set of inequalities is stored in
  // ineqOut and sizeOut where firstOut determines whether the first
  // coefficient in the output corresponds to this polyhedron
  // (firstOut=false) or to the other set of inequalities
  // (firstOut=true). The input set of inequalities must be
  // ordered by angle but may contain redundant constraints. The
  // output set is totally ordered by angle and may contain redundant
  // constraints.
  void calcResultants(bool commonThis,
		      const Interval<isZ>& boundThis,
		      const Inequality* const ineqOther[],
		      const Partitioning& partOther,
		      bool commonOther,
		      const Interval<isZ>& boundOther,
		      std::vector<Inequality*>& ineqOut,
		      bool firstOut) const {
    assert(poly);
    poly->calcResultants(commonThis,
			 boundThis,
			 ineqOther,
			 partOther,
			 commonOther,
			 boundOther,
			 ineqOut,
			 firstOut);
  };

  // Swap the variable of the polyhedron around.
  Polyhedron& swapVars() {
    if (poly->facet.size()==0) return *this;
    if (poly->twin) {
      assert(poly->twin->twin==poly);
      poly->twin->refCount++;
      poly->refCount--;
    } else {
      poly->twin=new PolyhedronImpl<isZ>(*poly);
      poly->twin->twin=poly;
      poly->twin->swapVars();
      poly->refCount--;
      assert(poly->twin->refCount==1);
    };
    poly=poly->twin;
    assert(poly->twin->refCount>=0);
    assert(poly->refCount>0);
    return *this;
  };


  // Check if all inequalities in this polyhedron include the other
  // polyhedron which is cut off at the given bounds.
  bool includes(const Interval<isZ>& xBound,
		const Interval<isZ>& yBound,
		const Polyhedron& other) const {
    assert(poly);
    assert(other.poly);
    return poly->includes(xBound, yBound, *other.poly);
  };

  // Return the number of dimensions of the space that this polyhedron
  // describes. Returns 0, 1 or 2.
  inline int dimensionality(const Interval<isZ>& x,
			    const Interval<isZ>& y) const {
    assert(poly);
    return poly->dimensionality(x, y);
  };

  // Check if this polyhedron contains an equality.
  inline bool isEquality() const {
    assert(poly);
    return poly->isEquality();
  };

  // Widen this polyhedron with respect to the other (smaller)
  // polyhedron. This polyhedron must include the other
  // polyhedron. The bounds of this polyhedron must have been widened
  // already. The extrapolate parameter determines how many times the
  // change should be applied and has to be greater than one (for
  // extrapolation) or 0 for standard widening.
  inline void widen(Interval<isZ>& xBound,
		    Interval<isZ>& yBound,
		    const Interval<isZ>& xBoundOther,
		    const Interval<isZ>& yBoundOther,
		    const Polyhedron& other,
		    const mpz_class& extrapolate=0) {
    assert(poly);
    assert(other.poly);
    if (poly==other.poly) return;
    assert(includes(xBoundOther, yBoundOther, other));
    makeUnique();
    //int dimThis = dimensionality(xBound, yBound);
    //int dimOther = other.dimensionality(xBoundOther, yBoundOther);
    // Extrapolate the change from the previous polyhedron to this
    // polyhedron (which must include this polyhedron). Inequalities
    // in this polyhedron whose angle corresponds to an inequality in
    // the other polyhedron are moved by "extrapolate" times their
    // relative displacement. If "extrapolate" is zero, the inequality
    // in this polyhedron is removed (or displaced by an infinite
    // amount). Inequalities in this polyhedron for which no
    // inequality in the other polyhedron exists that has the same
    // angle are kept if the dimensionality of this and the other
    // polyhedron is different, otherwise they are discarded. Ups: for
    // the moment, we remove all inequalities, otherwise we might
    // extrapolation doesn't extrapolate properly (i.e. the bounds get
    // reduced and then full widening is applied, leading to a
    // precision loss).
    poly->widen(*other.poly, extrapolate, false /*dimThis!=dimOther*/);
    // Moving some constraints out means that they might have become
    // redundant and that their intersection points with adjacent
    // inequalities are no longer integeral. Simply empty the
    // polyhedron and re-insert all inequalities again.
    std::vector<Inequality*> tempFacets;
    tempFacets.swap(poly->facet);
    poly->part=Partitioning();
    poly->insertInequalities(xBound, yBound, tempFacets);
  };

  // Stretch the polyhedron in x and y direction.
  inline const Polyhedron stretch(Mult mx, Mult my) const {
    if (mx==0 && my==0) return *this;
    assert(poly);
    if (poly->facet.size()==0) return *this;
    Polyhedron p = Polyhedron(*this);
    p.makeUnique();
    if (mx>0) p.poly->stretchX(mx);
    if (my>0) p.poly->stretchY(my);
    return p;
  };

  // Return the maximum integral value in this polyhedron that lies in
  // the direction given by the coefficients of the inequality. If the
  // found value is finite, the function returns true and set the
  // constant of the inequality result. In case the polyhedron extends
  // to infinity towards the given direction, the function returns
  // false.
  inline bool linOpt(const Interval<isZ>& xBound,
		     const Interval<isZ>& yBound,
		     Inequality& goal) const {
    assert(poly);
    return poly->linOpt(xBound, yBound, goal);
  };

  // Propagate upper and lower bounds implied by the inequalities in
  // this polyhedron to the given bounds. This function is needed to
  // narrow the bounds after widening the bounds and the inequalities
  // separately.
  void tightenBounds(Interval<isZ>& xBound, Interval<isZ>& yBound) {
    assert(poly);
    makeUnique();
    poly->tightenBounds(xBound, yBound);
  };

  friend std::ostream& operator<<(std::ostream& stream,
				  const Polyhedron<isZ>& p) {
    assert(p.poly);
    return stream << *p.poly;
  };

  friend std::ostream& output<isZ>(std::ostream& stream,
				   const Polyhedron<isZ>& p,
				   const DomVar varNameX,
				   const DomVar varNameY);
    
  void showAngles() {
    using namespace std;
    cout << "e_y<e_x";
    for(size_t current=0; current<poly->facet.size(); current++)
      cout << "\t" << current;
    cout << endl;
    for(size_t y=0; y<poly->facet.size(); y++) {
      cout << y;
      for(size_t x=0; x<poly->facet.size(); x++) {
	bool smaller = (*poly->facet[y]<=*poly->facet[x]);
	cout << "\t" << smaller;
      };
      cout << endl;
    };
  };

  // Query the number of inequalities contained in this polyhedron.
  inline size_t getNoOfInequalities() const {
    assert(poly);
    return poly->facet.size();
  };

  // Get a pointer to the idx'th inequality. The pointer belongs to
  // the polyhedron and may only be used until the next modification
  // of this polyhedron.
  const Inequality* operator[](int idx) const;

  // Query the polyhedron if the content is shared.
  int isShared() const { assert(poly); return poly->refCount; };

private:
  // Ensure that the underlying PolyhedronImpl object referenced only by
  // this class. This function is called before any destructive actions
  // are taken.
  void makeUnique();

public:
  void sane() const {
    assert(poly);
    assert(poly->facet.size()==poly->part[total]);
    for(size_t i=0; i<poly->part[total]; i++) {
      assert(sgn(poly->facet[i]->getA())!=0);
      assert(sgn(poly->facet[i]->getB())!=0);
    }
  };
};

};

#endif // __POLYHEDRON_H
