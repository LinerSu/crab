/* Derived from the TVPI library (original file: tvpi.cpp).
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
// tvpi.cpp - Implementation of the two-variables-per-inequality domain.
//
// The planar projections are kept in a triangular matrix. Each
// projection is accessed by a macro polyhedron(i,j) where i<j. In
// such a projection i corresponds to the x variable while j
// corresponds to the y axis.

#include <crab/domains/tvpi_impl/tvpi.hh>
#include <crab/domains/tvpi_impl/planar.hh>
#include <crab/domains/tvpi_impl/interval.hh>
#include <crab/domains/tvpi_impl/polyhedron.hh>
#include <crab/domains/tvpi_impl/memory.hh>
// affine.hh not needed: we use DenseTvpi directly, not the Domain class
#include <stdlib.h>
#include <vector>

#ifdef DEBUG_MEMORY
#include <string>
#endif // DEBUG_MEMORY

using namespace std;

#undef DEBUG_JOIN
#undef DEBUG_CLOSURE
#undef DEBUG_INCLUDES
#undef DEBUG_INCLUDES_FAIL
#undef DEBUG_AUGMENT
#undef DEBUG_WIDEN
#undef DEBUG_NEWVAR
#undef DEBUG_LINOPT
#undef DEBUG_APPROXIMATE
#undef DEBUG_SUBSTITUTE

#undef WIDEN_CHECK_ENTAILMENT
// TODO:

// - we could triple the code for calculating the resultants for two
// hops and perhaps gain some speed: right now we check if either of
// the way x,j-j,y and x,k-k,y is possible. We can split the test into
// three, if both yield constraints, if x,j yields constraints or if
// x,k yields constraints. That way we could avoid recalulating a
// whole row for one of the detours.

// - add two variable-choice parameters to the convex hull of two
// polyhedra and redo augment and update.

// Min and max functions.
#define min(x,y) (x<y ? x : y)
#define max(x,y) (x>y ? x : y)

// Define a macro to access the polyhedra array for a polyhedron.
#define polyhedronIndex(i,j) \
 (assert(i<j), \
  ((j)*((j)-1)/2+(i)))

#define polyhedronArraySize(el) polyhedronIndex(el-2,el-1)+1

#ifdef NDEBUG
#define polyhedron(i,j) polyhedra[polyhedronIndex(i,j)]
#else
#define polyhedron(i,j) \
  polyhedra.at(polyhedronIndex(i,j))
#endif

// isIntegral is true for the integer TVPI domain (equivalent to template param isZ)
#define isIntegral isZ

namespace Tvpi {

// Create a new Two-Variable-Per-Inequality polyhedron.
template<bool isZ>
DenseTvpi<isZ>::DenseTvpi(size_t dimension) {
  size_t arraySize=1;
  if (dimension>2) arraySize=polyhedronArraySize(dimension);
  // Allocate memory for arraySize entries.
  polyhedra.reserve(arraySize);
  bounds.reserve(dimension);
};

// Take the convex hull of two TVPI systems. The convex hull of a TVPI
// system is calculated by running the planar convex hull algorithm on
// each projection. The perm parameter maps each variable to a
// variable in the other domain. It may not be zero and each position
// must be a valid variable in the other domain.
template<bool isZ>
void DenseTvpi<isZ>::join(const DenseTvpi& other, TvpiVar* perm,
			  Mult* otherMult) {
  size_t size=bounds.size();
  assert(size<=other.bounds.size());
  assert(perm);
  assert(perm[0]!=invalidTvpiVar);

  // Note to self: Do not insert code here that shows the permutation. 
  // Enable DEBUG_JOIN in the caller instead!

  // Calculate the convex hull of the polyhedra.
  size_t proj=0;
  for(size_t y=1; y<size; y++) {
    TvpiVar otherY = perm[y];
    assert(otherY!=invalidTvpiVar);
    assert(0<=otherY);
    assert(otherY<other.bounds.size());
    assert(otherMult[otherY]!=invalidMult);
    for(size_t x=0; x<y; x++) {
      TvpiVar otherX = perm[x];
      assert(0<=otherX);
      assert(otherX!=otherY);
      assert(otherX<other.bounds.size());
      assert(otherMult[otherX]!=invalidMult);
#ifdef DEBUG_JOIN
      cerr << "joining " << (otherX>otherY ? "flipped " :"")
	   << "polyhedron at prev(" << otherX << ", " << otherY
	   << "), 2^" << (int) otherMult[otherX] << " * " 
	   << other.bounds[otherX].interval << ", 2^"
	   << (int) otherMult[otherY] << " * "
	   << other.bounds[otherY].interval  << endl;
      cerr << (otherX>otherY ? other.polyhedron(otherY,otherX) :
	       other.polyhedron(otherX,otherY));
      cerr <<"with *this index "
	   << proj << ", " << bounds[x].interval
	   << ", " << bounds[y].interval << endl;
      cerr << polyhedra[proj];
#endif // DEBUG_JOIN
      if (otherX<otherY) {
	polyhedra[proj]=
	  Polyhedron<isZ>(other.polyhedron(otherX,otherY).
			  stretch(otherMult[otherX],otherMult[otherY]),
			  other.bounds[otherX].interval.
			  stretch(otherMult[otherX]),
			  other.bounds[otherY].interval.
			  stretch(otherMult[otherY]),
			  polyhedra[proj],
			  bounds[x].interval,
			  bounds[y].interval);
      } else {
	polyhedra[proj]=
	  Polyhedron<isZ>(other.polyhedron(otherY,otherX).
			  stretch(otherMult[otherY],otherMult[otherX]),
			  other.bounds[otherY].interval.
			  stretch(otherMult[otherY]),
			  other.bounds[otherX].interval.
			  stretch(otherMult[otherX]),
			  polyhedra[proj],
			  bounds[x].interval,
			  bounds[y].interval,
			  true);
      };

#ifdef DEBUG_JOIN
      cerr << "resulting in " << Interval<isZ>(bounds[x].interval,
					       other.bounds[otherX].interval.
					       stretch(otherMult[otherX]))
	   << ", " << Interval<isZ>(bounds[y].interval,
				    other.bounds[otherY].interval.
				    stretch(otherMult[otherY])) << endl;
      cerr << polyhedra[proj];
#endif // DEBUG_JOIN

      proj++;
    }
  };
  // Calculate the larger of the intervals.
  for(size_t i=0; i<size; i++)
    bounds[i].interval=Interval<isZ>(bounds[i].interval,
				     other.bounds[perm[i]].interval.
				     stretch(otherMult[perm[i]]));
 };

  // Check if the state space of this polyhedron includes that of the
  // other polyhedron. For each variable in this polyhedron, perm
  // denotes the variable in the other polyhedron to which this
  // variable should be compared to. If an index is invalidTvpiVar,
  // then the corresponding variable in this domain is not compared at
  // all.
  template<bool isZ>
  bool DenseTvpi<isZ>::includes(const DenseTvpi& other, TvpiVar* perm) const {
    size_t size=bounds.size();
    assert(perm);

#ifdef DEBUG_INCLUDES
    cerr << "inclusion check: this(" << size << "):\n"
	 << *this << "other(" << other.bounds.size() << "):\n" << other;
#endif // DEBUG_INCLUDES

    size_t proj=0;
    for(size_t y=0; y<size; y++) {
      if (perm[y]!=invalidTvpiVar) {
	assert(bounds[y].xRef==other.bounds[perm[y]].xRef);
	if (!bounds[y].interval.includes(other.bounds[perm[y]].interval)) {
#ifdef DEBUG_INCLUDES_FAIL
	  cerr << "inclusion: false due to bounds of p_"
	       << y << " : " << bounds[y].interval << " does not include "
	       << other.bounds[perm[y]].interval << endl;
#endif // DEBUG_INCLUDES_FAIL
	  return false;
	}
      };
      for(size_t x=0; x<y; x++, proj++) if (perm[x]!=invalidTvpiVar) {
	size_t otherX = perm[x];
	size_t otherY = perm[y];
	if (otherX<otherY) {
	  // Make sure the other polyhedron does not contain any
	  // redundant inequalities.
	  ((DenseTvpi) other).polyhedron(otherX,otherY).
	    enforceBounds(other.bounds[otherX].interval,
			  other.bounds[otherY].interval);
	  if (!polyhedra[proj].includes(other.bounds[otherX].interval,
					other.bounds[otherY].interval,
					other.polyhedron(otherX,otherY))) {
#ifdef DEBUG_INCLUDES_FAIL
	    cerr << "inclusion: false due to p_" << x << ": "
		 << bounds[x].interval << ", p_" << y << " : "
		 << bounds[y].interval << " and " << endl
		 << polyhedra[proj] << " does not include p_"
		 << otherX << ": " << other.bounds[otherX].interval
		 << ", p_" << otherY << " : "
		 << other.bounds[otherY].interval << endl
		 << other.polyhedron(otherX,otherY) << endl;
#endif // DEBUG_INCLUDES_FAIL
	    return false;
	  }
	} else {
	  assert(otherX!=otherY);
	  ((DenseTvpi) other).polyhedron(otherY,otherX).
	    enforceBounds(other.bounds[otherY].interval,
			  other.bounds[otherX].interval);
	  Polyhedron<isZ> p = other.polyhedron(otherY, otherX);
	  if (!polyhedra[proj].includes(other.bounds[otherX].interval,
					other.bounds[otherY].interval,
					p.swapVars())) {
#ifdef DEBUG_INCLUDES_FAIL
	    cerr << "inclusion: false due to p_" << x << ": "
		 << bounds[x].interval << ", p_" << y << " : "
		 << bounds[y].interval << " and " << endl
		 << polyhedra[proj] << " does not include flipped p_"
		 << otherX << ": " << other.bounds[otherX].interval
		 << ", p_" << otherY << " : "
		 << other.bounds[otherY].interval << endl
		 << other.polyhedron(otherY,otherX) << endl;
#endif // DEBUG_INCLUDES_FAIL
	    return false;
	  }
	}
      }
    }
    assert(proj==polyhedra.size());
#ifdef DEBUG_INCLUDES
    cerr << "inclusion held" << endl;
#endif // DEBUG_INCLUDES
    return true;
  }
    

  // Widen this TVPI system with repect to the other given TVPI
  // system. The perm parameter maps each variable in this domain to a
  // variable in the other domain. It may not be zero and each
  // position must be a valid variable in the other domain. The
  // extrapolate value denotes the number of steps that should be
  // extrapolated.
  template<bool isZ>
  void DenseTvpi<isZ>::widen(const DenseTvpi<isZ>& other, TvpiVar* perm,
			     const mpz_class extrapolate) {
    assert(bounds.size()<=other.bounds.size());
    size_t size = bounds.size();

    assert(perm);
#ifdef DEBUG_WIDEN
    cerr << "widening: apply changes of this" << endl
	 << *this << extrapolate << " times to other" << endl << other;
#endif // DEBUG_WIDEN
  
  // To widen this polyhedron, it should define a superset of the
  // other polyhedron.
#ifdef WIDEN_CHECK_ENTAILMENT
  assert(includes(other, perm));
#endif // WIDEN_CHECK_ENTAILMENT

  // Widen the intervals of this polyhedron.
  for(size_t i=0; i<size; i++) {
    assert(perm[i]>=0);
    assert(perm[i]<other.bounds.size());
    bounds[i].interval.widen(other.bounds[perm[i]].interval, extrapolate);
#ifdef DEBUG_WIDEN
    cerr << "var p_" << i << " with p_" << perm[i] << ": "
	 << bounds[i].interval << endl;
#endif // DEBUG_WIDEN
  };

  size_t proj=0;
  for(size_t y=1; y<size; y++)
    for(size_t x=0; x<y; x++, proj++) {
      size_t otherX = perm[x];
      size_t otherY = perm[y];

#ifdef DEBUG_WIDEN
      cerr << "widening: this, x=" << x << ", y=" << y
	   << " with other, x=" << otherX << ", y=" << otherY
	   << ", i.e. proj=" << proj << endl;
#endif // DEBUG_WIDEN

	// Thirdly, widen the inequalities in this polyhedron.
      if (otherX<otherY) {
	polyhedra[proj].widen(bounds[x].interval,
			      bounds[y].interval,
			      other.bounds[otherX].interval,
			      other.bounds[otherY].interval,
			      other.polyhedron(otherX, otherY), extrapolate);
      } else {
	Polyhedron<isZ> p = other.polyhedron(otherY, otherX);
	polyhedra[proj].widen(bounds[x].interval,
			      bounds[y].interval,
			      other.bounds[otherX].interval,
			      other.bounds[otherY].interval,
			      p.swapVars(), extrapolate);
      };
#ifdef DEBUG_WIDEN
      cerr << "resulting in " << polyhedra[proj] << endl;
#endif // DEBUG_WIDEN

    };

  // We must have widened all projections by now.
  assert(proj==polyhedra.size());

}

// Return the maximum integral value in this domain that lies in the
// direction given by the linear expression.
template<bool isZ>
bool DenseTvpi<isZ>::linOpt(const vector<TVPIComponent>& comps,
			    mpz_class& value) const {
  switch (comps.size()) {
  case 0: {
    value=0;
    return true;
  }; break;
  case 1: return linOpt(comps[0], value); break;
  case 2: return linOpt(comps[0], comps[1], value); break;
  case 3: {
    mpz_class val1, val2;
    bool isFinite=false;
    if (linOpt(comps[0], val1) && linOpt(comps[1], comps[2], val2)) {
      value=val1+val2;
      isFinite=true;
    };
#ifdef DEBUG_LINOPT
    cerr << "0,1-2 - value found: ";
    if (isFinite) cerr << value; else cerr << "none";
    cerr << endl;
#endif // DEBUG_LINOPT
    if (linOpt(comps[0], comps[1], val1) && linOpt(comps[2], val2)) {
      if (isFinite) {
	mpz_class sum=val1+val2;
	if (value<sum) mpz_swap(value.get_mpz_t(), sum.get_mpz_t());
      } else {
	value=val1+val2;
	isFinite=true;
      };
    };
#ifdef DEBUG_LINOPT
    cerr << "0-1,2 - value found: ";
    if (isFinite) cerr << value; else cerr << "none";
    cerr << endl;
#endif // DEBUG_LINOPT
    if (linOpt(comps[0], comps[2], val1) && linOpt(comps[1], val2)) {
      if (isFinite) {
	mpz_class sum=val1+val2;
	if (value<sum) mpz_swap(value.get_mpz_t(), sum.get_mpz_t());
      } else {
	value=val1+val2;
	isFinite=true;
      };
    };
#ifdef DEBUG_LINOPT
    cerr << "0-2,1 - value found: ";
    if (isFinite) cerr << value; else cerr << "none";
    cerr << endl;
#endif // DEBUG_LINOPT
    return isFinite;
  }; break;
  case 4: {
    mpz_class val1, val2;
    bool isFinite=false;
    if (linOpt(comps[0], comps[1], val1) && linOpt(comps[2], comps[3], val2)) {
      value=val1+val2;
      isFinite=true;
    };
    if (linOpt(comps[0], comps[2], val1) && linOpt(comps[1], comps[3], val2)) {
      if (isFinite) {
	mpz_class sum=val1+val2;
	if (value<sum) mpz_swap(value.get_mpz_t(), sum.get_mpz_t());
      } else {
	value=val1+val2;
	isFinite=true;
      };
    };
    if (linOpt(comps[0], comps[3], val1) && linOpt(comps[1], comps[2], val2)) {
      if (isFinite) {
	mpz_class sum=val1+val2;
	if (value<sum) mpz_swap(value.get_mpz_t(), sum.get_mpz_t());
      } else {
	value=val1+val2;
	isFinite=true;
      };
    };
  }; break;
  };
  // Just pretend we'd do more than for vars and discover there's no
  // bound.
  return false;
}

  // Query the maximum value of a term. Returns true if a maximum
  // exists and sets value to this maximum.
  template<bool isZ>
    bool DenseTvpi<isZ>::linOpt(const TVPIComponent& single,
				mpz_class& value) const {
    assert(sgn(single.coefficient)!=0);
    if (sgn(single.coefficient)>0) {
      if (bounds[single.variable].interval.upperIsFinite()) {
	value = bounds[single.variable].interval.getUpper();
	value = value * single.coefficient;
	if (single.mult!=invalidMult)
	  mpz_mul_2exp(value.get_mpz_t(), value.get_mpz_t(), single.mult);
#ifdef DEBUG_LINOPT
	cerr << "maximum of " << single.coefficient
	     << bounds[single.variable].interval << " is " << value
	     << " (using upper)" << endl;
#endif // DEBUG_LINOPT
	return true;
      } else return false;
    } else {
      assert(sgn(single.coefficient)<0);
      if (bounds[single.variable].interval.lowerIsFinite()) {
	value = bounds[single.variable].interval.getLower();
	value = value * single.coefficient;
	if (single.mult!=invalidMult)
	  mpz_mul_2exp(value.get_mpz_t(), value.get_mpz_t(), single.mult);
#ifdef DEBUG_LINOPT
	cerr << "maximum of " << single.coefficient
	     << bounds[single.variable].interval << " is " << value
	     << " (using lower)" << endl;
#endif // DEBUG_LINOPT
	return true;
      } else return false;
    }
  }

  // Query the maximum value of two terms. Returns true if a maximum
  // exists and sets value to this maximum.
  template<bool isZ>
    bool DenseTvpi<isZ>::linOpt(const TVPIComponent& first,
				const TVPIComponent& second,
				mpz_class& value) const {
    assert(sgn(first.coefficient)!=0);
    assert(sgn(second.coefficient)!=0);
    Inequality e(first.coefficient, second.coefficient, 0);
    if (first.mult!=invalidMult)
      mpz_mul_2exp(e.getA().get_mpz_t(), e.getA().get_mpz_t(), first.mult);
    if (second.mult!=invalidMult)
      mpz_mul_2exp(e.getB().get_mpz_t(), e.getB().get_mpz_t(), second.mult);
    mpz_class gcd;
    mpz_gcd(gcd.get_mpz_t(),
	    e.getA().get_mpz_t(),
	    e.getB().get_mpz_t());
    if (gcd>1) {
      mpz_divexact(e.getA().get_mpz_t(), e.getA().get_mpz_t(),
		   gcd.get_mpz_t());
      mpz_divexact(e.getB().get_mpz_t(), e.getB().get_mpz_t(),
		   gcd.get_mpz_t());
    };
    assert(first.variable<second.variable);
    bool isFinite =
      polyhedron(first.variable,
		 second.variable).linOpt(bounds[first.variable].interval,
					 bounds[second.variable].interval, e);
    mpz_swap(value.get_mpz_t(),e.getC().get_mpz_t());
    if (gcd>1) value*=gcd;
#ifdef DEBUG_LINOPT
    if (isFinite) 
      cerr << "maximum of " << first.coefficient << "p_" << first.variable
	   << " + " << second.coefficient << "p_" << second.variable
	   << " is " << value << endl;
#endif // DEBUG_LINOPT
    return isFinite;
  }

  // Ensure that the polyhedron has room for at least n more variables.
  template<bool isZ>
    void DenseTvpi<isZ>::makeRoomForVariables(size_t headroom) {
    size_t goalSize=bounds.size()+headroom;
    if (goalSize>2) {
      bounds.reserve(goalSize);
      polyhedra.reserve(polyhedronArraySize(goalSize));
    };
  }


// Remove all variables which are not mentioned in the given array.
template<bool isZ>
void DenseTvpi<isZ>::projectOnto(size_t sizeVarSeq,
				 TvpiVar varSeq[],
				 size_t headroom) {
  assert(sizeVarSeq<=bounds.size());
  // Create an empty new bounds vector that has the required space.
  vector<PolyBound> newBounds;
  newBounds.reserve(sizeVarSeq+headroom);

  // Fill the vector. Note that if intervals are implemented as
  // copy-by-value objects, then this should be improved by a hack
  // into the STL library since right now we are creating sizeVarSeq
  // new intervals before we are deleting all the previous ones.
  for(size_t i=0; i<sizeVarSeq; i++)
    newBounds.push_back(bounds[varSeq[i]]);

  // Create a new space to store the projections, but don't fill it
  // yet.
  Polyhedra newPolyhedra;
  size_t newArraySize=(sizeVarSeq+headroom>2 ?
		       polyhedronArraySize(sizeVarSeq+headroom) : 1);
  newPolyhedra.reserve(newArraySize);
  // For each entry in the new space, pick the right projection from
  // the old system, possibly swapping the variables in that
  // projection.
  for(size_t y=1; y<sizeVarSeq; y++) {
    for(size_t x=0; x<y; x++) {
      if (varSeq[x]==varSeq[y])
	throw IllegalArgument("DenseTvpi::projectOnto",
			      "given map is not a function");
      size_t source;
      if (varSeq[x]<varSeq[y]) {
	// The projection can be copied verbatim.
	source=polyhedronIndex(varSeq[x],varSeq[y]);
      } else {
	// Sine x<y for all loop iterations but in the new system
	// varSeq[x]<varSeq[y] does not hold, the variables in the
	// projection have to be swapped.
	source=polyhedronIndex(varSeq[y],varSeq[x]);
	polyhedra[source].swapVars();
      };
      assert(polyhedronIndex(x,y)==newPolyhedra.size());
      newPolyhedra.push_back(polyhedra[source]);
    }
  };
  swap(newBounds,bounds);
  swap(newPolyhedra,polyhedra);
};

// Update a variable with the variable added most recently. 
template<bool isZ>
void DenseTvpi<isZ>::update(TvpiVar var, bool additive)
  throw (IllegalArgument) {
  size_t curSize=bounds.size();
  if (var>=curSize-1)
    throw IllegalArgument("DenseTvpi::update",
			  "variable number out of range");
  if (curSize<=2) {
    assert(curSize==2);
    // There are only two variables in the system. The result will be
    // a lower and an upper bound on the only variable left.
    if (additive)
      bounds[0].interval=Interval<isZ>(bounds[0].interval, bounds[1].interval);
    else
      bounds[0].interval=bounds[1].interval;
    bounds[0].xRef=bounds[1].xRef;
    bounds.resize(1);
    polyhedra.resize(0);
  } else {
    // The last row of the triangular matrix resides between these
    // boundaries.
#ifndef NDEBUG
    size_t upper=polyhedronArraySize(curSize);
#endif
    // After this line curSize is the number of the last variable in
    // the old system.
    curSize--;
    size_t lower=polyhedronArraySize(curSize);

    // Replacing a variable corresponds to copying the last row of the
    // triangular matrix to the row and column var. Since the
    // valuation of the source variable is always expressed in the
    // second variable, all polyhedra which wind up in the column have
    // to have their two variables exchanged. Start with copying the
    // row part.
    for (size_t i=0; i<var; i++) {
      polyhedron(i,var)=(additive ?
			 Polyhedron<isZ>(polyhedron(i,var),
				    bounds[i].interval,
				    bounds[var].interval,
				    polyhedra[lower],
				    bounds[i].interval,
				    bounds[curSize].interval) :
			 polyhedra[lower]);
      lower++;
    };
    // Skip the projection that talks about the relationship between
    // the two operands of this update operation.
    lower++;
    // Copy the remainder of the column into the row of the (var)
    // variable. Exchange x and y coordinates before copying.
    for (size_t i=var+1; i<curSize; i++) {
      polyhedra[lower].swapVars();
      polyhedron(var,i)=(additive ?
			 Polyhedron<isZ>(polyhedron(var,i),
				    bounds[var].interval,
				    bounds[i].interval,
				    polyhedra[lower],
				    bounds[i].interval,
				    bounds[curSize].interval) :
			 polyhedra[lower]);
      lower++;
    };
    // We should have taken care of all entries in the last row.
    assert(lower==upper);

    // Update the bound on the goal variable.
    if (additive)
      bounds[var].interval=Interval<isZ>(bounds[curSize].interval,
					 bounds[var].interval);
    else
      bounds[var].interval=bounds[curSize].interval;
    bounds[var].xRef=bounds[curSize].xRef;
    // Remove the last row from the set of variables.
    polyhedra.resize(polyhedronArraySize(curSize));

    // Remove the last interval from the bounds vector.
    bounds.resize(curSize);
  }
}

  template<bool isZ>
  DomVar DenseTvpi<isZ>::removeLast() {
    assert(bounds.size());
    TvpiVar last=bounds.size()-1;
    DomVar res = bounds[last].xRef;
    if (last>1) {
      polyhedra.resize(polyhedronArraySize(last));
    } else {
      polyhedra.clear();
    };
    bounds.resize(last);
    return res;
  }

// Assign the convex hull of all values from the target and the source
// to the target.
template<bool isZ>
void DenseTvpi<isZ>::augment(TvpiVar source, TvpiVar target)
  throw (IllegalArgument) {
  if ((target>=bounds.size()) || (source>=bounds.size()))
    throw IllegalArgument("DenseTvpi::augment",
			  "variable number out of range");
  if (target==source) return;

  size_t curSize=bounds.size();

  // Take the pair-wise convex hull between all i pairs of projections
  // (i,source) and (i,target) and store the result in (i,target).
  // Since each projections are stored such that the first index is
  // always strictly smaller than the second, there are three ranges:
  // (1) if both i<source and i<target, (2) either i>source or
  // i>target and (3) both i>source and i>target. The case i=source is
  // skipped, finally the case i=target updates the (source,target)
  // projection.

  // Calculate ranges.
  TvpiVar lower=min(target,source);
  TvpiVar higher=max(target,source);

  // Join the projections in range (1) where the running index is
  // lower than both source and target.

#ifdef DEBUG_AUGMENT
  cerr << "augment: first range 0.." << lower << " (="
       << (lower==source ? "source" : "target") << ")" << endl;
#endif // DEBUG_AUGMENT
  for (TvpiVar i=0; i<lower; i++)
    polyhedron(i,target)= Polyhedron<isZ>(polyhedron(i,target),
				     bounds[i].interval,
				     bounds[target].interval,
				     polyhedron(i,source),
				     bounds[i].interval,
				     bounds[source].interval);

  // Join the projections in range (2) where the variables of the
  // source polyhedron have to be swapped around before calculating
  // the convex hull with the target.

  // For the (source,target) projection, the source variables is in an
  // equality relation with itself. Create a polyhedron for this.
  Polyhedron<isZ> p;
  vector<Inequality*> ineq;
  ineq.push_back(new Inequality(-1,1,0));
  ineq.push_back(new Inequality(1,-1,0));
  p.addInequalitySet(bounds[source].interval,
		     bounds[source].interval,
		     ineq);
  assert( (bounds[source].interval.isSingleton() ?
	   p.getNoOfInequalities()==0 : true) );

  if (source<target) {
#ifdef DEBUG_AUGMENT
    cerr << "augment: source,target projection"  << endl;
#endif // DEBUG_AUGMENT
    // Handle the special projection (source,target).
    polyhedron(source,target)=Polyhedron<isZ>(polyhedron(source,target),
					 bounds[source].interval,
					 bounds[target].interval,
					 p,
					 bounds[source].interval,
					 bounds[source].interval);
    // Handle the range (2).
    for (TvpiVar i=source+1; i<target; i++) {
      Polyhedron<isZ> p = polyhedron(source,i);
      polyhedron(i,target)= Polyhedron<isZ>(polyhedron(i,target),
				       bounds[source].interval,
				       bounds[target].interval,
				       p.swapVars(),
				       bounds[target].interval,
				       bounds[source].interval);
    }
  } else { // target<source
#ifdef DEBUG_AUGMENT
    cerr << "augment: target,source projection"  << endl;
#endif // DEBUG_AUGMENT
    // Handle the special projection (target,source).
    polyhedron(target,source)=Polyhedron<isZ>(polyhedron(target,source),
					 bounds[target].interval,
					 bounds[source].interval,
					 p,
					 bounds[source].interval,
					 bounds[source].interval);
    for (TvpiVar i=target+1; i<source; i++) {
      Polyhedron<isZ> p = polyhedron(i,source);
      polyhedron(target,i)= Polyhedron<isZ>(polyhedron(target,i),
				       bounds[target].interval,
				       bounds[i].interval,
				       p.swapVars(),
				       bounds[source].interval,
				       bounds[i].interval);
    }
  };

#ifdef DEBUG_AUGMENT
  cerr << "augment: last range " << higher+1 << ".." << curSize-1
       << " (=" << (higher==source ? "source" : "target") << "+1 .. end)"
       << endl;
#endif // DEBUG_AUGMENT
  // Join the projections in range (3) where the running index is
  // higher than both source and target.
  for (TvpiVar i=higher+1; i<curSize; i++)
    polyhedron(target,i)= Polyhedron<isZ>(polyhedron(target,i),
				     bounds[target].interval,
				     bounds[i].interval,
				     polyhedron(source,i),
				     bounds[source].interval,
				     bounds[i].interval);

  // Extend the target interval by the source interval.
  bounds[target].interval = Interval<isZ>(bounds[target].interval,
					  bounds[source].interval);
}

// Set the target to the value of the source.
template<bool isZ>
void DenseTvpi<isZ>::update(TvpiVar source, TvpiVar target)
  throw (IllegalArgument) {
  if ((target>=bounds.size()) || (source>=bounds.size()))
    throw IllegalArgument("DenseTvpi::update",
			  "variable number out of range");
  if (target==source) return;

  size_t curSize=bounds.size();

  // Copy all i projections (i,source) to (i,target).  Since each
  // projection is stored such that the first index is always
  // strictly smaller than the second, there are three ranges: (1) if
  // both i<source and i<target, (2) either i>source or i>target and
  // (3) both i>source and i>target. The case i=source is skipped,
  // finally the case i=target updates the (source,target) projection.

  // Calculate ranges.
  TvpiVar lower=min(target,source);
  TvpiVar higher=max(target,source);

  // Copy the projections in range (1) where the running index is
  // lower than both source and target.
  for (TvpiVar i=0; i<lower; i++)
    polyhedron(i,target)= polyhedron(i,source);

  // Copy the projections in range (2) where the variables of the
  // source polyhedron have to be swapped around.

  // For the (source,target) projection, the two variables are in an
  // equality relation. Create a polyhedron for this.
  Polyhedron<isZ> p;
  vector<Inequality*> ineq;
  ineq.push_back(new Inequality(-1,1,0));
  ineq.push_back(new Inequality(1,-1,0));
  p.addInequalitySet(bounds[source].interval,
		     bounds[source].interval,
		     ineq);
  assert( (bounds[source].interval.isSingleton() ?
	   p.getNoOfInequalities()==0 : true) );

  if (source<target) {
    // Handle the special projection (source,target).  
    polyhedron(source,target)=p;

    // Handle the range (2).
    for (TvpiVar i=source+1; i<target; i++) {
      Polyhedron<isZ> p = polyhedron(source,i);
      polyhedron(i,target) = p.swapVars();
    }
  } else { // target<source
    polyhedron(target,source)=p;

    for (TvpiVar i=target+1; i<source; i++) {
      Polyhedron<isZ> p = polyhedron(i,source);
      polyhedron(target,i) = p.swapVars();
    }
  };

  // Join the projections in range (3) where the running index is
  // higher than both source and target.
  for (TvpiVar i=higher+1; i<curSize; i++)
    polyhedron(target,i) = polyhedron(source, i);

  // Set the target interval to the source interval.
  bounds[target].interval = bounds[source].interval;
}


// Insert a new variable into the polyhedron. There will be no bound
// on the created variable.
template<bool isZ>
TvpiVar DenseTvpi<isZ>::createVariable() {
  TvpiVar varId=bounds.size();
#ifdef DEBUG_NEWVAR
  cerr << "Tvip.createVariable: created unbounded var p_" << varId << endl;
#endif // DEBUG_NEWVAR
  size_t newSize=varId+1;
  bounds.resize(newSize);
  // Allocate new empty polyhedra if there are at least 2 dimensions.
  if (newSize>1) polyhedra.resize(polyhedronArraySize(newSize));
  return varId;
}

// Insert a new variable into the polyhedron. The specified argument
// contains the initial value of the variable.
template<bool isZ>
TvpiVar DenseTvpi<isZ>::createVariable(long value) {
  TvpiVar varId = bounds.size();
#ifdef DEBUG_NEWVAR
  cerr << "Tvip.createVariable: created new var p_" << varId
       << " = " << value << endl;
#endif // DEBUG_NEWVAR
  size_t newSize=varId+1;
  bounds.push_back(PolyBound(Interval<isZ>(value)));
  // Allocate new empty polyhedra if there are at least 2 dimensions.
  if (newSize>1) polyhedra.resize(polyhedronArraySize(newSize));
  return varId;
}

// Insert a new variable into the polyhedron. The specified argument
// contains the initial value of the variable.
template<bool isZ>
TvpiVar DenseTvpi<isZ>::createVariable(mpz_class value, DomVar xRef) {
  TvpiVar varId = bounds.size();
#ifdef DEBUG_NEWVAR
  cerr << "Tvip.createVariable: created new var p_" << varId
       << " = " << value << " corresp. to x_" << xRef << endl;
#endif // DEBUG_NEWVAR
  size_t newSize=varId+1;
  bounds.push_back(PolyBound(Interval<isZ>(value), xRef));
  // Allocate new empty polyhedra if there are at least 2 dimensions.
  if (newSize>1) polyhedra.resize(polyhedronArraySize(newSize));
  return varId;
}

// Insert a new variable into the polyhedron. The specified argument
// contains the initial interval of the variable.
template<bool isZ>
TvpiVar DenseTvpi<isZ>::createVariable(Interval<isZ>& value, DomVar xRef) {
  TvpiVar varId = bounds.size();
#ifdef DEBUG_NEWVAR
  cerr << "Tvip.createVariable: created new var p_" << varId
       << " = " << value << " corresp. to x_" << xRef << endl;
#endif // DEBUG_NEWVAR
  size_t newSize=varId+1;
  bounds.push_back(PolyBound());
  bounds.back().interval.swap(value);
  bounds.back().xRef=xRef;
  // Allocate new empty polyhedra if there are at least 2 dimensions.
  if (newSize>1) polyhedra.resize(polyhedronArraySize(newSize));
  return varId;
}

  // Ask for the xRefs of all relational variables. The first
  // passed-in array must have n elements, the second one (n*n-n)/2
  // elements where n is the number returned by relVarsInDomain.
  template<bool isZ>
  void DenseTvpi<isZ>::queryXRefs(DomVar* xRefPtr, int* relPtr) const {
    assert(xRefPtr);
    assert(relPtr);
    for (typename vector<PolyBound>::const_iterator iter=bounds.begin();
	 iter!=bounds.end(); iter++) *xRefPtr++=iter->xRef;
    for (typename Polyhedra::const_iterator iter=polyhedra.begin();
	 iter!=polyhedra.end(); iter++)
      *relPtr++=iter->getNoOfInequalities();
  };

// Propagate the changed bounds of the given variable to the other
// bounds. Returns false if the system turned unsatisfiable.
template<bool isZ>
bool DenseTvpi<isZ>::propagateBounds(TvpiVar var) {
  size_t curSize=bounds.size();
  assert(var<curSize);
  assert(bounds[var].interval.hasChanged());

  if (bounds[var].interval.isEmpty()) return false;

  // Propagate the changed bound to all other bounds.
  for (size_t i=0; i<curSize; i++) {
    if (i==var) continue;
    if (i<var)
      polyhedron(i,var).propagateYBounds(bounds[i].interval,
					 bounds[var].interval);
    else
      polyhedron(var,i).propagateXBounds(bounds[var].interval,
					 bounds[i].interval);
    if (bounds[i].interval.isEmpty()) return false;
  };
  bounds[var].interval.clearChanged();
  return true;
}

template<bool isZ>
vector<Inequality*> DenseTvpi<isZ>::jResultants;
template<bool isZ>
vector<Inequality*> DenseTvpi<isZ>::kResultants;
template<bool isZ>
vector<Inequality*> DenseTvpi<isZ>::ineqVia;

// Insert a set of inequalities over two specified variables and close
// the system. The inequalities must be sorted by angle if sorted is
// true. The return value is false if the system became unsatisfiable.
template<bool isZ>
Result DenseTvpi<isZ>::addInequalitySet(TvpiVar j,
					TvpiVar k,
					vector<Inequality*>& ineqNew,
					bool sorted)
  throw (IllegalArgument) {

  size_t sizeNew=ineqNew.size();
  if (sizeNew==0) return resRedundant;

  size_t curSize=bounds.size();

  if (j>=curSize || k>=curSize)
    throw IllegalArgument("DenseTvpi::addInequalities",
			  "variable number out of range");

  if (j==k)
    throw IllegalArgument("DenseTvpi::addInequalities",
			  "cannot add inequality over the same variable");

  // Ensure that the j is smaller than k.
  if (j>k) {
    for (size_t i=0; i<sizeNew; i++)
      ineqNew[i]->swapVars();
    TvpiVar temp=j;
    j=k;
    k=temp;
    sorted=false;
  }

  // In case the user didn't bother sorting the inequalities, sort
  // them now.
  if (!sorted) sortAndRemoveQuasiSyntacticRed(0, ineqNew);

  // Reset the changed flag of the intervals.
  for(size_t i=0; i<bounds.size(); i++) bounds[i].interval.clearChanged();

  // Insert the new inequalities in the polyhedron that contains the
  // variables j and k. The smaller variable is always found on the
  // x-axis.
  size_t startNew=0;
  polyhedron(j,k).addInequalitySet(bounds[j].interval, bounds[k].interval,
				   ineqNew, startNew);
  sizeNew=ineqNew.size()-startNew;

  Partitioning part;
  if (sizeNew) part = Partitioning(sizeNew, &ineqNew[startNew]);

#ifdef DEBUG_CLOSURE
  cerr << "Adding " << sizeNew << " new non-redundant inequalities:" << endl;
  for(size_t i=startNew; i<ineqNew.size(); i++) cerr << *ineqNew[i] << endl;
#endif // DEBUG_CLOSURE

  // Check if the polyhedron is still satisfiable. Return if this is
  // not the case.
  if (bounds[j].interval.isEmpty() || bounds[k].interval.isEmpty())
    return resUnsatisfiable; else {

    // Return here if none of the new inequalities changed anything.
    if (sizeNew==0) {
      Result res=resRedundant;
      // There might be no non-redundant inequalities at this point,
      // but the bounds might have changed.
      if (bounds[j].interval.hasChanged()) {
	res=resChanged;
	if (!propagateBounds(j)) return resUnsatisfiable;
      };
      if (bounds[k].interval.hasChanged()) {
	res=resChanged;
	if (!propagateBounds(k)) return resUnsatisfiable;
      };
      return res;
    };

    // Update all projections pairs with distance one, i.e. those that
    // share one variable with {j,k}. The resultants of combining i,k
    // and k,j and the resultants of combining i,j and j,k are stored
    // in two arrays. The inequalities are stored consecutive with a
    // second array indexShareX keeping track on where the set of
    // inequalities for the next projection starts. After calculating
    // the resultants the inequalities are added to the projections
    // i,j and i,k and redundant elements are removed. The resulting
    // two arrays of inequality sets are later used to update the rest
    // of the triangle.

    // Allocate arrays to safe the partitioning and size of the
    // resultants of the ith projection. Create an array that is one
    // element larger. The last element will keep the size of the
    // whole set of resultants.
    size_t indexShareJ[curSize+1];
    size_t indexShareK[curSize+1];
    
    Partitioning partShareJ[curSize];
    Partitioning partShareK[curSize];

    // Clear the static lists of distance-one resultants.
    jResultants.clear();
    kResultants.clear();

    // Calculate the resultants that share one variable with the
    // updated polyhedron.
    size_t nonRedResultants=0;
    for(size_t i=0; i<curSize; i++) {
      // Remember the index at which the inequalities for i start.
      indexShareJ[i]=jResultants.size();
      indexShareK[i]=kResultants.size();

      // Don't calculate resultants with the projection j,k itself.
      if (i==j) continue;
      if (i==k) continue;
	
#ifdef DEBUG_CLOSURE
      cerr << "Calculating resultants for (i,j) and (j,k) where i=x_"
	   << bounds[i].xRef << ", j=x_" << bounds[j].xRef
	   << " and k=x_" << bounds[k].xRef << endl;
#endif // DEBUG_CLOSURE

      // Propagate the bound from j to i and remove redundant
      // inequalities in the projection to minimize the number of
      // resultants.
      if (bounds[j].interval.hasChanged()) if (i<j) {
	polyhedron(i,j).propagateYBounds(bounds[i].interval,
					 bounds[j].interval); 
	polyhedron(i,j).enforceBounds(bounds[i].interval,
				      bounds[j].interval);
      } else {
	polyhedron(j,i).propagateXBounds(bounds[j].interval,
					 bounds[i].interval);
	polyhedron(j,i).enforceBounds(bounds[j].interval,
				      bounds[i].interval);
      };
    
      // Calculate the resultants by combining the inequalities of the
      // projections i,j and j,k to get a set of constraints over i,k.
      (i<j ? polyhedron(i,j) : polyhedron(j,i)).
	calcResultants(j>i,
		       bounds[i].interval,
		       &ineqNew[startNew],
		       part,
		       false, // j>k
		       bounds[k].interval,
		       jResultants,
		       i>k);
	
#ifdef DEBUG_CLOSURE
      for(size_t e=indexShareJ[i]; e<jResultants.size(); e++) {
	jResultants[e]->output(cerr, bounds[i<k ? i : k].xRef,
			       bounds[i<k ? k : i].xRef);
	cerr << endl;
      };
	  
      cerr << "Calculating resultants for (i,k) and (k,j) where i=x_"
	   << bounds[i].xRef << ", j=x_" << bounds[j].xRef
	   << " and k=x_" << bounds[k].xRef << endl;
#endif // DEBUG_CLOSURE
      
      // Propagate the bound from k to i and remove redundancies.
      if (bounds[j].interval.hasChanged()) if (i<k) {
	polyhedron(i,k).propagateYBounds(bounds[i].interval,
					 bounds[k].interval); 
	polyhedron(i,k).enforceBounds(bounds[i].interval,
				      bounds[k].interval);
      } else {
	polyhedron(k,i).propagateXBounds(bounds[k].interval,
					 bounds[i].interval);
	polyhedron(k,i).enforceBounds(bounds[k].interval,
				      bounds[i].interval);
      };
      
      // Calculate the resultants by combining the inequalities of the
      // projections i,k and k,j to get a set of constraints over i,j.
      (i<k ? polyhedron(i,k) : polyhedron(k,i)).
	calcResultants(k>i,
		       bounds[i].interval,
		       &ineqNew[startNew],
		       part,
		       true, // k>j
		       bounds[j].interval,
		       kResultants,
		       i>j);
#ifdef DEBUG_CLOSURE
      for(size_t e=indexShareK[i]; e<kResultants.size(); e++) {
	kResultants[e]->output(cerr, bounds[i<j ? i : j].xRef,
			       bounds[i<j ? j : i].xRef);
	cerr << endl;
      };
#endif // DEBUG_CLOSURE
	
      // Add the resultants to the i,k projection which is then up-to-date.
      (i<k ? polyhedron(i,k) : polyhedron(k,i)).
	addInequalitySet(bounds[min(i,k)].interval,
			 bounds[max(i,k)].interval,
			 jResultants,
			 indexShareJ[i]);

      // Count the number of non-redundant inequalities.
      size_t jSize=jResultants.size()-indexShareJ[i];
      nonRedResultants+=jSize;
      
      // Update the partitioning information for this projection.
      if (jSize)
	partShareJ[i]=Partitioning(jSize, &jResultants[indexShareJ[i]]);

      // Dto. for i,j.
      (i<j ? polyhedron(i,j) : polyhedron(j,i)).
	addInequalitySet(bounds[min(i,j)].interval,
			 bounds[max(i,j)].interval,
			 kResultants,
			 indexShareK[i]);

      size_t kSize=kResultants.size()-indexShareK[i];
      nonRedResultants+=kSize;
      if (kSize) 
	partShareK[i]=Partitioning(kSize, &kResultants[indexShareK[i]]);

#ifdef DEBUG_CLOSURE
      cerr << "Propagating " << jSize << " and " << kSize
	   << " non-redundant inequalities." << endl;
      if (nonRedResultants) {
	for(size_t e=indexShareJ[i]; e<jResultants.size(); e++) {
	  jResultants[e]->output(cerr, bounds[i<k ? i : k].xRef,
				 bounds[i<k ? k : i].xRef);
	  cerr << endl;
	};
	for(size_t e=indexShareK[i]; e<kResultants.size(); e++) {
	  kResultants[e]->output(cerr, bounds[i<j ? i : j].xRef,
				 bounds[i<j ? j : i].xRef);
	  cerr << endl;
	};
      };
#endif // DEBUG_CLOSURE
    };
      
    // The new inequalities are not used beyond this point.
    if (bounds[j].interval.isEmpty() || bounds[k].interval.isEmpty())
      return resUnsatisfiable;

    // Remember the size of the whole array.
    indexShareJ[curSize]=jResultants.size();
    indexShareK[curSize]=kResultants.size();

    // If updating the projections that share one variable
    // ("distance-one" projections) with j,k changed the system
    // (nonRedResultants>0), the remaining projections might change
    // as well. In this case create partitioning information for all
    // distance-one projections and calculate and insert resultants
    // for all distance-two projections.
    if (nonRedResultants) {

#ifdef DEBUG_CLOSURE
      cerr << "Calculate shortest paths with one via." << endl;
#endif // DEBUG_CLOSURE

      // Now close the remaining projections p_{x,y} by analyzing if
      // there is a shorter path (tighter constraints) from x to y
      // by either going via j (p_{x,j} \join p_{j,y}) or via k
      // (p_{x,k} \join p_{k,y}).
      for(size_t x=0; x<curSize; x++) {
	if (x==j) continue;
	if (x==k) continue;

	// No new inequalities for the xth distance-one projection
	// means no possible resultants, so skip this index.
	if (indexShareJ[x+1]-indexShareJ[x]==0 &&
	    indexShareK[x+1]-indexShareK[x]==0) continue;

#ifdef DEBUG_CLOSURE
	cerr << "Adding resultants to column p_" << x 
	     << ", i.e. x_" << bounds[x].xRef << endl;
#endif // DEBUG_CLOSURE

	for(size_t y=x+1; y<curSize; y++) {
	  if (y==j) continue;
	  if (y==k) continue;

	  // Store the pair-wise resultants in ineqVia.
	  ineqVia.clear();
	
	  // The set jResultants[indexShareJ[x]] expresses
	  // inequalities over x and k (j was elimintated). Insert the
	  // resultants arising by elimintating k in x,k and k,y at
	  // the beginning of ineqVia.
#ifdef DEBUG_CLOSURE
	  cerr << "Calculating resultants for (x,k) and (k,y) where x=x_"
	       << bounds[x].xRef << ", k=x_" << bounds[k].xRef
	       << " and y=x_" << bounds[y].xRef << endl;
#endif // DEBUG_CLOSURE
	  if (partShareJ[x][total] && partShareJ[y][total])
	    calcResultantsOfSets(&jResultants[indexShareJ[x]],
				 partShareJ[x],
				 x<k,
				 bounds[x].interval,
				 &jResultants[indexShareJ[y]],
				 partShareJ[y],
				 y<k,
				 bounds[y].interval,
				 ineqVia);
	  // The set kResultants[indexShareK[x]] expresses
	  // inequalities over x and j (k was elimintated). Insert the
	  // resultants arising by elimintating j in x,j and j,y at
	  // the first free entry in ineqVia.
#ifdef DEBUG_CLOSURE
	  cerr << "Calculating resultants for (x,j) and (j,y) where x=x_"
	       << bounds[x].xRef << ", j=x_" << bounds[j].xRef
	       << " and y=x_" << bounds[y].xRef << endl;
#endif // DEBUG_CLOSURE
	  if (partShareK[x][total] && partShareK[y][total])
	    calcResultantsOfSets(&kResultants[indexShareK[x]],
				 partShareK[x],
				 x<j,
				 bounds[x].interval,
				 &kResultants[indexShareK[y]],
				 partShareK[y],
				 y<j,
				 bounds[y].interval,
				 ineqVia);

	  // Sort the set of inequalities and insert them into the
	  // appropriate projection.
	  sortAndRemoveQuasiSyntacticRed(0, ineqVia);

#ifdef DEBUG_CLOSURE
	  if (!ineqVia.empty()) {
	    cerr << "Created " << ineqVia.size() 
		 << " resultant(s) with row p_" << y
		 << ", i.e. x_" << bounds[y].xRef << ":" << endl;
	    for (size_t i=0; i<ineqVia.size(); i++) {
	      ineqVia[i]->output(cerr, bounds[x].xRef, bounds[y].xRef);
	      cerr << endl;
	    }
	  };
#endif // DEBUG_CLOSURE

	  polyhedron(x,y).addInequalitySet(bounds[x].interval,
					   bounds[y].interval,
					   ineqVia);
		
	  // This distance-two projection is completed. If any
	  // inequalities were added, check for satisfiability.
	  if (bounds[y].interval.isEmpty()) return resUnsatisfiable;
	}
	// The xth projections j,x and x,k are no longer used.
	if (bounds[x].interval.isEmpty()) return resUnsatisfiable;
      }
    }
  };
  // Each time inequalities are inserted into a particular projection
  // the satisfiablity is checked and the function returns prematurely
  // if the system is unsatisfiable. If the system has stabilized and is
  // still feasible, this block is executed.

  return resChanged;
};
  
  

  // Create TVPI inequalities that approximate the general inequality.
  template<bool isZ>
  Result DenseTvpi<isZ>::
  approximateInequality(const vector<LinComponent>& comps,
			mpz_class constant,
			bool isEquality) {
    // Require at least one variable.
    assert(comps.size()>0);
    
#ifdef DEBUG_APPROXIMATE
    cerr << "approximating:";
    for(std::vector<LinComponent>::const_iterator termIter=comps.begin();
	termIter!=comps.end(); termIter++) {
      cerr << " ";
      if (termIter!=comps.begin() &&
	  termIter->coefficient>=0) cerr << "+";
      cerr << termIter->coefficient << " p_" << termIter->variable << "(=x_"
	   << xRef(termIter->variable) << ")";
    };
    if (constant>=0) cerr << "+";
    cerr << constant;
    if (isEquality) cerr << " = 0"; else cerr << " <= 0";
    cerr << endl;
#endif // DEBUG_APPROXIMATE

    // Clear the flags of the bounds so that a redundant bound does not
    // return 'Changed'.
    for(vector<LinComponent>::const_iterator x=comps.begin();
	x!=comps.end(); x++) bounds[x->variable].interval.clearChanged();

    vector<Inequality*> eqs;
    eqs.reserve(2);

    bool allRedundant = true;

    vector<TVPIComponent> otherTerms;
    otherTerms.reserve(comps.size()-1);

    // Now create inequalities for all pairs of variables in the
    // linear expression.
    for(vector<LinComponent>::const_iterator x=comps.begin();
	x!=comps.end(); x++) {
      for(vector<LinComponent>::const_iterator y=x+1; y!=comps.end(); y++) {
      
#ifdef DEBUG_APPROXIMATE
	cerr << "infering bounds for "
	     << x->coefficient << " p_" << x->variable << " + "
	     << y->coefficient << " p_" << y->variable << endl;
#endif // DEBUG_APPROXIMATE

	for(vector<LinComponent>::const_iterator i=comps.begin();
	    i!=comps.end(); i++) {
	  if (i==x) continue;
	  if (i==y) continue;
	  otherTerms.push_back(*i);
	};

	mpz_class termValue;
        bool isFinite;
	
	// We separated x and y to get ts + cx x + cy y + c = 0. The
	// strongest information we can gleam from this equation
	// (inequality) is [min(ts), max(ts)] + cx x + cy y + c = 0. 
	// Bringing the interval to the other side negates it and thus
	// cx x + cy y +c = [-max(ts) , -min(ts)]. Hence, it follows
	// that -max(ts) <= cx x + cy y + c <= -min(ts). Calculate
	// max(ts) and min(ts)=-max(-ts) where -ts is ts with every
	// coefficient negated.
	if (isEquality) {
	  // Calculate max(ts).
	  isFinite = linOpt(otherTerms, termValue);
	
#ifdef DEBUG_APPROXIMATE
	  cerr << "equality: maximum of ";
	  for(vector<TVPIComponent>::iterator i=otherTerms.begin();
	      i!=otherTerms.end(); i++)
	    cerr << "+" << i->coefficient << " p_" << i->variable;
	  if (isFinite) cerr << " is " << termValue; 
	  else cerr << " does not exist";
	  cerr << endl;
#endif // DEBUG_APPROXIMATE

	  // Intersect with -max(ts) <= cx x + cy y + c which is
	  // equivalent to -cx x -cy y <= max(ts)+c
	  if (isFinite) {
	    eqs.push_back(new Inequality(-x->coefficient, -y->coefficient,
					 termValue+constant));
	    eqs.back()->canonicalize(isIntegral);
	  }
	};

	// Negate the coefficients in ts in order to calculate
	// min(ts) = -max(-ts).
	for(vector<TVPIComponent>::iterator i=otherTerms.begin();
	    i!=otherTerms.end(); i++) i->coefficient=-i->coefficient;

	// Calculate -min(ts).
	isFinite = linOpt(otherTerms, termValue);

#ifdef DEBUG_APPROXIMATE
	cerr << "maximum value of ";
	for(vector<TVPIComponent>::iterator i=otherTerms.begin();
	    i!=otherTerms.end(); i++)
	  cerr << "+" << i->coefficient << " p_" << i->variable;
	if (isFinite) cerr << " is " << termValue; 
	else cerr << " does not exist";
	cerr << endl;
#endif // DEBUG_APPROXIMATE

	// Intersect with cx x + cy y + c <= -min(ts) which is
	// equivalent to cx x + cy y <= -min(ts)-c.
	if (isFinite) {
	  eqs.push_back(new Inequality(x->coefficient, y->coefficient,
				       termValue-constant));
	  eqs.back()->canonicalize(isIntegral);
	};

#ifdef DEBUG_APPROXIMATE
	  cerr << "adding ";
	  for(vector<Inequality*>::const_iterator i=eqs.begin();
	      i!=eqs.end(); i++)
	    (*i)->output(cerr,x->variable, y->variable) << ", ";
#endif // DEBUG_APPROXIMATE

	Result res = addInequalitySet(x->variable, y->variable, eqs);

#ifdef DEBUG_APPROXIMATE
	cerr << ": result " << res << endl;
#endif // DEBUG_APPROXIMATE

	eqs.clear();
	otherTerms.clear();

	if (res==resUnsatisfiable) return res;
	if (res==resChanged) allRedundant=false;
      }
    };

    // Now create bounds for each variable in the linear expression.
    for(vector<LinComponent>::const_iterator x=comps.begin();
	x!=comps.end(); x++) {
      
#ifdef DEBUG_APPROXIMATE
      cerr << "infering bound for "
	   << x->coefficient << " p_" << x->variable << endl;
#endif // DEBUG_APPROXIMATE

      for(vector<LinComponent>::const_iterator i=comps.begin();
	  i!=comps.end(); i++) {
	if (i==x) continue;
	otherTerms.push_back(*i);
      };

      mpz_class termValue;
      bool isFinite;
	
      // We separated out x to get ts + cx x + c = 0. The strongest
      // information we can gleam from this equation (inequality) is
      // [min(ts), max(ts)] + cx x + c = 0. Bringing the interval to
      // the other side gives cx x + c = [-max(ts), -min(ts)], that is,
      // -max(ts) <= cx x + c <= max(-ts) (see above).
      if (isEquality) {
	// Calculate max(ts).
	isFinite = linOpt(otherTerms, termValue);
	
#ifdef DEBUG_APPROXIMATE
	cerr << "equality: maximum of ";
	for(vector<TVPIComponent>::iterator i=otherTerms.begin();
	    i!=otherTerms.end(); i++)
	  cerr << "+" << i->coefficient << " p_" << i->variable;
	if (isFinite) cerr << " is " << termValue; 
	else cerr << " does not exist";
	cerr << endl;
#endif // DEBUG_APPROXIMATE
	
	// Given max(ts), an upper bound on x can be inferred by using
	// cx x + c <= [-max(ts), -min(ts)].
	if (isFinite) {
	  // If cx is positive, then x >= (-max(ts)-c)/cx.
	  if (sgn(x->coefficient)>=0) {
	    assert(sgn(x->coefficient)!=0);
	    mpq_class newBound = mpq_class((-termValue)-constant,
					   x->coefficient);
	    newBound.canonicalize();
#ifdef DEBUG_APPROXIMATE
	    cerr << newBound << " <= p_" << x->variable << endl;
#endif // DEBUG_APPROXIMATE
	    bounds[x->variable].interval.updateLower(newBound);
	  } else {
	    // Since cx is negative, -x >= (-max(ts)-c)/-cx, hence
	    // x <= max(ts)+c/-cx
	    mpq_class newBound = mpq_class(termValue+constant,
					   -x->coefficient);
	    newBound.canonicalize();
#ifdef DEBUG_APPROXIMATE
	    cerr << "p_" << x->variable << " <= " << newBound << endl;
#endif // DEBUG_APPROXIMATE
	    bounds[x->variable].interval.updateUpper(newBound);
	  }
	}
      };
      
      // Negate the coefficients in ts in order to calculate
      // min(ts) = -max(-ts).
      for(vector<TVPIComponent>::iterator i=otherTerms.begin();
	  i!=otherTerms.end(); i++) i->coefficient=-i->coefficient;
      
      // Calculate -min(ts).
      isFinite = linOpt(otherTerms, termValue);
      
#ifdef DEBUG_APPROXIMATE
      cerr << "maximum value of ";
      for(vector<TVPIComponent>::iterator i=otherTerms.begin();
	  i!=otherTerms.end(); i++)
	cerr << "+" << i->coefficient << " p_" << i->variable;
      if (isFinite) cerr << " is " << termValue; 
      else cerr << " does not exist";
      cerr << endl;
#endif // DEBUG_APPROXIMATE
      
      // Given -min(ts), we can infer a bound from the
      // constraint cx x + c = [-max(ts) , -min(ts)] which is
      // equivalent to cx x <= -min(ts)-c.
      if (isFinite) {
	// If cx is positive, then x <= -min(ts)-c/cx.
	if (sgn(x->coefficient)>=0) {
	  assert(sgn(x->coefficient)!=0);
	  mpq_class newBound = mpq_class(termValue-constant,
					 x->coefficient);
	  newBound.canonicalize();
#ifdef DEBUG_APPROXIMATE
	  cerr << "p_" << x->variable << " <= " << newBound << endl;
#endif // DEBUG_APPROXIMATE
	  bounds[x->variable].interval.updateUpper(newBound);
	} else {
	  // Since cx is negative, enforce -x <= -min(ts)-c/-cx
	  // which is equivalent to x >= min(ts)+c/-cx.
	  mpq_class newBound = mpq_class(constant-termValue,
					 -x->coefficient);
	  newBound.canonicalize();
#ifdef DEBUG_APPROXIMATE
	  cerr << newBound << " <= p_" << x->variable << endl;
#endif // DEBUG_APPROXIMATE
	  bounds[x->variable].interval.updateLower(newBound);
	}
      };

#ifdef DEBUG_APPROXIMATE
      bool isRedundant = true;
#endif // DEBUG_APPROXIMATE
      
      bool unSat = bounds[x->variable].interval.isEmpty();
      if (!unSat && bounds[x->variable].interval.hasChanged()) {
#ifdef DEBUG_APPROXIMATE
	isRedundant = false;
#endif // DEBUG_APPROXIMATE
	allRedundant = false;
	if (!propagateBounds(x->variable)) unSat=true;
      };

#ifdef DEBUG_APPROXIMATE
      cerr << "result " << (unSat ? "unsatisfiable" :
			    (isRedundant ? "redundant" : "changed"))
	   << endl;
#endif // DEBUG_APPROXIMATE

      if (unSat) return resUnsatisfiable;
    
      eqs.clear();
      otherTerms.clear();
    
    }

    return (allRedundant ? resRedundant : resChanged);
  }

  // Multiply the value of the given variable by 2^m.
  template<bool isZ>
  void DenseTvpi<isZ>::stretch(TvpiVar var, Mult m) {
    assert(var<bounds.size());
    if (m==0) return;
    bounds[var].interval=bounds[var].interval.stretch(m);
    for (TvpiVar i=0; i<var; i++)
      polyhedron(i,var)=polyhedron(i,var).stretch(0,m);
    for (TvpiVar i=var+1; i<bounds.size(); i++)
      polyhedron(var,i)=polyhedron(var,i).stretch(m,0);
  }

  // Calculate the minimal multiplicity of the given variable, using
  // the equality relationships in the domain.
  
  // This function is unfinished.
  template<bool isZ>
  Mult DenseTvpi<isZ>::seekMultFromEquality(TvpiVar var, Mult* mults) {
    assert(var<bounds.size());
    Mult m = mults[var];
    for (TvpiVar x=0; x<bounds.size(); x++) {
      if (x==var) continue;
      if (x<var) {
	size_t idx = polyhedronIndex(x,var);
	if (polyhedra[idx].isEquality()) {
	  // A polyhedron consisting of two coinciding inequalities must have
	  // the first inequality in the east or north quadrant and, hence,
	  // the coefficient in front of the y variable must be positive.
	}
      } else { //y<x
      }
      
    }
    assert(false);
    return 0;
  }

  // Substitute the variable y for another if they are in an equality
  // relationship.
  template<bool isZ>
  TvpiVar DenseTvpi<isZ>::substEquality(TvpiVar y,
					mpz_class& k,
					mpz_class& c,
					mpz_class& m) const {
    assert(sgn(k)!=0);
    for (TvpiVar x=0; x<y; x++) {
#ifdef DEBUG_SUBSTITUTE
      cerr << "checking if p_" << y << " is equal to p_" << x
	   << ": " << ( polyhedron(x,y).isEquality() ? "yes" : "no" )
	   << endl;
#endif // DEBUG_SUBSTITUTE
      size_t idx = polyhedronIndex(x,y);
      if (polyhedra[idx].isEquality()) {
	// A polyhedron consisting of two coinciding inequalities must have
	// the first inequality in the east or north quadrant and, hence,
	// the coefficient in front of the y variable must be positive.
	assert(sgn(polyhedra[idx][0]->getB())>=0);
	if (mpz_divisible_p(k.get_mpz_t(), 
			    polyhedra[idx][0]->getB().get_mpz_t())) {
	  mpz_class divisor;
	  mpz_divexact(divisor.get_mpz_t(), 
		       k.get_mpz_t(), polyhedra[idx][0]->getB().get_mpz_t());
	  k = -polyhedra[idx][0]->getA()*divisor;
	  c += polyhedra[idx][0]->getC()*divisor;
	  m = 1;
#ifdef DEBUG_SUBSTITUTE
	  cerr << "divisor is " << divisor
	       << ", adjusting constant by " << polyhedra[idx][0]->getC()
	       << " * divisor to " << c << endl;
#endif // DEBUG_SUBSTITUTE
	} else {
  	  m =  polyhedra[idx][0]->getB();
	  k = -polyhedra[idx][0]->getA();
	  c = (c*m)+polyhedra[idx][0]->getC();
#ifdef DEBUG_SUBSTITUTE
	  cerr << " multiplier is " << m 
	       << ", constant changed to c*m + " << polyhedra[idx][0]->getC()
	       << endl;
#endif // DEBUG_SUBSTITUTE
        };
	return x;
      }
    };
    return y;
  }

// Retrieve a projection from the TVPI system (for visualization).
template<bool isZ>
const Polyhedron<isZ>& DenseTvpi<isZ>::getProjection(TvpiVar x,
						     TvpiVar y) const {
  assert(x<bounds.size());
  assert(y<bounds.size());
  assert(x<y);
  return polyhedron(x,y);
}

template<bool isZ>
void DenseTvpi<isZ>::output(std::ostream& stream, Mult *mults) const {
  assert((bounds.size()<=1 && polyhedra.size()==0) ||
	 polyhedronArraySize(bounds.size())==polyhedra.size());

  bool showXRef = false;
  for(size_t i=0; i<bounds.size(); i++) if (bounds[i].xRef!=0) showXRef = true;

  for(TvpiVar y=0; y<bounds.size(); y++) {
    if (showXRef)
      if (bounds[y].xRef<0) 
	stream << "p" << y << " t_" << (-bounds[y].xRef); else
	stream << "p" << y << " x_" << bounds[y].xRef;
    else
      stream << "p" << y;
    stream << " = " << bounds[y].interval;
    if (mults)
      if (mults[y]==invalidMult) stream << ", mult n/a";
	else stream << ", mult " << (int) mults[y];
    // Domain::getVariableTop / getVariableType omitted (not available without affine.hh)
    stream << endl;
    for(TvpiVar x=0; x<y; x++) {
      size_t projection=y*(y-1)/2+x;
      const Polyhedron<isZ> p=polyhedra[projection];
      if (showXRef)
	Tvpi::output<isZ>(stream, p, bounds[x].xRef, bounds[y].xRef);
      else
	Tvpi::output<isZ>(stream, p, x, y);
    }
  }
};

template<bool isZ>
std::ostream& DenseTvpi<isZ>::showDist(std::ostream& stream) const {
  size_t curSize=bounds.size();


  bool showXRef = false;
  for(size_t i=0; i<curSize; i++) if (bounds[i].xRef!=0) showXRef = true;

  // Draw horizontal legend.
  stream.width(3);
  stream << bounds.size();
  stream << "  |";
  for(size_t x=0; x<curSize; x++) {
    stream << "x_";
    stream.setf(stream.left);
    stream.width(3);
    if (showXRef) stream << bounds[x].xRef; else stream << x;
    stream.unsetf(stream.left);
    stream << "|";
  };
  stream << endl;
  // Draw table.
  for(size_t y=0; y<curSize; y++) {
    stream << "x_";
    stream.setf(stream.left);
    stream.width(3);
    if (showXRef) stream << bounds[y].xRef; else stream << y;
    stream.unsetf(stream.left);
    stream << "|";
    for(size_t x=0; x<y; x++) {
      stream.width(2);
      stream << polyhedron(x,y).getNoOfInequalities()
	     << "*";
      stream.width(2);
      stream << polyhedron(x,y).isShared() << "|";
    };
    stream << endl;
  }
  return stream;
}

  template<bool isZ>
  void DenseTvpi<isZ>::sane() const {
    size_t noOfVars = bounds.size();
    if (noOfVars<=1) return;
    assert(polyhedra.size()==polyhedronArraySize(noOfVars));
    for (size_t i=0; i<polyhedra.size(); i++)
      polyhedra[i].sane();
  };

template class DenseTvpi<false>;
template class DenseTvpi<true>;

}; // namespace Tvpi
