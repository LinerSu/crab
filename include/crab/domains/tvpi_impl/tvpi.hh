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
// tvpi.h - A two-variables-per-inequality polyhedron.
//

#ifndef __TVPI_H
#define __TVPI_H

#include "common.hh"
#include "lincomponent.hh"
#include "polyhedron.hh"
#include "planar.hh"
#include <vector>



// TODO

// projectOnto: The function could be faster if the given variable
// sequence starts with 0,1,2,..., i.e. the first dimensions remain
// unchanged. Alternatively, it could return a new Tvpi system, since
// that comes free if the whole matrix is copied anyway.

// addInequalitySet: This function is only ever called with a fixed
// sized array. Check whether we can change the costly vector back to
// an array. This might not be possible with the shrinking around the
// integral hull algorithm.

namespace Tvpi {

const static TvpiVar invalidTvpiVar = (TvpiVar) -1;

// Declare a function type for the callback parameter in
// 'approximateInequality'.
typedef bool (*TvpiResultFunc) (TvpiVar,
				TvpiVar,
				size_t,
				Inequality*[]);

template<bool isZ>
class DenseTvpi {

  // A sorted container of polyhedra.
  typedef std::vector<Polyhedron<isZ> > Polyhedra;

  // The array holding pointer to planar polyhedra. The array is
  // treated as a left-aligned triangle. Each pair (i,j) with i<j is
  // located at position i*(i-1)/2+j.
  Polyhedra polyhedra;

  // An internal structure to store an interval and a back-reference.
  struct PolyBound {
    PolyBound() : interval(Interval<isZ>()), xRef(0) {};
    PolyBound(Interval<isZ> i) : interval(i), xRef(0) {};
    PolyBound(Interval<isZ> i, DomVar v) : interval(i), xRef(v) {};
    Interval<isZ> interval;
    DomVar xRef;
  };

  // Upper and lower bounds for each dimension.
  std::vector<PolyBound> bounds;

 public:
  // Create a new Two-Variable-Per-Inequality polyhedron. The optional
  // argument specifies the expected number of variables. Variables
  // have to be added explicitly by calling createVariable.
  DenseTvpi(size_t dimension=2);

  // Swap the content of this relational domain with the given
  // argument.
  void swapWith(DenseTvpi& other) {
    polyhedra.swap(other.polyhedra);
    bounds.swap(other.bounds);
  };

  // Update this TVPI system with the convex hull of this and the
  // other TVPI systems. The convex hull of a TVPI system is
  // calculated by running the planar convex hull algorithm on each
  // projection. The perm parameter maps each variable in this domain
  // to a variable in the other domain. It may not be zero and each
  // position must be a valid variable in the other domain. The
  // otherMult array maps every variable in the other domain to a
  // multiplicity value by which that variable has to be adjusted by.
  void join(const DenseTvpi& other, TvpiVar* perm, Mult* otherMult);

  // Check if the state space of this polyhedron includes that of the
  // other polyhedron. For each variable in this polyhedron, perm
  // denotes the variable in the other polyhedron to which this
  // variable should be compared to. If an index is invalidTvpiVar,
  // then the corresponding variable in this domain is not compared at
  // all.
  bool includes(const DenseTvpi& other, TvpiVar* perm) const;

  // Widen this TVPI system with repect to the other given TVPI
  // system. The perm parameter maps each variable in this domain to a
  // variable in the other domain. It may not be zero and each
  // position must be a valid variable in the other domain. The
  // extrapolate value denotes the number of steps that should be
  // extrapolated.
  void widen(const DenseTvpi& other, TvpiVar* perm,
	     const mpz_class extrapolate);

  // Return the maximum integral value in this domain that lies in the
  // direction given by the linear expression. If the found value is
  // finite, the function returns true and set the given value to the
  // maximum. In case the polyhedron extends to infinity towards the
  // given direction, the function returns false. The linear
  // components must be ascending with respect to the variables.
  bool linOpt(const std::vector<TVPIComponent>& comps,
	      mpz_class& value) const;

private:
  // Query the maximum value of a term. Returns true if a maximum
  // exists and sets value to this maximum.
  bool linOpt(const TVPIComponent& single,
	      mpz_class& value) const;

  // Query the maximum value of two terms. Returns true if a maximum
  // exists and sets value to this maximum.
  bool linOpt(const TVPIComponent& first,
	      const TVPIComponent& second,
	      mpz_class& value) const;

public:
  // Ensure that the polyhedron has room for at least n more
  // variables.  This function can be called before a set of variables
  // is created using createVariable to avoid copying the internal
  // data structure several times.
  void makeRoomForVariables(size_t headroom=0);

  // Remove all variables which are not mentioned in the given
  // array. The resulting polyhedron has sizeVarSeq active variables
  // where the ith dimension of the new polyhedron corresponds to the
  // variable v_i of the previous polyhedron, assuming v_i is the ith
  // array element. The optional argument headroom gives the number of
  // variables that can be created by createVariable without copying
  // the internal data structures.
  void projectOnto(size_t sizeVarSeq,
		   TvpiVar varSeq[],
		   size_t headroom=0);

  // Remove the last variable form the domain. Returns the xRef field
  // of the removed bound.
  DomVar removeLast();

  // Update a variable with variable that was added last. This method
  // is a special case of the projectOnto function in that it replaces
  // the projections of variable var with those of the last
  // variable. The optional argument additive can be set to true in
  // which case the variable will be updated with the join of the old
  // and the new value. This operation reduces the dimension of the
  // system by one.
  void update(TvpiVar var, bool additive = false) throw (IllegalArgument);

  // Assign the convex hull of all values from the target and the
  // source to the target. This operation corresponds to the update
  // function where "additive" is set to true. In contrast to update
  // not only the last variable but any variable can be used to
  // augment the target variable. Furthermore, the source variable is
  // not altered in any way, in particular, the dimensionality of the
  // system doesn't change.
  void augment(TvpiVar source, TvpiVar target)
    throw (IllegalArgument);

  // Set the target to the value of the source. This operation
  // corresponds to the update function where "additive" is set to
  // false. In contrast to update not only the last variable but any
  // variable can be used to augment the target variable. Furthermore,
  // the source variable is not altered in any way, in particular, the
  // dimensionality of the system doesn't change.
  void update(TvpiVar source, TvpiVar target)
    throw (IllegalArgument);

  // Query the number of variables in this system.
  TvpiVar size() const {
    return bounds.size();
  };

  // Insert a new variable into the polyhedron. There will be no bound
  // on the created variable.
  TvpiVar createVariable();

  // Insert a new variable into the polyhedron. The specified argument
  // contains the initial value of the variable.
  TvpiVar createVariable(long value);

  // Insert a new variable into the polyhedron. The specified argument
  // contains the initial value of the variable.
  TvpiVar createVariable(mpz_class value, DomVar xRef = 0);

  // Insert a new variable into the polyhedron. The specified argument
  // contains the initial interval of the variable.
  TvpiVar createVariable(Interval<isZ>& value, DomVar xRef = 0);

  inline
  TvpiVar createVariable(const Interval<isZ>& value, DomVar xRef = 0) {
    Interval<isZ> i = value;
    return createVariable(i,xRef);
  };

  // Query the cross reference of a TvpiVariable.
  DomVar xRef(TvpiVar var) const { 
    assert(var<bounds.size());
    return bounds[var].xRef;
  };

  // Set the cross reference of a TvpiVar.
  void setXRef(TvpiVar var, DomVar xRef) {
    assert(var<bounds.size());
    bounds[var].xRef=xRef;
  };
    
  // Ask for the xRefs of all relational variables. The first
  // passed-in array must have n elements, the second one (n*n-n)/2
  // elements where n is the number returned by size().
  void queryXRefs(DomVar* xRefPtr, int* relPtr) const;

private:
  // Propagate the changed bounds of the given variable to the other
  // bounds. Returns false if the system turned unsatisfiable.
  bool propagateBounds(TvpiVar var);

  // Scratch space for resultants.
  static std::vector<Inequality*> jResultants;
  static std::vector<Inequality*> kResultants;
  static std::vector<Inequality*> ineqVia;

public:
  // Insert a set of inequalities over the two specified variables. 
  // The variables can be in any order, but preferably xVar should be
  // lower. The inequalities in the input vector have to be heap
  // allocated. The ownership is transfered to this function. If the
  // inequalities all have different angles and are ordered the
  // optional flag sorted can be set to true. An exception is thrown
  // if the variables are out of bounds or they are the same. Note
  // that adding inequalities to an unsatisfiable polyhedron will
  // result in an arbitrary return value.
  Result addInequalitySet(TvpiVar xVar,
			  TvpiVar yVar,
			  std::vector<Inequality*>& ineqNew,
			  bool sorted=false)
    throw (IllegalArgument);

  // Same as above, but takes an array.
  bool addInequalities(TvpiVar xVar,
		       TvpiVar yVar,
		       size_t ineqSize,
		       Inequality* ineqNew[],
		       bool sorted=false)
    throw (IllegalArgument) {
    std::vector<Inequality*> ineqVec;
    ineqVec.reserve(ineqSize);
    for (size_t i=0; i<ineqSize; i++)
      ineqVec.push_back(ineqNew[i]);
    return addInequalitySet(xVar, yVar, ineqVec, sorted);
  }

  // Calculate a "dimensionality" of this polyhedron. The result is
  // simply the sum of the dimensionality of all planar polyhedra. The
  // idea is that comparing the dimensionality of this TVPI system
  // with another, gives an indication if widening should be applied. 
  // This should always be correct assuming that the dimensionality of
  // a sequence of polyhedra grows monotonically.
  int dimensionality() const {
    int dim=0;
    size_t size=bounds.size();

    typename Polyhedra::const_iterator iter=polyhedra.begin();
    for (TvpiVar y=0; y<size; y++)
      for (TvpiVar x=0; x<y; x++, iter++) {
	assert(iter!=polyhedra.end());
	dim+=iter->dimensionality(bounds[x].interval, bounds[y].interval);
      };
    return dim;
  }

  // Approximate an n-dimensional inequality. The given linear
  // inequality is approximated with TVPI constraints. The input is
  // given as a vector of variables and coefficients together with a
  // constant. The inequality a_1 x_1 + ... + a_n x_n + c <=0 will be
  // approximated by sets of inequalities over two variables, if
  // isEquality is true, replace "<=" with "=". No a_i may be zero, no
  // x_i may appear twice. The inequalities that result from the
  // approximation are inserted into the domain.
  Result approximateInequality(const std::vector<LinComponent>& comps,
			       mpz_class constant,
			       bool isEquality = false);
  
  // Update the value of one variable in terms of another (or the
  // same).  I.e. execute a1 x1 = a2 x2 + [cMin..cMax] . If x1 is
  // invalidVariable then a fresh variable is created. In every case
  // the function returns the variable x1 as result.
  TvpiVar updateVariable(mpz_class a1, TvpiVar x1,
			  mpz_class a2, TvpiVar x2,
			  mpz_class cMin, mpz_class cMax) {
    assert(cMin<=cMax);
    std::vector<Inequality*> set;
    set.push_back(new Inequality(-a2, a1, cMax));
    set.push_back(new Inequality(a2, -a1, -cMin));
    TvpiVar temp=createVariable();
#ifndef NDEBUG
    bool satisfiable=
#endif
      addInequalitySet(x2, temp, set, true);
    assert(satisfiable);
    if (x1==invalidTvpiVar) return temp;
    update(x1);
    return x1;
  };

  // Intersect the range of the variable with the given bound. Returns
  // false if the domain turned unsatisfiable.
  bool intersectBound(TvpiVar var, Interval<isZ>& i) {
    assert(var<bounds.size());
    bounds[var].interval.intersect(i);
    bool res = true;
    if (bounds[var].interval.hasChanged()) res = propagateBounds(var);
    return res;
  };

  // Multiply the value of the given variable by 2^m.
  void stretch(TvpiVar var, Mult m);

  // Calculate the minimal multiplicity of the given variable, using
  // the equality relationships in the domain.
  Mult seekMultFromEquality(TvpiVar var, Mult* mults);

  // Query the interval of a specific variable.
  const Interval<isZ>& getInterval(TvpiVar var) const {
    assert(var<bounds.size());
    return bounds[var].interval;
  };

  // Query the maximum value of a variable.
  mpq_class maxValue(TvpiVar var) const {
    assert(var<bounds.size());
    return bounds[var].interval.getUpper();
  };

  // Query the minimum value of a variable.
  mpq_class minValue(TvpiVar var) const {
    assert(var<bounds.size());
    return bounds[var].interval.getLower();
  };

  // Check if the variable has an upper bound.
  bool upperIsFinite(TvpiVar var) const {
    assert(var<bounds.size());
    return !bounds[var].interval.upperIsInfinite();
  };

  // Check if the variable has a lower bound.
  bool lowerIsFinite(TvpiVar var) const {
    assert(var<bounds.size());
    return !bounds[var].interval.lowerIsInfinite();
  };

  // Substitute the variable y in ts + ky + c = 0 for another if the
  // other variable is in an equality relationship with y. The
  // function takes, k and c and returns a variable x that has a
  // smaller index than y. The variables k and c are set such that ts
  // + kx + c = 0. In case the initial k is greater than one, a variable
  // x might be found that is in equality relation but cannot be substituted
  // without a k<1. In this case, m is set to a multiplication factor
  // with which the whole equation (except k and c) has to be multiplied
  // before the substitution can take place.
  // The varialbe y is returned if no equal
  // variable could be found. In this case, k and c remain unchanged.
  TvpiVar substEquality(TvpiVar y,
			mpz_class& k,
			mpz_class& c,
			mpz_class& m) const;

  // Retrieve a projection from the TVPI system (for
  // visualization). The variables must be strictly increasing.
  const Polyhedron<isZ>& getProjection(TvpiVar x, TvpiVar y) const;

  void output(std::ostream& stream, Mult *mults) const;

  friend std::ostream& operator<<(std::ostream& stream,
				  const DenseTvpi& t) {
    t.output(stream, 0);
    return stream;
  };

  // Draw a little picture showing where inequalities are and if the
  // polyhedra are shared.
  std::ostream& showDist(std::ostream& stream) const;

  void sane() const;
};

};
#endif // __TVPI_H
