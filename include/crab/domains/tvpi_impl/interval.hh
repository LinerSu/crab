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
/* interval.hh
 * Operations on intervals.
 */

#ifndef __INTERVAL_H
#define __INTERVAL_H

#include "common.hh"
#include "memory.hh"
#include "tvpiexception.hh"
#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <assert.h>
#include <limits.h>


// Note:

// The bounds have a flag denoting if they are infinite. That is a
// misconception since there is only one flag for +\infty and one flag
// for -\infty which are clearly different values. The flag in the
// bound should probably be undefined or something. This shows in the
// isSingleton test of the interval since here a simple equality test
// between the upper and lower bounds doesn't work if both bounds are
// infinite.

namespace Tvpi {

  // Calculate the number of zero bits of a number. Returns
  // invalidMult if c is zero.
  Mult calcMult(const mpz_class& c);

  // The maximum value for Mult is defined here. It serves as a
  // special token indicating the arbitrary multiplicity of zero.
  const Mult invalidMult = UCHAR_MAX;

template<bool> class Interval;

// Declare functions to output intervals.
std::ostream& operator<< (std::ostream& stream,
			  const Interval<false>& i);
std::ostream& operator<< (std::ostream& stream,
			  const Interval<true>& i);

// Augment the integer class with a flag that says that the value is
// actually infinity.
struct mpzi_class {
  mpz_class value;
  bool isInfinite;
  mpzi_class(): isInfinite(true) {};
  mpzi_class(mpz_class v): value(v), isInfinite(false) {};
};

template<bool isZ>
struct Bound {
  mpq_class bound;
  Bound(): bound(1,0) {};
  Bound(mpz_class val): bound(val, 1) {};
  Bound(mpq_class val): bound(val) {};

  bool isInfinite() const { return sgn(bound.get_den())==0; };

  void negate() {
    mpq_neg(bound.get_mpq_t(), bound.get_mpq_t());
  };

  bool operator<(Bound<isZ>& other) const { 
    return (other.isInfinite() ? true :
	    (isInfinite() ? false :
	     bound<other.bound));
  };

  bool operator==(Bound<isZ>& other) {
    if (isInfinite() && other.isInfinite()) return true;
    if (isInfinite()) return false;
    if (other.isInfinite()) return false;
    return bound==other.bound;
  };

  // Widen this interval with respect to the other.
  void widen(const Bound<isZ>& other, const mpz_class& extrapolate) {
    assert(sgn(bound.get_den())!=0);
    if (sgn(extrapolate)==0) {
      if (bound!=other.bound) {
	bound.get_num()=1;
	bound.get_den()=0;
      }
    } else {
      mpq_class dist = bound-other.bound;
      bound+=dist*mpq_class(extrapolate);
    }
  };

  // Multiply a bound by 2^m.
  inline void stretch(Mult m) {
    if (!isInfinite()) {
      mpz_mul_2exp(bound.get_num().get_mpz_t(),
		   bound.get_num().get_mpz_t(), (unsigned long int) m);
      bound.canonicalize();
    }
  };

  // Subtract another Bound.
  inline Bound<isZ> operator-=(const Bound<isZ>& rhs) {
    if (sgn(bound.get_den())!=0) {
      if (sgn(rhs.bound.get_den())!=0) 
	bound-=rhs.bound;
      else
	bound.get_den()=0;
    };
    return *this;
  };
  
  // Add another Bound.
  inline Bound<isZ> operator+=(const Bound<isZ>& rhs) {
    if (sgn(bound.get_den())!=0) {
      if (sgn(rhs.bound.get_den())!=0) 
	bound+=rhs.bound;
      else
	bound.get_den()=0;
    };
    return *this;
  };

  // Multiply a bound with a constant.
  inline Bound<isZ> operator*(const mpz_class& num) const {
    if (sgn(bound.get_den())!=0) return Bound<isZ>(mpq_class(bound*num));
    return Bound<isZ>();
  };

  // Swap with another bound.
  void swap(Bound<isZ>& b) {
    mpq_swap(bound.get_mpq_t(), b.bound.get_mpq_t());
  };

  // Check if the bound is zero.
  bool isZero() const {
    if (sgn(bound.get_den())==0) return false; else return sgn(bound)==0;
  }

};

template<>
struct Bound<true> {
  mpzi_class bound;
  Bound(): bound() {};
  Bound(mpz_class val) : bound(val) {};

  bool isInfinite() const { return bound.isInfinite; };

  void negate() {
    mpz_neg(bound.value.get_mpz_t(), bound.value.get_mpz_t());
  };

  bool operator<(Bound<true>& other) { 
    return (other.isInfinite() ? true :
	    (isInfinite() ? false :
	     bound.value<other.bound.value));
  };

  bool operator==(Bound<true>& other) {
    if (isInfinite() && other.isInfinite()) return true;
    if (isInfinite()) return false;
    if (other.isInfinite()) return false;
    return bound.value==other.bound.value;
  };

  // Widen this interval with respect to the other.
  void widen(const Bound<true>& other, const mpz_class& extrapolate) {
    assert(!bound.isInfinite);
    if (sgn(extrapolate)==0) {
      if (bound.value!=other.bound.value) bound.isInfinite=true;
    } else {
      mpz_class dist = bound.value-other.bound.value;
      bound.value+=dist*extrapolate;
    }
  };

  // Multiply a bound by 2^m.
  inline void stretch(Mult m) {
    if (!bound.isInfinite)
      mpz_mul_2exp(bound.value.get_mpz_t(),
		   bound.value.get_mpz_t(), (unsigned long int) m);
  };
    
  // Subtract another Bound.
  inline Bound<true> operator-=(const Bound<true>& rhs) {
    if (!bound.isInfinite) {
      if (!rhs.bound.isInfinite) 
	bound.value-=rhs.bound.value;
      else
	bound.isInfinite=true;
    };
    return *this;
  };
  
  // Add another Bound.
  inline Bound<true> operator+=(const Bound<true>& rhs) {
    if (!bound.isInfinite) {
      if (!rhs.bound.isInfinite) 
	bound.value+=rhs.bound.value;
      else
	bound.isInfinite=true;
    };
    return *this;
  };

  // Multiply a bound with a constant.
  inline Bound<true> operator*(const mpz_class& num) const {
    if (!bound.isInfinite) return Bound<true>(bound.value*num);
    return Bound<true>();
  };

  // Swap with another bound.
  void swap(Bound<true>& b) {
    mpz_swap(bound.value.get_mpz_t(), b.bound.value.get_mpz_t());
    bool tmp = bound.isInfinite;
    bound.isInfinite = b.bound.isInfinite;
    b.bound.isInfinite = tmp;
  };

  // Check if the bound is zero.
  bool isZero() const {
    if (bound.isInfinite) return false; else return sgn(bound.value)==0;
  }

};


template<bool isZ> 
class IntervalImpl {
  friend class Interval<isZ>;

  // The upper and lower bound.
  Bound<isZ> lower;
  Bound<isZ> upper;

  // The number of Interval classes that reference this instance.
  size_t refCount;

  // A flag tracking if this interval has changed recently.
  bool changed;

  IntervalImpl(mpz_class l, mpz_class u):
    lower(l), upper(u), refCount(1), changed(false) {};

  IntervalImpl(Bound<isZ> l, Bound<isZ> u):
    lower(l), upper(u), refCount(1), changed(false) {};

public:
  // Create an unbounded interval.
  IntervalImpl() : lower(), upper(), refCount(1), changed(false) {};

  // Copy constructor. Creates a deep copy of the interval.
  IntervalImpl(IntervalImpl& i);

  // Compare the upper bound with another interval. Returns a positive
  // value if this upper bound is greater than that of the other, 0 if
  // their are equal, a negative number otherwise.
  int compareUpper(const IntervalImpl& other) const;

  // Same for the lower bound. Returns a positive number if this lower
  // bound is greater than the other.
  int compareLower(const IntervalImpl& other) const;

  // Returns true if the interval has no upper bound.
  bool upperIsInfinite() const;

  // Returns true if the interval has no lower bound.
  bool lowerIsInfinite() const;

  // Returns true if the interval has an upper bound.
  bool upperIsFinite() const;

  // Returns true if the interval has a lower bound.
  bool lowerIsFinite() const;

  // Return true if the other interval is smaller.
  bool includes(const IntervalImpl& i) const;

  // Return true if the value is in the interval.
  bool includes(const mpz_class& v) const;

  // Compare interval for equality.
  bool operator==(const IntervalImpl& other) const;

  // Divide the interval by 2^m, m>1. Returns false if the interval
  // becomes empty.
  bool enforceMult(Mult m);

  // Invert an interval.
  inline void negate() {
    upper.negate();
    lower.negate();
    lower.swap(upper);
  };

  // Return the size of the interval.
  Bound<isZ> getSize();

  declMemDbg;

  friend std::ostream& operator<< (std::ostream& stream,
				   const Interval<false>& i);
  friend std::ostream& operator<< (std::ostream& stream,
				   const Interval<true>& i);

  friend Mult calcMult(const Interval<true>& i);
};


template<bool isZ> 
class Interval {

  // The (possibly shared) information on this interval.
  IntervalImpl<isZ>* inter;

  // Wrap an IntervalImpl pointer.
  Interval(IntervalImpl<isZ>* i) : inter(i) {};

public:
  // Create an unbounded interval.
  Interval();

  // Create a new reference to this interval.
  Interval(const Interval& i);

  ~Interval();

  Interval& operator=(const Interval& other);

  // Calculate an interval that includes both i1 and i2.
  Interval(const Interval& i1, const Interval& i2);

  // Create an interval containing a specific value.
  Interval(long value) {
    mpz_class val(value);
    inter=new IntervalImpl<isZ>(val, val);
  };

  // Create an interval containing a specific range.
  Interval(long lower, long upper) {
    assert(lower<upper);
    mpz_class low(lower);
    mpz_class upp(upper);
    inter=new IntervalImpl<isZ>(low, upp);
  };

  // Create an interval containing a specific value.
  Interval(mpz_class val) {
    inter=new IntervalImpl<isZ>(val, val);
  };

  // Access the upper value for reading.
  mpq_class getUpper() const;

  // Access the lower value for reading.
  mpq_class getLower() const;

  // Join in another interval.
  void join(const Interval& other);

  // Widen this interval with respect to the other.
  inline void widen(const Interval& other, const mpz_class& extrapolate=0) {
    assert(inter);
    // Don't touch anything if the pointers are the same.
    if (inter==other.inter) return;
    if (*inter==*other.inter) return;
    makeUnique();
    if (inter->upperIsFinite())
      inter->upper.widen(other.inter->upper, extrapolate);
    if (inter->lowerIsFinite())
      inter->lower.widen(other.inter->lower, extrapolate);
  };

  // Tighten the upper bound of this interval, return true if the
  // interval changed. In this case, the changed flag is set.
  bool updateUpper(mpq_class bound);

  // Tighten the lower bound of this interval, return true if the
  // interval changed. In this case, the changed flag is set.
  bool updateLower(mpq_class bound);

  // Check if this interval has changed.
  inline bool hasChanged() { return inter->changed; };
  
  // Clear the changed flag.
  inline void clearChanged() { inter->changed=false; };

  // Check if the upper bound is smaller than the lower bound.
  bool isEmpty() const;
  
  // Returns true if the interval has no upper bound.
  bool upperIsInfinite() const { return inter->upperIsInfinite(); };

  // Returns true if the interval has no lower bound.
  bool lowerIsInfinite() const { return inter->lowerIsInfinite(); };

  // Returns true if the interval has an upper bound.
  bool upperIsFinite() const { return inter->upperIsFinite(); };

  // Returns true if the interval has a lower bound.
  bool lowerIsFinite() const { return inter->lowerIsFinite(); };

  // Return true if the other interval is smaller.
  bool includes(const Interval& i) const {
    assert(inter);
    assert(i.inter);
    return inter->includes(*i.inter);
  };

  // Return true if the other interval is smaller.
  bool includes(const mpz_class& v) const {
    assert(inter);
    return inter->includes(v);
  };

  void intersect(const Interval& other) {
    if (other.upperIsFinite()) {
      if (upperIsInfinite() || other.inter->upper<inter->upper) {
	makeUnique();
	inter->upper=other.inter->upper;
	inter->changed=true;
      }
    };
    if (other.lowerIsFinite()) {
      if (lowerIsInfinite() || inter->lower<other.inter->lower) {
	makeUnique();
	inter->lower=other.inter->lower;
	inter->changed=true;
      }
    }
  };

  // Compare interval for equality.
  bool operator==(const Interval& other) const { 
    if (this->inter==other.inter) return true;
    else return (*inter)==(*other.inter);
  };

  // Compare interval for disequality.
  bool operator!=(const Interval& other) const {
    return !(*this==other);
  };

  // Invert an interval.
  inline void negate() {
    assert(inter);
    makeUnique();
    inter->negate();
  };

  inline void swap(Interval& other) {
    IntervalImpl<isZ>* tmp=inter;
    inter=other.inter;
    other.inter=tmp;
  };

  // Return an interval which is multiplied by 2^m.
  inline const Interval<isZ> stretch(Mult m) const {
    if (m==0) return *this;
    assert(inter);
    assert(inter->refCount>0);
    IntervalImpl<isZ>* i = new IntervalImpl<isZ>(*inter);
    i->upper.stretch(m);
    i->lower.stretch(m);
    return Interval<isZ>(i);
  };

  // Divide an interval by 2^m. Return false if interval
  // became unsatisfiable.
  inline bool enforceMult(Mult m) {
    assert(inter);
    if (m==0) return true;
    makeUnique();
    return inter->enforceMult(m);
  };

  // Subtract another Interval point-wise.
  inline Interval<isZ> operator-=(const Interval<isZ>& rhs) {
    makeUnique();
    // [5..7] - [5..7] = [5..7] + (-[5..7]) = [5..7] + [-7..-5] = [-2..2]
    inter->upper-=rhs.inter->lower;
    inter->lower-=rhs.inter->upper;
    return *this;
  };
  
  // Add another Interval point-wise.
  inline Interval<isZ> operator+=(const Interval<isZ>& rhs) {
    makeUnique();
    inter->upper+=rhs.inter->upper;
    inter->lower+=rhs.inter->lower;
    return *this;
  };

  // Multiply an interval with a constant point-wise.
  inline Interval<isZ> operator*(const mpz_class& num) const {
    assert(inter);
    assert(inter->refCount>0);
    IntervalImpl<isZ>* i;
    if (sgn(num)<0) {
      i = new IntervalImpl<isZ>(Bound<isZ>(inter->upper*num), 
				Bound<isZ>(inter->lower*num));
    } else {
      i = new IntervalImpl<isZ>(Bound<isZ>(inter->lower*num), 
				Bound<isZ>(inter->upper*num));
    };
    return Interval(i);
  };

  // Return the size of the interval.
  Bound<isZ> getSize() {
    return inter->getSize();
  };

  friend std::ostream& operator<< (std::ostream& stream,
				   const Interval<false>& i);
  friend std::ostream& operator<< (std::ostream& stream,
				   const Interval<true>& i);

  bool isSingleton() const {
    assert(inter);
    if (inter->upperIsInfinite() || inter->lowerIsInfinite()) return false;
    return inter->upper==inter->lower;
  };

  bool isZero() const {
    assert(inter);
    return (inter->upper.isZero() && inter->lower.isZero());
  };
    
  friend Mult calcMult(const Interval<true>& i);

private:
  // Make a deep copy if this interval is shared.
  void makeUnique();

};

}


#endif // __INTERVAL_H
