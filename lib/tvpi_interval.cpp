/* Derived from the TVPI library (original file: interval.cpp).
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
/* interval.cpp
 * Operations on intervals.
 */

#include <crab/domains/tvpi_impl/interval.hh>
#include <assert.h>

#undef DEBUG_UPDATE

#ifdef DEBUG_MEMORY
#include <iostream>
#include <string.h>
#endif

using std::cerr;
using std::endl;

namespace Tvpi {

  Mult calcMult(const mpz_class& c) {
    unsigned long log2 = mpz_scan1(c.get_mpz_t(), 0);
    if (log2>(unsigned int) invalidMult) return invalidMult;
    return (Mult) log2;
  };


// Copy constructor. Creates a deep copy of the interval.
template<bool isZ>
IntervalImpl<isZ>::IntervalImpl(IntervalImpl<isZ>& i) {
  upper.bound=i.upper.bound;
  lower.bound=i.lower.bound;
  refCount=1;
}

// Compare the upper bound with another interval. Returns a positive
// value if this upper bound is greater than that of the other, 0 if
// their are equal, a negative number otherwise. Infinities compare as
// equal.
template<>
int IntervalImpl<false>::compareUpper(const IntervalImpl<false>& other) const {
  if (sgn(upper.bound.get_den())==0) {
    if (sgn(other.upper.bound.get_den())==0) {
      // both are infinity (assumption: +\infty)
      return 0;
    } else {
      // the other is less than infinity
      return 1;
    };
  } else {
    if (sgn(other.upper.bound.get_den())==0) {
      // this number is finite, the other is not
      return -1;
    } else {
      // both are finite
      return cmp(upper.bound, other.upper.bound);
    };
  };
}

template<>
int IntervalImpl<true>::compareUpper(const IntervalImpl<true>& other) const {
  if (upper.bound.isInfinite) {
    if (other.upper.bound.isInfinite) {
      // both are infinity (assumption: +\infty)
      return 0;
    } else {
      // the other is less than infinity
      return 1;
    };
  } else {
    if (other.upper.bound.isInfinite) {
      // this number is finite, the other is not
      return -1;
    } else {
      // both are finite
      return cmp(upper.bound.value, other.upper.bound.value);
    };
  };
}

// Same for the lower bound. Returns a positive number if this lower
// bound is greater than the other.
template<>
int IntervalImpl<false>::compareLower(const IntervalImpl<false>& other) const {
  if (sgn(lower.bound.get_den())==0) {
    if (sgn(other.lower.bound.get_den())==0) {
      // both are infinity (assumption: -\infty)
      return 0;
    } else {
      // the other is greater than minus infinity
      return -1;
    };
  } else {
    if (sgn(other.lower.bound.get_den())==0) {
      // this number is finite, the other is not
      return 1;
    } else {
      // both are finite
      return cmp(lower.bound, other.lower.bound);
    };
  };
}

template<>
int IntervalImpl<true>::compareLower(const IntervalImpl<true>& other) const {
  if (lower.bound.isInfinite) {
    if (other.lower.bound.isInfinite) {
      // both are infinity (assumption: -\infty)
      return 0;
    } else {
      // the other is greater than minus infinity
      return -1;
    };
  } else {
    if (other.lower.bound.isInfinite) {
      // this number is finite, the other is not
      return 1;
    } else {
      // both are finite
      return cmp(lower.bound.value, other.lower.bound.value);
    };
  };
}

// Returns true if the interval has no upper bound.
template<bool isZ>
bool IntervalImpl<isZ>::upperIsInfinite() const {
  return upper.isInfinite();
}


// Returns true if the interval has no lower bound.
template<bool isZ>
bool IntervalImpl<isZ>::lowerIsInfinite() const {
  return lower.isInfinite();
}


// Returns true if the interval has no upper bound.
template<bool isZ>
bool IntervalImpl<isZ>::upperIsFinite() const {
  return !upper.isInfinite();
}


// Returns true if the interval has no lower bound.
template<bool isZ>
bool IntervalImpl<isZ>::lowerIsFinite() const {
  return !lower.isInfinite();
}



// Return true if this interval includes the other.
template<bool isZ>
bool IntervalImpl<isZ>::includes(const IntervalImpl<isZ>& i) const {
  if (compareLower(i)>0) return false;
  if (compareUpper(i)<0) return false;
  return true;
}

template<>
bool IntervalImpl<false>::includes(const mpz_class& v) const {
  if (!upper.isInfinite() && mpq_class(v)>upper.bound) return false;
  if (!lower.isInfinite() && mpq_class(v)<lower.bound) return false;
  return true;
}

template<>
bool IntervalImpl<true>::includes(const mpz_class& v) const {
  if (!upper.isInfinite() && v>upper.bound.value) return false;
  if (!lower.isInfinite() && v<lower.bound.value) return false;
  return true;
}

  // Divide the interval by 2^m, m>1. Returns false if the interval
  // becomes empty.
  template<>
  bool IntervalImpl<false>::enforceMult(Mult m) {
    if (!upper.isInfinite())
      mpq_div_2exp(upper.bound.get_mpq_t(),
		   upper.bound.get_mpq_t(),
		   (unsigned long int) m);
    if (!lower.isInfinite())
      mpq_div_2exp(lower.bound.get_mpq_t(),
		   lower.bound.get_mpq_t(),
		   (unsigned long int) m);
    return true;
  }

  template<>
  bool IntervalImpl<true>::enforceMult(Mult m) {
    int infiniteBounds = 2;
    if (!upper.bound.isInfinite) {
      mpz_fdiv_q_2exp(upper.bound.value.get_mpz_t(),
		      upper.bound.value.get_mpz_t(),
		      (unsigned long int) m);
      infiniteBounds--;
    };
    if (!lower.bound.isInfinite) {
      mpz_cdiv_q_2exp(lower.bound.value.get_mpz_t(),
		      lower.bound.value.get_mpz_t(),
		      (unsigned long int) m);
      infiniteBounds--;
    };
    if (infiniteBounds==0) return (lower.bound.value<=upper.bound.value);
    return true;
  }


// Return the size of the interval.
template<>
Bound<false> IntervalImpl<false>::getSize() {
  if (upper.isInfinite() || lower.isInfinite()) return Bound<false>(); else
    return Bound<false>(mpq_class(upper.bound-lower.bound));
}

// Return the size of the interval.
template<>
Bound<true> IntervalImpl<true>::getSize() {
  if (upper.isInfinite() || lower.isInfinite()) return Bound<true>(); else
    return Bound<true>(upper.bound.value-lower.bound.value);
}

// Compare interval for equality.
template<>
bool IntervalImpl<false>::operator==(const IntervalImpl& other) const {
  return mpq_equal(lower.bound.get_mpq_t(), other.lower.bound.get_mpq_t()) 
    &&   mpq_equal(upper.bound.get_mpq_t(), other.upper.bound.get_mpq_t());
}

template<>
bool IntervalImpl<true>::operator==(const IntervalImpl& other) const {
  if (lower.bound.isInfinite!=other.lower.bound.isInfinite) return false;
  if (upper.bound.isInfinite!=other.upper.bound.isInfinite) return false;
  if (!lower.bound.isInfinite &&
      lower.bound.value!=other.lower.bound.value) return false;
  if (!upper.bound.isInfinite &&
      upper.bound.value!=other.upper.bound.value) return false;
  return true;
}

defMemDbg(template<bool isZ>,IntervalImpl<isZ>,i,I)



// Interval: a reference-counting wrapper around IntervalImpl.
  
// see polyhedron.cpp:Polyhedron
IntervalImpl<false>* unboundedQIntervalImpl = ::new IntervalImpl<false>();
IntervalImpl<true>* unboundedZIntervalImpl = ::new IntervalImpl<true>();
  
// Create an unbounded interval.
template<>
Interval<false>::Interval() {
  // Make a shallow copy of the unbounded interval.
  inter=unboundedQIntervalImpl;
  inter->refCount++;
}

template<>
Interval<true>::Interval() {
  inter=unboundedZIntervalImpl;
  inter->refCount++;
}

// Create a new reference to this interval.
template<bool isZ>
Interval<isZ>::Interval(const Interval& i) {
  inter=i.inter;
  inter->refCount++;
}

template<bool isZ>
Interval<isZ>::~Interval() {
  assert(inter);
  assert(inter->refCount>0);
  if (!--inter->refCount) delete inter;
}

template<bool isZ>
Interval<isZ>& Interval<isZ>::operator=(const Interval<isZ>& other) {
  assert(inter);
  assert(inter->refCount>0);
  assert(other.inter);
  other.inter->refCount++;
  if (!--inter->refCount) delete inter;
  inter=other.inter;
  return *this;
}

// Calculate an interval that includes both i1 and i2.
template<bool isZ>
Interval<isZ>::Interval(const Interval<isZ>& i1, const Interval<isZ>& i2) {
  assert(i1.inter);
  assert(i2.inter);
  if (i1.inter==i2.inter) {
    // This is the copy constructor.
    inter=i1.inter;
    inter->refCount++;
  } else {
    // Assume that both intervals are different
    int cmpL=i1.inter->compareLower(*i2.inter);
    int cmpU=i1.inter->compareUpper(*i2.inter);
    if (cmpL<0) {
      if (cmpU<0) inter=new IntervalImpl<isZ>(i1.inter->lower,
					      i2.inter->upper); 
      else {
	inter=i1.inter;
	inter->refCount++;
      };
    } else if (cmpL==0) {
      if (cmpU<0) {
	inter=i2.inter;
	inter->refCount++;
      } else {
	inter=i1.inter;
	inter->refCount++;
      };
    } else {
      assert(cmpL>0);
      if (cmpU>0) inter=new IntervalImpl<isZ>(i2.inter->lower,
					      i1.inter->upper); 
      else {
	inter=i2.inter;
	inter->refCount++;
      }
    }
  }
}

// Access the upper value for reading.
template<>
mpq_class Interval<false>::getUpper() const {
 return inter->upper.bound;
}

template<>
mpq_class Interval<true>::getUpper() const {
 return mpq_class(inter->upper.bound.value,
		  inter->upper.bound.isInfinite ? 0 : 1);
}

// Access the lower value for reading.
template<>
mpq_class Interval<false>::getLower() const {
 return inter->lower.bound;
}

template<>
mpq_class Interval<true>::getLower() const {
 return mpq_class(inter->lower.bound.value,
		  inter->lower.bound.isInfinite ? 0 : 1);
}


// Join the other interval into this.
template<bool isZ>
void Interval<isZ>::join(const Interval<isZ>& other) {
  assert(inter);
  assert(other.inter);
  if (inter==other.inter) return;
  int resU = inter->compareUpper(*other.inter);
  int resL = inter->compareLower(*other.inter);
  if (resU>=0 && resL<=0) return;
  makeUnique();
  if (resU<0) inter->upper=other.inter->upper;
  if (resL>0) inter->lower=other.inter->lower;
}

// Update the upper bound of this interval. If the given bound
// causes a change, the `changed' flag is set.
template<>
bool Interval<false>::updateUpper(mpq_class bound) {
  assert(inter);
  if (!inter->upperIsInfinite() && inter->upper.bound<=bound) return false;
  makeUnique();
  mpq_swap(inter->upper.bound.get_mpq_t(),bound.get_mpq_t());
  inter->changed=true;
  return true;
}

template<>
bool Interval<true>::updateUpper(mpq_class bound) {
  assert(inter);
  mpz_class newBound;
  mpz_fdiv_q(newBound.get_mpz_t(),
	     bound.get_num().get_mpz_t(), bound.get_den().get_mpz_t());
  if (!inter->upperIsInfinite() && inter->upper.bound.value<=newBound)
    return false;
#ifdef DEBUG_UPDATE
  cerr << "updating upper bounds from ";
  if (inter->upper.bound.isInfinite) cerr << "oo";
  else cerr << inter->upper.bound.value;
  cerr << " to " << newBound << ", i.e. " << bound << endl;
#endif // DEBUG_UPDATE
  makeUnique();
  inter->upper.bound.isInfinite=false;
  mpz_swap(inter->upper.bound.value.get_mpz_t(),newBound.get_mpz_t());
  inter->changed=true;
  return true;
}

// Update the lower bound of this interval. If the given bound
// causes a change, the `changed' flag is set.
template<>
bool Interval<false>::updateLower(mpq_class bound) {
  assert(inter);
  if (!inter->lowerIsInfinite() && inter->lower.bound>=bound) return false;
  makeUnique();
  mpq_swap(inter->lower.bound.get_mpq_t(),bound.get_mpq_t());
  inter->changed=true;
  return true;
}

template<>
bool Interval<true>::updateLower(mpq_class bound) {
  assert(inter);
  mpz_class newBound;
  mpz_cdiv_q(newBound.get_mpz_t(),
	     bound.get_num().get_mpz_t(), bound.get_den().get_mpz_t());
  if (!inter->lowerIsInfinite() && inter->lower.bound.value>=newBound)
    return false;
#ifdef DEBUG_UPDATE
  cerr << "updating lower bounds from "; 
  if (inter->lower.bound.isInfinite) cerr << "-oo"; 
  else cerr << inter->lower.bound.value;
  cerr << " to " << newBound << ", i.e. " << bound << endl;
#endif // DEBUG_UPDATE
  makeUnique();
  inter->lower.bound.isInfinite=false;
  mpz_swap(inter->lower.bound.value.get_mpz_t(),newBound.get_mpz_t());
  inter->changed=true;
  return true;
}

// Check if the upper bound is smaller than the lower bound.
template<>
bool Interval<false>::isEmpty() const {
  if (inter->upperIsInfinite()) return false;
  if (inter->lowerIsInfinite()) return false;
  return (inter->upper.bound<inter->lower.bound);
}

template<>
bool Interval<true>::isEmpty() const {
  if (inter->upperIsInfinite()) return false;
  if (inter->lowerIsInfinite()) return false;
  return (inter->upper.bound.value<inter->lower.bound.value);
}

std::ostream& operator<< (std::ostream& stream,
			  const Interval<false>& i) {
  using namespace std;
  stream << "[";
  if (i.inter->lowerIsInfinite()) 
    stream << "-infty"; else stream << i.inter->lower.bound;
  stream << "..";
  if (i.inter->upperIsInfinite()) 
    stream << "+infty"; else stream << i.inter->upper.bound;
  stream << "]";
  return stream;
}

std::ostream& operator<< (std::ostream& stream,
			  const Interval<true>& i) {
  using namespace std;
  stream << "[";
  if (i.inter->lowerIsInfinite()) stream << "-infty";
  else stream << i.inter->lower.bound.value;
  stream << "..";
  if (i.inter->upperIsInfinite()) stream << "+infty"; 
  else stream << i.inter->upper.bound.value;
  stream << "]";
  return stream;
}


  Mult calcMult(const Interval<true>& i) {
    assert(i.inter);
    if (!i.isSingleton()) return 0;
    if (sgn(i.inter->upper.bound.value)==0) return 0;
    unsigned long log2 = mpz_scan1(i.inter->upper.bound.value.get_mpz_t(), 0);
    if (log2>255) return 255;
    return (Mult) log2;
  }

// Make a deep copy if this interval is shared.
template<bool isZ>
void Interval<isZ>::makeUnique() {
  assert(inter);
  assert(inter->refCount>0);
  if (inter->refCount==1) return;
  // Make a deep copy.
  inter->refCount--;
  assert(inter->refCount>0);
  inter=new IntervalImpl<isZ>(*inter);
}

template struct Bound<false>;
template struct Bound<true>;
template class IntervalImpl<false>;
template class IntervalImpl<true>;
template class Interval<false>;
template class Interval<true>;

} // namespace Tvpi
