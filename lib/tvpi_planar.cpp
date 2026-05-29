/* Derived from the TVPI library (original file: planar.cpp).
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
/* planar.c
 * Operations on  planar polyhedra.
 */

#include <crab/domains/tvpi_impl/planar.hh>
#include <crab/domains/tvpi_impl/tvpiexception.hh>
#include <gmp.h>
#include <assert.h>

#undef DEBUG_INTERSECT
#undef DEBUG_ENTAIL
#undef DEBUG_CUT
#undef DEBUG_UPDATE
#undef DEBUG_EXTREME
#undef DEBUG_REFINE

#ifdef DEBUG_MEMORY
#include <iostream>
using std::cerr;
#endif

#ifdef NDEBUG

#undef DEBUG_INTERSECT
#undef DEBUG_ENTAIL
#undef DEBUG_CUT
#undef DEBUG_UPDATE
#else

#include <iostream>
using std::cerr;
using std::endl;
using std::iostream;

#endif // NDEBUG

namespace Tvpi {

  // Set the coefficients that should be approximated. This is a
  // no-op of the coefficients are the same as they were before the
  // call.
  size_t Convergents::initialize(mpz_class aa, mpz_class bb) {
    mpz_abs(aa.get_mpz_t(), aa.get_mpz_t());
    mpz_abs(bb.get_mpz_t(), bb.get_mpz_t());
    if (a==aa && b==bb) return approx.size();

    a=aa;
    b=bb;

    mpz_class curN = a;
    mpz_class curD = b;
    approx.clear();
    if (sgn(curD)==0) return 0;

    approx.push_back(Approximation());
    // Calculate q_0 and q_1.
    mpz_tdiv_qr(approx[0].term.get_mpz_t(),
		curN.get_mpz_t(),
		curN.get_mpz_t(),
		curD.get_mpz_t());
    mpz_swap(curD.get_mpz_t(), curN.get_mpz_t());

    // Set A_0, and B_0.
    approx[0].A=approx[0].term;
    approx[0].B=1;
    if (sgn(curD)==0) return 1;

    approx.push_back(Approximation());
    mpz_tdiv_qr(approx[1].term.get_mpz_t(),
		curN.get_mpz_t(),
		curN.get_mpz_t(),
		curD.get_mpz_t());
    mpz_swap(curD.get_mpz_t(), curN.get_mpz_t());
    
    // Set A_1, and B_1.
    approx[1].A=approx[0].term*approx[1].term+1;
    approx[1].B=approx[1].term;

    // Last convergent already caculated?
    while (sgn(curD)!=0) {
      
      approx.push_back(Approximation());
      
      // Indices m, m-1 and m-2.
      std::vector<Convergents::Approximation>::iterator m=approx.end()-1;
      std::vector<Convergents::Approximation>::const_iterator m1, m2;
      m1=approx.end()-2;
      m2=approx.end()-3;
      
      // Calculate q_i.
      mpz_tdiv_qr(m->term.get_mpz_t(),
		  curN.get_mpz_t(),
		  curN.get_mpz_t(),
		  curD.get_mpz_t());
      mpz_swap(curD.get_mpz_t(), curN.get_mpz_t());
      
      // Calculate A_m and B_m.
      m->A=m->term*m1->A+m2->A;
      m->B=m->term*m1->B+m2->B;
    };
    
    return approx.size();
  }      

  // Find a maximum convergent that has an x value that is lower
  // than the parameter.
  mpq_class Convergents::seekClosest(mpz_class val) {
    assert(val>0);
    size_t idx=2;
    size_t last = approx.size()-1;
    // Find the first convergent that has a larger denominator than val.
    while (true) {
      if (idx>last) goto takePrevious;
      if (approx[idx].B>val) break;
      idx+=2;
    };
    approx.resize(idx+1);

    {
      // Check if there are intermediate convergents.
      mpz_class conv = val-approx[idx-2].B;
      mpz_fdiv_q(conv.get_mpz_t(), conv.get_mpz_t(),
		 approx[idx-1].B.get_mpz_t());
      if (sgn(conv)==0) goto takePrevious;

      assert(conv>0 && conv<approx[idx].term);
      approx[idx].term=conv;
      approx[idx].A=conv*approx[idx-1].A+approx[idx-2].A;
      approx[idx].B=conv*approx[idx-1].B+approx[idx-2].B;
      goto takeCurrent;
    };

    // There are none. Return the previous convergent (which is
    // smaller than val).
  takePrevious:
    idx-=2;
    approx.resize(idx+1);

  takeCurrent:
    assert(approx[idx].B<=val);
    a=approx[idx].A;
    b=approx[idx].B;
    // Ensure that we don't trigger a bug in the mpq_class assignment
    // operator.
    assert(sgn(b)>=0);
    return mpq_class(a, b);
  }

  // Given a and b, solve the equation alpha * a + gamma * b = 1.
  void solveDiophantine(const mpz_class& a, const mpz_class& b,
		        mpz_class& alpha, mpz_class& gamma) {
    // Continued fraction for the coefficients <a,b> of an inequality in the
    // original system. These are used to solve ax+by=1 for x and y.
    static Convergents original;
    alpha = 0; gamma = 0;
    size_t size = original.initialize(a, b);
    switch (size) {
      case 0 : if (sgn(a)<0) alpha=-1; else alpha=1; break;
      case 1 : if (sgn(b)<0) gamma=-1; else gamma=1; break;
      default : {
	alpha=original[size-2]->B;
	gamma=original[size-2]->A;
	// Davenport sais: |other.a|*alpha - |other.b|*gamma = 1 if
	// size-1 is even or ... = -1 if size-1 is odd. Adjust alpha
	// and gamma so that other.a*alpha + other.b*gamma = 1.
	if (size & 1) alpha=-alpha; else gamma=-gamma;
	if (sgn(a)<0) alpha=-alpha;
	if (sgn(b)<0) gamma=-gamma;
      }
    };

#ifdef DEBUG_CUT
    cerr << "orig convergents: " << original << endl;
    cerr << "coefficients a = " << a << ", b  = " << b
	 << ", alpha = " << alpha << ", gamma = " << gamma << endl;
#endif // DEBUG_CUT
    assert(a*alpha + b*gamma == 1);
  }

// Create a ray from the direction vector of an inequality. The
// second argument details towards which direction the ray should
// point.
Ray::Ray(const Inequality& e, RayDirection dir) {
  switch (dir) {
  case dirWithArrow: {
    dx=-e.b;
    dy=e.a;
  }; break;
  case dirToFeasible: {
    dx=-e.a;
    dy=-e.b;
  }; break;
  case dirAgainstArrow: {
    dx=e.b;
    dy=-e.a;
  }; break;
  };
};
      
void Ray::ensureLength(mpz_class l) {
  mpz_class factor = sqrt(dx*dx+dy*dy);
  mpz_cdiv_q(factor.get_mpz_t(), l.get_mpz_t(), factor.get_mpz_t());
  dx=dx*factor;
  dy=dy*factor;
}

// Create an arbitrary point that saturates the given inequality.
Point::Point(const Inequality& e) {
  if (sgn(e.a)==0) {
    // Bug in gmp: Assigning a fraction with a negative denominator
    // allocates a negative amount of memory. Workaround: assign the
    // fraction as integers and then canonicalize.
    y.get_num()=e.c;
    y.get_den()=e.b;
    y.canonicalize();
    return;
  };
  if (sgn(e.b)==0) {
    x.get_num()=e.c;
    x.get_den()=e.a;
    x.canonicalize();
    return;
  };
  mpq_class cDivB(e.c,e.b);
  cDivB.canonicalize();
  mpq_class cDivA(e.c,e.a);
  cDivA.canonicalize();
  if (cmp(cDivB, cDivA)<0) y=cDivB; else x=cDivA;
}

// Calculate the intersection point of the boundaries of the
// halfspaces induced by the inequalities.
Point::Point(const Inequality& e1,
	     const Inequality& e2) : x(0), y(0) {
  mpz_class detAB = det(e1.a,e1.b,e2.a,e2.b);
#ifdef DEBUG_INTERSECT
  cerr << "Intersection point of " << e1 << " and " << e2;
#endif // DEBUG_INTERSECT
  // Work around a bug in the assignment operator of the mpq_class
  // that can't handle negative denominators.
  mpq_class resX;
  resX.get_num()=det(e1.c,e1.b,e2.c,e2.b);
  resX.get_den()=detAB;
  resX.canonicalize();
  x=resX;
  mpq_class resY;
  resY.get_num()=det(e1.a,e1.c,e2.a,e2.c);
  resY.get_den()=detAB;
  resY.canonicalize();
  y=resY;
#ifdef DEBUG_INTERSECT
  cerr << " is <" << x << "," << y << ">" << endl; 
#endif // DEBUG_INTERSECT
}

int Point::ordering(const Point* p1, const Point* p2) {
  // In case two points are colinear with respect the the origin, we
  // need to ensure a total order in order not to confuse the Graham
  // scan. We chose to sort lexicographically be smaller x and larger
  // y.
  int res = calcDeterminant(p1, p2);
  if (res!=0) return res;
  res = sgn(p2->x-p1->x);
  if (res!=0) return res;
  return sgn(p1->y-p2->y);
}

// Compare two points for being counter-clockwise with regard to this
// point. Return 1 if the point p1, this, p2 are counter-clockwise.
int Point::calcBend(Point& p1, Point& p2) {
  return sgn(det(p2.x-x, p2.y-y, p1.x-x, p1.y-y));
}

// Increase the given value to the maximum of the |x| and |y|.
void Point::updateLength(mpz_class& len) const {
  mpz_class v;
  mpz_cdiv_q(v.get_mpz_t(), x.get_num().get_mpz_t(), x.get_den().get_mpz_t());
  v=abs(v);
  if (v>len) len=v;
  mpz_cdiv_q(v.get_mpz_t(), y.get_num().get_mpz_t(), y.get_den().get_mpz_t());
  v=abs(v);
  if (v>len) len=v;
}

// Check if this point is within the two given intervals. The
// direction argument gives a hint to the function which direction
// it needs to check.
template<bool isZ>
bool Point::inBox(const Interval<isZ>& xBound,
		  const Interval<isZ>& yBound,
		  Direction dir) const {
  switch (dir) {
  case east:  return xBound.upperIsInfinite() || x<=xBound.getUpper();
  case north: return yBound.upperIsInfinite() || y<=yBound.getUpper();
  case west:  return xBound.lowerIsInfinite() || x>=xBound.getLower();
  case south: return yBound.lowerIsInfinite() || y>=yBound.getLower();
  default: {
    if (xBound.upperIsFinite() && x>xBound.getUpper()) return false;
    if (yBound.upperIsFinite() && y>yBound.getUpper()) return false;
    if (xBound.lowerIsFinite() && x<xBound.getLower()) return false;
    if (yBound.lowerIsFinite() && y<yBound.getLower()) return false;
    return true;
  }
  }
}


// Check if this point within the square (-s,-s) (s,s)
bool Point::inBox(mpq_class s) const {
  return (-s<=x && x<=s && -s<=y && y<=s);
}

// A false from this function is a sufficient indicator that the line
// between p1 and p2 does not intersect with a square of size s.
bool outsideSquare(mpq_class s, Point& p1, Point& p2) {
  return ((p1.x>s && p2.x>s) ||
	  (p1.x<-s && p2.x<-s) ||
	  (p1.y>s && p2.y>s) ||
	  (p1.y<-s && p2.y<-s));
};


// Function that takes three rationals aa, bb and cc as arguments and
// divides these by the common factor of aa, bb (and cc if isZ is false).
Inequality::Inequality(const Point& p1, const Point& p2, bool isZ)
    throw (IllegalArgument) {
  mpq_class aa=p2.y-p1.y;
  mpq_class bb=p1.x-p2.x;
  mpq_class cc=aa*p1.x+bb*p1.y;
  if (sgn(aa)==0 && sgn(bb)==0) 
    throw new IllegalArgument("Inequality::Inequality(Point&, Point&)",
			      "points must be distinct");
  // calculate the least common multiple of the two coefficients
  mpz_class l;
  mpz_lcm(l.get_mpz_t(), aa.get_den_mpz_t(), bb.get_den_mpz_t());
  // calculate the greatest common divisor of the two coefficients
  mpz_class d;
  mpz_gcd(d.get_mpz_t(), aa.get_num_mpz_t(), bb.get_num_mpz_t());
  // ...and the constant, in case this is the rational domain
  if (!isZ) {
    mpz_lcm(l.get_mpz_t(), l.get_mpz_t(), cc.get_den_mpz_t());
    mpz_gcd(d.get_mpz_t(), d.get_mpz_t(), cc.get_num_mpz_t());
  };
  // the inequality with minimal coefficients can now be deduced by
  // multiplying the numerator with l and the denominator with d, yielding
  // a whole number for each rational
  mpz_mul(aa.get_num_mpz_t(), aa.get_num_mpz_t(), l.get_mpz_t());
  mpz_mul(aa.get_den_mpz_t(), aa.get_den_mpz_t(), d.get_mpz_t());
  mpz_divexact(a.get_mpz_t(), aa.get_num_mpz_t(), aa.get_den_mpz_t());
  mpz_mul(bb.get_num_mpz_t(), bb.get_num_mpz_t(), l.get_mpz_t());
  mpz_mul(bb.get_den_mpz_t(), bb.get_den_mpz_t(), d.get_mpz_t());
  mpz_divexact(b.get_mpz_t(), bb.get_num_mpz_t(), bb.get_den_mpz_t());
  mpz_mul(cc.get_num_mpz_t(), cc.get_num_mpz_t(), l.get_mpz_t());
  mpz_mul(cc.get_den_mpz_t(), cc.get_den_mpz_t(), d.get_mpz_t());
  if (isZ)
    mpz_fdiv_q(c.get_mpz_t(), cc.get_num_mpz_t(), cc.get_den_mpz_t());
  else 
    mpz_divexact(c.get_mpz_t(), cc.get_num_mpz_t(), cc.get_den_mpz_t());
}



// Copy an inequality.
Inequality::Inequality(const Inequality& other, bool reverse) {
  c=other.c;
  if (reverse) {
    a=other.b;
    b=other.a;
  } else {
    a=other.a;
    b=other.b;
  }
};


Inequality::Inequality(long aa, long bb, long cc) {
  assert(aa!=0 || bb!=0);
  a=mpz_class(aa);
  b=mpz_class(bb);
  c=mpz_class(cc);
}

// Create a new inequality by eliminating a common variable from the
// two given inequalities.
Inequality::Inequality(bool isZ,
		       const Inequality& e1, bool commonVar1,
		       const Inequality& e2, bool commonVar2) {
  assert(commonVar1 ? sgn(e1.b)!=0 : sgn(e1.a)!=0);
  assert(commonVar2 ? sgn(e2.b)!=0 : sgn(e2.a)!=0);
  assert((commonVar1 ? sgn(e1.b) : sgn(e1.a)) !=
   (commonVar2 ? sgn(e2.b) : sgn(e2.a)));
  mpz_class common1 = abs(commonVar1 ? e1.b : e1.a);
  mpz_class common2 = abs(commonVar2 ? e2.b : e2.a);
  a=(commonVar1 ? common2*e1.a : common2*e1.b);
  b=(commonVar2 ? common1*e2.a : common1*e2.b);
  c=common2*e1.c+common1*e2.c;
  assert (sgn(a)!=0 || sgn(b)!=0);
  mpz_class d;
  mpz_gcd(d.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
  if (!isZ) mpz_gcd(d.get_mpz_t(), d.get_mpz_t(), c.get_mpz_t());
  // Don't divide if the divisor is 1.
  if (mpz_cmpabs_ui(d.get_mpz_t(),1)==0) return;
  mpz_divexact(a.get_mpz_t(), a.get_mpz_t(), d.get_mpz_t());
  mpz_divexact(b.get_mpz_t(), b.get_mpz_t(), d.get_mpz_t());
  if (isZ)
    mpz_fdiv_q(c.get_mpz_t(), c.get_mpz_t(), d.get_mpz_t());
  else
    mpz_divexact(c.get_mpz_t(), c.get_mpz_t(), d.get_mpz_t());
}

// Create an inequality that saturates the point p and has its
// feasible space on the side given by dir.
Inequality::Inequality(const Point& p, FaceDirection dir) {
  switch (dir) {
  case dirUp : {
    b=-p.y.get_den();
    c=-p.y.get_num();
  }; break;
  case dirLeft : {
    a=p.x.get_den();
    c=p.x.get_num();
  }; break;
  case dirDown : {
    b=p.y.get_den();
    c=p.y.get_num();
  }; break;
  case dirRight : {
    a=-p.x.get_den();
    c=-p.x.get_num();
  }; break;
  };
};

template<>
bool Inequality::exactlyPi<false>(const Inequality& other) const {
  int diff=(other.calcDirection()-calcDirection()+4) % 4;
  if (diff!=2) return false;
  return (other.a*b==a*other.b);
}
  
template<>
bool Inequality::exactlyPi<true>(const Inequality& other) const {
  int diff=(other.calcDirection()-calcDirection()+4) % 4;
  if (diff!=2) return false;
  // Over Z, the coefficients have a gcd of 1. Hence they must be equal.
  return (mpz_cmpabs(a.get_mpz_t(), other.a.get_mpz_t())==0 &&
	  mpz_cmpabs(b.get_mpz_t(), other.b.get_mpz_t())==0);
}
  
// Check if the halfspace defined by this inequality is a superset of
// the intersection of the other two inequalities.
bool Inequality::includes(const Inequality& other1, 
			  const Inequality& other2) const {
#ifdef DEBUG_ENTAIL
  cerr << "Test {" << other1 << ", " << other2 << "} |= " << *this << endl;
#endif // DEBUG_ENTAIL
  mpz_class d=other1.a*other2.b-other2.a*other1.b;
#ifdef DEBUG_ENTAIL
  cerr << "d = " << d << endl;
#endif // DEBUG_ENTAIL
  if (sgn(d)==0) {
    // The other inequalities are parallel or anti-parallel. In both
    // cases we can check if either of the inequalities entails this
    // inequality.
    return includes(other1) || includes(other2);
  };
  mpq_class lambda1;
  // Do not use the constructor of mpq_class since it cannot handle
  // negative denominators.
  lambda1.get_num()=a*other2.b - other2.a*b;
  lambda1.get_den()=d;
  // this is necessary to turn -x/-y into x/y
  lambda1.canonicalize();
#ifdef DEBUG_ENTAIL
  cerr << "l1 = " << lambda1 << endl;
#endif // DEBUG_ENTAIL
  // lambda 1 is 0 iff other1 is parallel to this inequality
  if (sgn(lambda1)<0) return false;
  mpq_class lambda2;
  lambda2.get_num()=other1.a*b - a*other1.b;
  lambda2.get_den()=d;
  lambda2.canonicalize();
#ifdef DEBUG_ENTAIL
  cerr << "l2 = " << lambda2 << endl;
#endif // DEBUG_ENTAIL
  if (sgn(lambda2)<0) return false;
#ifdef DEBUG_ENTAIL
  cerr << "l1*c1+l2*c2 = " << mpq_class(lambda1*other1.c+lambda2*other2.c);
  cerr << ", c= " << c << endl;
#endif // DEBUG_ENTAIL
  return (cmp(mpq_class(lambda1*other1.c+lambda2*other2.c),c)<=0);
}

// Check if the halfspace defined by this inequality is a superset of
// the halfspace of the other inequality. The relation is non-strict.
bool Inequality::includes(const Inequality& other) const {
  // They must be parallel.
  if (sgn(other.a*b-a*other.b)!=0) return false;
  // They may not be antiparallel.
  if (sgn(other.a)*sgn(a)<0) return false;
  if (sgn(other.b)*sgn(b)<0) return false;
  //  Note that the direction of the less-than operator is reversed
  //  when the inequality is multiplied by a negative a.
  switch (sgn(other.a)) {
  case 1: return (a*other.c <= other.a*c); break;
  case -1: return (a*other.c >= other.a*c); break;
  default: {
    switch (sgn(other.b)) {
    case 1: return (b*other.c <= other.b*c); break;
    case -1: return (b*other.c >= other.b*c); break;
    default: return (sgn(other.c)<0 ||
		     (sgn(c)>=0 && sgn(a)==0 && sgn(b)==0)); break;
    }; break;
  };
  };
}

// Check if the halfspace defined by this inequality is a superset of
// the space defined by the intersection of the other inequality and
// the bounding box given by the two intervals on x and y.
template<bool isZ>
bool Inequality::includes(const Inequality& other,
			  const Interval<isZ>& xBound,
			  const Interval<isZ>& yBound) const {
  if (xBound.isEmpty() || yBound.isEmpty()) return true;

  // This test does not work if this and the other inequality are not
  // in the same quadrant. In case this and the other inequality lie
  // in a different quadrant inclusion could be checked with the box
  // alone . This is not implemented since it never happens within this
  // library.
  Direction dir=calcDirection();
  assert(dir==other.calcDirection());

  // If both inequalities are parallel, we can use the entailment
  // check between inequalities to determine if this inequalitiy
  // defines a superspace.
  int angle=this->calcDeterminant(other);
  if (angle==0) return this->includes(other);

  if (angle<0) {
    // The other inequality has a smaller angle than this inequality. 
    // Intersect the other inequality with the next bound in
    // counter-clockwise direction. This inequality is entailed
    // if it entails the intersection point.
    switch (dir) {
    case east: {
      if (yBound.upperIsFinite()) 
	return includes(Point(other.calculateX(yBound.getUpper()),
			      yBound.getUpper()));
      if (xBound.lowerIsFinite()) 
	return includes(Point(xBound.getLower(),
			      other.calculateY(xBound.getLower())));
    }; break;
    case north: {
      if (xBound.lowerIsFinite()) 
	return includes(Point(xBound.getLower(),
			      other.calculateY(xBound.getLower())));
      if (yBound.lowerIsFinite()) 
	return includes(Point(other.calculateX(yBound.getLower()),
			      yBound.getLower()));
    }; break;
    case west: {
      if (yBound.lowerIsFinite()) 
	return includes(Point(other.calculateX(yBound.getLower()),
			      yBound.getLower()));
      if (xBound.upperIsFinite()) 
	return includes(Point(xBound.getUpper(),
			      other.calculateY(xBound.getUpper())));
    }; break;
    case south: {
      if (xBound.upperIsFinite()) 
	return includes(Point(xBound.getUpper(),
			      other.calculateY(xBound.getUpper())));
      if (yBound.upperIsFinite()) 
	return includes(Point(other.calculateX(yBound.getUpper()),
			      yBound.getUpper()));
    }; break;
    default: assert(false); break;
    }
  } else {
    // This inequality has a smaller angle than the other inequality. 
    // Intersect the other inequality with the next bound in clockwise
    // direction. The other inequality is entailed if it entails the
    // intersection point.
    switch (dir) {
    case east: {
      if (xBound.upperIsFinite()) 
	return includes(Point(xBound.getUpper(),
			      other.calculateY(xBound.getUpper())));
      if (yBound.lowerIsFinite()) 
	return includes(Point(other.calculateX(yBound.getLower()),
			      yBound.getLower()));
    }; break;
    case north: {
      if (yBound.upperIsFinite()) 
	return includes(Point(other.calculateX(yBound.getUpper()),
			      yBound.getUpper()));
      if (xBound.upperIsFinite()) 
	return includes(Point(xBound.getUpper(),
			      other.calculateY(xBound.getUpper())));
    }; break;
    case west: {
      if (xBound.lowerIsFinite()) 
	return includes(Point(xBound.getLower(),
			      other.calculateY(xBound.getLower())));
      if (yBound.upperIsFinite()) 
	return includes(Point(other.calculateX(yBound.getUpper()),
			      yBound.getUpper()));
    }; break;
    case south: {
      if (yBound.lowerIsFinite()) 
	return includes(Point(other.calculateX(yBound.getLower()),
			      yBound.getLower()));
      if (xBound.lowerIsFinite()) 
	return includes(Point(xBound.getLower(),
			      other.calculateY(xBound.getLower())));
    }; break;
    default: assert(false); break;
    }
  };
  return false;
}

// Check if the halfspace defined by this inequality contains the
// bounding box.
template<>
bool Inequality::includes<false>(const Interval<false>& xBound,
				 const Interval<false>& yBound) const {
  // Make sure the inequality is not parallel to the axes.
  assert(a!=0);
  assert(b!=0);
  mpq_class x,y;
  if (sgn(a)>0) {
    if (xBound.upperIsInfinite()) return false;
    // Check against the upper bound of x.
    x=xBound.getUpper();
  } else {
    if (xBound.lowerIsInfinite()) return false;
    // Check against the lower bound of x.
    x=xBound.getLower();
  };
  if (sgn(b)>0) {
    if (yBound.upperIsInfinite()) return false;
    y=yBound.getUpper();
  } else {
    if (yBound.lowerIsInfinite()) return false;
    y=yBound.getLower();
  };
  return a*x+b*y<=c;
}      

template<>
bool Inequality::includes<true>(const Interval<true>& xBound,
				const Interval<true>& yBound) const {
  // Make sure the inequality is not parallel to the axes.
  //  assert(a!=0);
  //assert(b!=0);
  mpz_class x,y;
  if (sgn(a)>0) {
    if (xBound.upperIsInfinite()) return false;
    // Check against the upper bound of x.
    x=xBound.getUpper().get_num();
  } else {
    if (xBound.lowerIsInfinite()) return false;
    // Check against the lower bound of x.
    x=xBound.getLower().get_num();
  };
  if (sgn(b)>0) {
    if (yBound.upperIsInfinite()) return false;
    y=yBound.getUpper().get_num();
  } else {
    if (yBound.lowerIsInfinite()) return false;
    y=yBound.getLower().get_num();
  };
  return a*x+b*y<=c;
}      

template<bool isZ>
bool Inequality::tightenBounds(Interval<isZ>& xBound,
			       Interval<isZ>& yBound) const {
  bool changed = false;
  switch (calcDirection()) {
  case east:
    if (xBound.lowerIsFinite())
      if (yBound.updateUpper(calculateY(xBound.getLower()))) changed = true;
    break;
  case north:
    if (yBound.lowerIsFinite())
      if (xBound.updateLower(calculateX(yBound.getLower()))) changed = true;
    break;
  case west:
    if (xBound.upperIsFinite())
      if (yBound.updateLower(calculateY(xBound.getUpper()))) changed = true;
    break;
  case south:
    if (yBound.upperIsFinite())
      if (xBound.updateUpper(calculateX(yBound.getUpper()))) changed = true;
    break;
  default: break;
  };
  return changed;
}

int Inequality::compareWith(const Inequality& other) const {
  Direction dirThis = calcDirection();
  Direction dirOther = other.calcDirection();
  if (dirThis<dirOther)
    return (-1);
  else
    if (dirThis>dirOther)
      return 1;
    else
      // see the gmp docu for information on these casts
      return cmp(mpz_class(other.a*b), mpz_class(a*other.b));
}

// Calculate an integral cut between an inequality and a bound. The
// given direction determines if the bound represents the x or y
// direction and if it is the upper or lower value that is to be used.
Inequality* Inequality::calculateCut(const Interval<true>& xBound,
				     const Interval<true>& yBound,
				     Direction dir) const {
  if (xBound.isEmpty()) return NULL;
  if (yBound.isEmpty()) return NULL;
  switch (dir) {
  case east:
    return (xBound.upperIsInfinite() ? NULL :
	    calculateCut(Inequality(1, 0, xBound.getUpper().get_num())));
  case north:
    return (yBound.upperIsInfinite() ? NULL :
	    calculateCut(Inequality(0, 1, yBound.getUpper().get_num())));
  case west:
    return (xBound.lowerIsInfinite() ? NULL :
	    calculateCut(Inequality(-1, 0, -xBound.getLower().get_num())));
  case south:
    return (yBound.lowerIsInfinite() ? NULL :
	    calculateCut(Inequality(0, -1, -yBound.getLower().get_num())));
  default: {
    assert(false);
    return NULL;
  }
  }
}

  // Continued fraction for the coefficients of an inequality in the
  // rotated system, where the second inequality is an upper x bound.
  Convergents Inequality::rotated;

  // Calculate an integral cut. This inequality has to have a smaller
  // angle than the other. The generated cut will have an angle
  // inbetween the two given. The next cut can be calculated by
  // calling this function again on the last cut with the other
  // inequality as argument. The function returns NULL if there are no
  // more cuts.
  Inequality* Inequality::calculateCut(const Inequality& other) const {
    assert(lessThanPi(other));

    // There is no cut if the intersection point of the two given
    // inequalities is integral.
    {
      Point p = Point(*this, other);
      if (p.isIntegral()) return NULL;
#ifdef DEBUG_CUT
      cerr << "intersection point " << p << " between " << *this
	   << " and " << other << " not integral." << endl;
#endif // DEBUG_CUT
    };

    // Ensure that the gcd of the other inequality is 1. This should always
    // be the case but it is important since the calculation below otherwises
    // doesn't work.
#ifndef NDEBUG
    mpz_class g;
    mpz_gcd(g.get_mpz_t(), other.a.get_mpz_t(), other.b.get_mpz_t());
    assert(g==1);
#endif
    
    mpz_class alpha, beta, gamma, delta;
    {
      // Calculate the unimodular rotation matrix.
      mpz_class alpha0, gamma0;
      solveDiophantine(other.a, other.b, alpha0, gamma0);
      
      // Calculate k.
      mpz_class Delta = det(a, b, other.a, other.b);
      assert(Delta>0);
      mpz_class k = - (a*alpha0) - (b*gamma0);
      mpz_fdiv_q(k.get_mpz_t(), k.get_mpz_t(), Delta.get_mpz_t());

#ifdef DEBUG_CUT
      cerr << "k = " << a << "*alpha0 + " << b << "*gamma0 / " << Delta
	   <<" = " << k << endl;
#endif // DEBUG_CUT
      
      // Set the matrix coefficients.
      alpha = alpha0 + k*other.b;
      beta = other.b;
      gamma = gamma0 - k*other.a;
      delta = -other.a;
      assert( det(alpha, beta, gamma, delta) == 1 ||
	      det(alpha, beta, gamma, delta) == -1);
    };

#ifdef DEBUG_CUT
    cerr << "matrix is:\t" << alpha << "\t" << beta << endl;
    cerr << "          \t" << gamma << "\t" << delta << endl;
#endif // DEBUG_CUT

    // Calculate the rotated system.
    mpz_class t = a*alpha + b*gamma;
    mpz_class u = a*beta + b*delta;
    assert(   1 == other.a*alpha + other.b*gamma);
    assert(   0 == other.a*beta + other.b*delta);

#ifdef DEBUG_CUT
    cerr << "coeff are: t=" << t << "\tu=" << u << endl;
#endif // DEBUG_CUT
    
    assert(t<=0);
    assert(u>=0);

    // Calculate an integral point <x1, y1> on t x1 + u y1 = c
    // that also satisfies x1 <= other.c.
    mpz_class x1;
    mpz_class y1;
    {
      size_t size = rotated.initialize(t,u);
      // Calculate an integral point <x0,y0> on t x0 + u y0 = c.
      mpz_class x0;
      mpz_class y0;
      switch (size) {
      case 0 : if (sgn(t)<0) x0=c; else x0=-c; break;
      case 1 : if (sgn(u)<0) y0=c; else y0=-c; break;
      default : {
	y0 = rotated[size-2]->A*c;
	x0 = rotated[size-2]->B*c;
	// At this point: - t x0 - u y0 = c if size-1 is odd and - t x0
	// - u y0 = -c if last is even.
	if (!(size & 1)) { x0=-x0; y0=-y0; }
      }
      };

#ifdef DEBUG_CUT
      cerr << "rotated: " << rotated << endl;
      cerr << "x0 = " << x0 << ", y0 = " << y0 << endl;
#endif // DEBUG_CUT
      assert(t*x0+u*y0==c);

      // Move <x0,y0> along the equality t x0 + u y0 = c such that
      // x0 <= other.c.
      mpz_class ofs = other.c-x0;
      mpz_fdiv_q(ofs.get_mpz_t(), ofs.get_mpz_t(), u.get_mpz_t());
      x1=x0+(ofs*u);
      y1=y0-(ofs*t);
    };

#ifdef DEBUG_CUT
    cerr << "x_1 = " << x1 << ", y_1 = " << y1 << endl;
#endif // DEBUG_CUT

    mpq_class slope = rotated.seekClosest(other.c-x1);

#ifdef DEBUG_CUT
    cerr << "slope is " << slope << endl;
#endif // DEBUG_CUT
    
    Inequality* res =
      new Inequality(slope.get_num()*delta + slope.get_den()*gamma,
		     -(slope.get_num()*beta + slope.get_den()*alpha),
		     slope.get_den()*y1 - slope.get_num()*x1);

#ifdef DEBUG_CUT
    cerr << "cut between " << *this << " and " << other << " is " 
	 << *res << endl;
#endif // DEBUG_CUT
    return res;
  }

  defMemDbg(,Inequality,e,E)


// Calculate an extremal integral value of x in the intersection space
// of the two given inequalities. The inequalities must be increasing
// in angle. The Boolean flag indicates if an upper or a lower bound
// is to be calculated, that is, if isUpper is true the result is rounded
// down (!). The flag has no effect in the rational domain.
template<bool isZ>
mpq_class extremeX(bool isUpper, 
		   const Inequality& e1,
		   const Inequality& e2) {
  assert(e1.lessThanPi(e2));
  Point p = Point(e1, e2);
#ifdef DEBUG_EXTREME
  cerr << "extreme " << (isUpper ? "upper" : "lower")
       << " x value of point " << p;
#endif // DEBUG_EXTREME
  if (isZ) {
    if (isUpper)
      mpz_fdiv_q(p.getX().get_num().get_mpz_t(), p.getX().get_num().get_mpz_t(),
        p.getX().get_den().get_mpz_t());
    else
      mpz_cdiv_q(p.getX().get_num().get_mpz_t(), p.getX().get_num().get_mpz_t(),
        p.getX().get_den().get_mpz_t());
    p.getX().get_den()=1;
  };
#ifdef DEBUG_EXTREME
  cerr << " is " << p.getX() << endl;
#endif // DEBUG_EXTREME
  return p.getX();
}

template<bool isZ>
mpq_class extremeY(bool isUpper, 
		   const Inequality& e1,
		   const Inequality& e2) {
  assert(e1.lessThanPi(e2));
  Point p(e1, e2);
#ifdef DEBUG_EXTREME
  cerr << "extreme " << (isUpper ? "upper" : "lower")
       << " y value of point " << p;
#endif // DEBUG_EXTREME
  if (isZ) {
    if (isUpper)
      mpz_fdiv_q(p.getY().get_num().get_mpz_t(), p.getY().get_num().get_mpz_t(),
        p.getY().get_den().get_mpz_t());
    else
      mpz_cdiv_q(p.getY().get_num().get_mpz_t(), p.getY().get_num().get_mpz_t(),
        p.getY().get_den().get_mpz_t());
    p.getY().get_den()=1;
  };
#ifdef DEBUG_EXTREME
  cerr << " is " << p.getY() << endl;
#endif // DEBUG_EXTREME
  return p.getY();
};

// 
// Function to tighten bounds to the closest integral point according
// to the two inequalities.
template<>
void refineBoundsToZ<true>(Interval<true>& xBound,
			   Interval<true>& yBound,
			   const Inequality* currentF,
			   const Inequality* nextF) {
  Direction cDir = currentF->calcDirection();
  Direction nDir = nextF->calcDirection();

#ifdef DEBUG_REFINE
  cerr << "checking integral bound between (" << cDir << ")" << *currentF
       << " and (" << nDir << ")" << *nextF << endl;
#endif // DEBUG_REFINE

  // Find out how the two inequalities intersect with the bounds by calculating
  // an upper and a lower value at which the inequalities intersect with the bound
  // in direction dir.
  mpq_class lower, upper;
  switch (cDir) {
  case east: {
    // We assume that all bounds are rationally as tight as possible. Thus,
    // if the next bound is infinite, we must be looking at a ray.
    // There will be an integral point, so return immediately.
    if (yBound.upperIsInfinite()) return;
    upper = currentF->calculateX(yBound.getUpper());
    lower = nextF->calculateX(yBound.getUpper());
    // If the next inequality intersects with the counter-clockwise next bound,
    // two cases can arise. Firstly, both inequality intersect with different
    // bounds. In this case the corner where the two bounds intersect constitutes
    // an integral point which both inequalities entail. We may return.
    if (xBound.lowerIsFinite() && lower<xBound.getLower()) {
      if (upper>=xBound.getLower()) return;
      // In the second case, both inequalities intersect with the next bound.
      // Calculate the axis intersections there.
      upper = currentF->calculateY(xBound.getLower());
      lower = nextF->calculateY(xBound.getLower());
    }
  }; break;
  case north: {
    if (xBound.lowerIsInfinite()) return;
    upper = currentF->calculateY(xBound.getLower());
    lower = nextF->calculateY(xBound.getLower());
    if (yBound.lowerIsFinite() && lower<yBound.getLower()) {
      if (upper>=yBound.getLower()) return;
      lower = currentF->calculateX(yBound.getLower());
      upper = nextF->calculateX(yBound.getLower());
    }
  }; break;
  case west: {
    if (yBound.lowerIsInfinite()) return;
    lower = currentF->calculateX(yBound.getLower());
    upper = nextF->calculateX(yBound.getLower());
    if (xBound.upperIsFinite() && upper>xBound.getUpper()) {
      if (lower<=xBound.getUpper()) return;
      lower = currentF->calculateY(xBound.getUpper());
      upper = nextF->calculateY(xBound.getUpper());
    }
  }; break;
  case south: {
    if (xBound.upperIsInfinite()) return;
    lower = currentF->calculateY(xBound.getUpper());
    upper = nextF->calculateY(xBound.getUpper());
    if (yBound.upperIsFinite() && upper>yBound.getUpper()) {
      if (lower<=yBound.getUpper()) return;
      upper = currentF->calculateX(yBound.getUpper());
      lower = nextF->calculateX(yBound.getUpper());
    }
  }; break;
  default: assert(0); break;
  };

#ifdef DEBUG_REFINE
  cerr << "intersection with bound is between " << lower << " and " << upper;
#endif // DEBUG_REFINE

  // Given an upper and a lower intersection with one of the bounds, round these two
  // values towards each other. If they have their integral value are still in the
  // same order, then there is an integral point between them. Since a feasible
  // integral point exists on the bound, there is no need to search for one.
  mpz_fdiv_q(upper.get_num().get_mpz_t(),
	     upper.get_num().get_mpz_t(), upper.get_den().get_mpz_t());
  mpz_cdiv_q(lower.get_num().get_mpz_t(),
	     lower.get_num().get_mpz_t(), lower.get_den().get_mpz_t());

#ifdef DEBUG_REFINE
  cerr << ", rounded: " << lower.get_num() << " and " << upper.get_num() << endl;
#endif // DEBUG_REFINE

  if (upper.get_num() >= lower.get_num()) return;

  // Calculate cuts from the current inequality to the bound in the direction
  // of quadrant that currentF leads up to.
  Direction dir = (cDir == south ? east : (Direction) (((int) cDir)+1));
  Inequality* newCut = currentF->calculateCut(xBound, yBound, dir);
  // Since there is no feasible, integral point on the bound, there must be at
  // least one cut.
  assert(newCut);
  Inequality* cut = 0;
  do {
#ifdef DEBUG_REFINE
    cerr << "cut with bound from first inequality: " << *newCut << endl;
#endif // DEBUG_REFINE
    delete cut;
    cut = newCut;
    newCut = cut->calculateCut(xBound, yBound, dir);
  } while (newCut && !newCut->includes(*cut, *nextF));
  delete newCut;

  // At this point, cut contains an inequality with which we can calculate new
  // cuts with the next inequality.
  newCut = cut;
  cut = 0;
  while (newCut && newCut->calcDirection()==cDir) {
    delete cut;
    cut = newCut;
    newCut = cut->calculateCut(*nextF);

#ifdef DEBUG_REFINE
    if (newCut) cerr << "cut towards second inequality: ("
		     << newCut->calcDirection() << ")" << *newCut << endl;
#endif // DEBUG_REFINE
  };
  const Inequality* last = newCut;
  if (!last) last=nextF;
#ifdef DEBUG_REFINE
  cerr << "updating bound next to first inequality (dir=" << dir << ")" << endl;
#endif // DEBUG_REFINE
  switch (dir) {
  case east: xBound.updateUpper(extremeX<true>(true, *cut, *last)); break;
  case north: yBound.updateUpper(extremeY<true>(true, *cut, *last)); break;
  case west: xBound.updateLower(extremeX<true>(false, *cut, *last)); break;
  case south: yBound.updateLower(extremeY<true>(false, *cut, *last)); break;
  default: assert(0); break;
  };

  // Continue calculating cuts until we cross the direction to nextF.
  dir = (nDir == east ? south :  (Direction) (((int) nDir)-1));
  while (newCut && newCut->calcDirection()==dir) {
    delete cut;
    cut = newCut;
    newCut = cut->calculateCut(*nextF);

#ifdef DEBUG_REFINE
    if (newCut) cerr << "more cuts towards second inequality: " << *newCut << endl;
#endif // DEBUG_REFINE
  };
  last = newCut;
  if (!last) last=nextF;
  switch (nDir) {
  case east: xBound.updateUpper(extremeX<true>(true, *cut, *last)); break;
  case north: yBound.updateUpper(extremeY<true>(true, *cut, *last)); break;
  case west: xBound.updateLower(extremeX<true>(false, *cut, *last)); break;
  case south: yBound.updateLower(extremeY<true>(false, *cut, *last)); break;
  default: assert(0); break;
  };
  delete newCut;
  delete cut;
};

template<>
void refineBoundsToZ<false>(Interval<false>& xBound,
			    Interval<false>& yBound,
			    const Inequality* currentF,
			    const Inequality* nextF) {};

// Define both combinations of all template functions.

template bool Point::inBox<false>(const Interval<false>& xBound,
				  const Interval<false>& yBound,
				  Direction dir=total) const;
template bool Inequality::exactlyPi<false>(const Inequality& other) const;
template bool Inequality::areEquality<false>(const Inequality& other) const;
template bool Inequality::includes<false>(const Inequality& other,
					  const Interval<false>& xBound,
					  const Interval<false>& yBound) const;
template bool Inequality::includes<false>(const Interval<false>& xBound,
					  const Interval<false>& yBound) const;
template bool Inequality::tightenBounds(Interval<false>& xBound,
					Interval<false>& yBound) const;
template mpq_class extremeX<false>(bool isUpper,
				   const Inequality& e1, const Inequality& e2);
template mpq_class extremeY<false>(bool isUpper,
				   const Inequality& e1, const Inequality& e2);

template bool Point::inBox<true>(const Interval<true>& xBound,
				 const Interval<true>& yBound,
				 Direction dir=total) const;
template bool Inequality::exactlyPi<true>(const Inequality& other) const;
template bool Inequality::areEquality<true>(const Inequality& other) const;
template bool Inequality::includes<true>(const Interval<true>& xBound,
					 const Interval<true>& yBound) const;
template bool Inequality::includes<true>(const Inequality& other,
					 const Interval<true>& xBound,
					 const Interval<true>& yBound) const;
template bool Inequality::tightenBounds(Interval<true>& xBound,
					Interval<true>& yBound) const;
template mpq_class extremeX<true>(bool isUpper,
				  const Inequality& e1, const Inequality& e2);
template mpq_class extremeY<true>(bool isUpper,
				  const Inequality& e1, const Inequality& e2);



}; // namespace Tvpi
