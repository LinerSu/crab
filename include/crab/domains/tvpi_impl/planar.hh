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
/* planar.cpp
 * Operations on planar polyhedra.
 */

#ifndef __PLANAR_H
#define __PLANAR_H

#include <iostream>
#include <gmpxx.h>
#include "interval.hh"
#include "tvpiexception.hh"
#include "memory.hh"
#include "common.hh"

// Forward declaration of classes which are given access rights in Ray
// and Point.
namespace Tvpi {
  class Point;
  class Inequality;
};

namespace Tvpi {

  // A class to store the development of a continued fraction.
  class Convergents {
    
    mpz_class a, b;

  public:
    class Approximation {
      friend class Convergents;
    private:
      mpz_class term;

    public:
      mpz_class A;
      mpz_class B;

      Approximation() {};

      friend std::ostream& operator<<(std::ostream& stream,
				      const Approximation& a) {
	return stream << a.term << " (" << a.A << "/" << a.B << ")";
      }
    };

  private:
    // Convergents and terms calculated so far.
    std::vector<Approximation> approx;

  public:
    // Set the coefficients that should be approximated. The
    // coefficients are stipped off their sign. This is a no-op of the
    // coefficients are the same as they were before the call. The
    // return value denotes the number of convergents that were
    // calculated.
    size_t initialize(mpz_class aa, mpz_class bb);

    // Ask for the nth approximation. Returns NULL if a smaller n
    // would already be exact.
    inline const Approximation* operator[](size_t n) const {
      if (n<approx.size()) return &approx[n]; else return NULL;
    };

    // Find a maximum convergent that has an x value that is lower
    // than the parameter.
    mpq_class seekClosest(mpz_class val);

    friend std::ostream& operator<<(std::ostream& stream,
				    const Convergents& c) {
      stream << "approx " << c.a << "/" << c.b << ": ";
      for (size_t i=0; i<c.approx.size(); i++) {
	if (i) stream << ", ";
	stream << *c[i];
      };
      return stream;
    }

  };

// Given a and b, solve the equation alpha * a + gamma * b = 1.
void solveDiophantine(const mpz_class& a, const mpz_class& b,
		      mpz_class& alpha, mpz_class& gamma);

// A ray is a tuple indicating a direction. It is represented by two integers.
class Ray {
  friend class Inequality;
  friend class Point;
  mpz_class dx;
  mpz_class dy;
 public:
  Ray(mpz_class dxx, mpz_class dyy): dx(dxx), dy(dyy) {};
  Ray() {};

  // When creating a ray from the coefficients of an inequality, the
  // direction vector needs to be rotated by pi/2 (dirWithArrow), pi
  // (dirToFeasible) or 3/2 pi (dirAgainstArrow).
  enum  RayDirection { dirWithArrow, dirToFeasible, dirAgainstArrow };

  // Create a ray from the direction vector of an inequality. The
  // second argument details towards which direction the ray should
  // point.
  Ray(const Inequality& e, RayDirection dir);
      
  // Swap the x and y coordinates around.
  void swap() {
    mpz_swap(dx.get_mpz_t(), dy.get_mpz_t());
  };

  Ray& operator=(const Ray& r) {
    dx=r.dx;
    dy=r.dy;
    return *this;
  };

  // Check whether two rays oppose each other.  This test possibly
  // returns false if the rays are opposed but do not have the same
  // axis intercept.
  bool opposes(Ray& other) { 
    return (sgn(dx+other.dx)==0) && (sgn(dy+other.dy)==0);
  };

  // Ensure that the ray is at least l long
  void ensureLength(mpz_class l);

  friend std::ostream& operator<<(std::ostream& stream,
				  const Ray& r) {
    return stream << "[" << r.dx << "," << r.dy << "]";
  };
};

// A macro to calculate the determinant of a 2x2 matrix
#define det(a11,a12,a21,a22) (a11)*(a22)-(a12)*(a21)

// This is a single point in planar space. It is defined by two rationals.
class Point {
  friend class Inequality;
  mpq_class x;
  mpq_class y;
 public:
  inline Point() {};

  Point(mpq_class xx, mpq_class yy) : x(xx), y(yy) {};

  Point(const long xx, const long yy) {
    x = mpq_class(xx);
    y = mpq_class(yy);
  };
  // Create a new point by translating another by a ray. It is assumed that
  // the ray has the appropriate length.
  Point(const Ray& r, const Point& p) {
    x=p.x+r.dx;
    y=p.y+r.dy;
  };
  
  // Create an arbitrary point that saturates the given inequality.
  Point(const Inequality& e);

  // Intersect two inequalities. Note that the angle between these
  // inequalities may not be pi or 2pi.
  Point(const Inequality& e1, const Inequality& e2);

  // Compare two points for being counter-clockwise with respect to
  // the origin. This boils down to calculating the determinant of the
  // coefficients.
  static int ordering(const Point* p1, const Point* p2);

  // Check if two points are counter-clockwise (1), colinear (0) or
  // clockwise (-1) with respect to the origin.
  static inline int calcDeterminant(const Point* p1, const Point* p2) {
    return sgn(det(p2->x, p2->y, p1->x, p1->y));
  };

  // Check if this point has a smaller x coordiante.
  static int smallerX(const Tvpi::Point* p1, const Tvpi::Point* p2) {
    return cmp(p1->x, p2->x)>0;
  };

  // Check if this point has a smaller x coordiante.
  static int largerX(const Tvpi::Point* p1, const Tvpi::Point* p2) {
    return cmp(p1->x, p2->x)<0;
  };

  // Query the x value.
  mpq_class const& getX() const { return x; };
  mpq_class& getX() { return x; };

  // Query the y value.
  mpq_class const& getY() const { return y; };
  mpq_class& getY() { return y; };

  // Compare two points for being counter-clockwise with regard to
  // this point.
  int calcBend(Point& p1, Point& p2);

  // Increase the given value to the maximum of the |x| and |y|.
  void updateLength(mpz_class& len) const;

  // Calculate the (square of the) absolute distance from here to the
  // other point.
  mpq_class calcDistance(Point& p) const {
    mpq_class dx=abs(x-p.x);
    mpq_class dy=abs(y-p.y);
    return (dx*dx+dy*dy);
  };

  // Check if this point is within the two given intervals. The
  // direction argument gives a hint to the function which direction
  // it needs to check.
  template<bool isZ>
  bool inBox(const Interval<isZ>& xBound,
	     const Interval<isZ>& yBound,
	     Direction dir=total) const;

  // Check if this point within the square (-s,-s) (s,s)
  bool inBox(mpq_class s) const;

  // A true from this function is a sufficient indicator that the
  // line between p1 and p2 does not intersect with a square of size
  // s. Instead of testing wether the line between p1 and p2
  // intersects with the square of size s (as suggested in the convex
  // hull paper), we just check if the x (or y) coordinates are both
  // greater than s or both smaller than s. This test is much cheaper.
  friend bool outsideSquare(mpq_class s, Point& p1, Point& p2);
  
  // Check whether this point lies in Z.
  bool isIntegral() const {
    return (x.get_den()==1 && y.get_den()==1);
  };

  // Swap the x and y coordinates around.
  void swap() {
    mpq_swap(x.get_mpq_t(), y.get_mpq_t());
  };

  // Swap the content of this point with the other point.
  void swapWith(Point& other) {
      mpq_swap(x.get_mpq_t(), other.x.get_mpq_t());
      mpq_swap(y.get_mpq_t(), other.y.get_mpq_t());
  };

  bool operator!=(const Point& other) const {
    return (cmp(this->x, other.x)!=0 || cmp(this->y, other.y)!=0);
  };

  bool operator==(const Point& other) const {
    return (cmp(this->x, other.x)==0 && cmp(this->y, other.y)==0);
  };

  // Yield true if this point has a higher x coordinate than the
  // other. If both have the same, return true if the y coordinate is
  // lower.
  bool operator<(const Point& other) const {
    int res=cmp(x, other.x);
    if (res==0) return (cmp(y,other.y)<0);
    return (res>0);
  };

  Point& operator-=(const Point& other) {
    x-=other.x;
    y-=other.y;
    return *this;
  };

  Point& operator+=(const Point& other) {
    x+=other.x;
    y+=other.y;
    return *this;
  };
  
  friend std::ostream& operator<<(std::ostream& stream,
				  const Point& p) {
    return stream << "<" << p.x << "," << p.y << ">";
  };

};

// An inequality ax+by<=c defines a half-space which entails all point
// <x,y> that satisfy the inequality.
class Inequality {
  friend class Point;
  friend class Ray;

  mpz_class a;
  mpz_class b;
  mpz_class c;
  
  // The continued fraction of the coefficients of an inequality, used
  // in calculateCut.
  static Convergents rotated;

 public:
  // Constructor that takes the internal types as arguments.
  Inequality(mpz_class aa, mpz_class bb, mpz_class cc) : a(aa), b(bb), c(cc) {};


  inline void canonicalize(bool isZ) {
    mpz_class gcd;
    mpz_gcd(gcd.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    if (isZ) {
      if (gcd>1) {
	mpz_divexact(a.get_mpz_t(), a.get_mpz_t(), gcd.get_mpz_t());
	mpz_divexact(b.get_mpz_t(), b.get_mpz_t(), gcd.get_mpz_t());
	mpz_fdiv_q(c.get_mpz_t(), c.get_mpz_t(), gcd.get_mpz_t());
      }
    } else {
      mpz_gcd(gcd.get_mpz_t(), gcd.get_mpz_t(), c.get_mpz_t());
      if (gcd>1) {
	mpz_divexact(a.get_mpz_t(), a.get_mpz_t(), gcd.get_mpz_t());
	mpz_divexact(b.get_mpz_t(), b.get_mpz_t(), gcd.get_mpz_t());
	mpz_divexact(c.get_mpz_t(), c.get_mpz_t(), gcd.get_mpz_t());
      }
    }
  };

  // Default constructor. Creates an illegal inequality with all
  // coefficients zero.
  Inequality() {};

  // Copy an inequality. If the optional boolean flag is set, the two
  // coefficients are exchanged in the resulting inequality.
  Inequality(const Inequality& other, bool reverse=false);

  // Create a new inequality from two points. If p1, p2 and p3 are ordered
  // counter-clockwise, then p3 lies in the feasible space of the halfspace
  // defined by this inequality.
  Inequality(const Point& p1, const Point& p2, bool isZ=false)
    throw (IllegalArgument);

  // Create an inequality from the three coefficients a, b and c in
  // ax+by <= c. For testing only.
  Inequality(long a, long b, long c);

  // Create a new inequality by eliminating a common variable from the
  // two given inequalities. The two boolean variables determine
  // whether the first (false) or the second variable (true) is the
  // common variable. The the first variable of the result corresponds
  // to the non-common variable of the first inequality, the second
  // variable of the result corresponds to the non-common variable of
  // the second inequality.
  Inequality(bool isZ,
	     const Inequality& e1, bool commonVar1,
	     const Inequality& e2, bool commonVar2);

  // Directions of inequalities that are 0, pi/2, pi and 3pi/2
  // repectively.
  enum FaceDirection { dirUp, dirLeft, dirDown, dirRight };

  // Create an inequality that saturates the point p and has its
  // feasible space on the side given by dir.
  Inequality(const Point& p, FaceDirection dir);

  // Test if a point saturates this inequality.
  inline bool isSaturated(const Point& p) const {
    return (a*p.x+b*p.y==c);
  };
  
  // Test if a point is feasible in this inequality.
  inline bool includes(const Point& p) const {
    return (a*p.x+b*p.y<=c);
  };

  // Determine if the angle between this and the given inequality is
  // less than 180 degrees. Assume that the other inequality always
  // has a larger angle, that is, if two inequalities of the same
  // coefficients are compare, assume that the angle is 2pi and return
  // false.
  inline bool lessThanPi(const Inequality& other) const {
    int diff=(other.calcDirection()-calcDirection()+4) % 4;
    if (diff==0 || diff==2) return (other.a*b<a*other.b);
    return (diff<2);
  };
  
  // Test if two inequalities are exactly 180 degree apart.
  template<bool isZ>
  bool exactlyPi(const Inequality& other) const;

  // Return the sign of the determinant. The result is positive if the
  // other inequality has a greater angle than this inequality, assuming
  // they lie in the same quadrant.
  inline int calcDeterminant(const Inequality& other) const {
    return sgn(a*other.b-other.a*b);
  };

  // Check if two inequalities constitute an equality.
  template<bool isZ>
  inline bool areEquality(const Inequality& other) const {
    return exactlyPi<isZ>(other) && c==-other.c;
  };

  // Check that both coefficients are non-zero.
  bool isTVPI() { return sgn(a)!=0 && sgn(b)!=0; };

  // Check if the halfspace defined by this inequality is a superset
  // of the intersection of the other two inequalities. It is assumed
  // that the other inequalities define a satisfiable system and that
  // they are non-redundant. The relation is interpreted non-strictly.
  bool includes(const Inequality& other1, const Inequality& other2) const;
  
  // Check if the halfspace defined by this inequality is a superset of
  // the halfspace of the other inequality. The relation is non-strict.
  bool includes(const Inequality& other) const;

  // Check if the halfspace defined by this inequality is a superset
  // of the space defined by the intersection of the other inequality
  // and the bounding box given by the two intervals on x and y.
  template<bool isZ>
  bool includes(const Inequality& other,
		const Interval<isZ>& xBound,
		const Interval<isZ>& yBound) const;

  // Check if the halfspace defined by this inequality contains the
  // bounding box.
  template<bool isZ>
  bool includes(const Interval<isZ>& xBound,
		const Interval<isZ>& yBound) const;

  // Tighten the next bound in counter-clockwise direction with
  // respect to this inequality.
  template<bool isZ>
  bool tightenBounds(Interval<isZ>& xBound,
		     Interval<isZ>& yBound) const;

  // Swap the coefficients of the two variables.
  inline void swapVars() {
    mpz_swap(a.get_mpz_t(), b.get_mpz_t());
  };

  // Check if an inequality is a tautology (both coefficients are zero
  // and the constant is greater than 0). The constant is actually not
  // checked.
  inline bool isTautology() {
    return (a==0 && b==0);
  };

  // Stretch this inequality in x direction.
  inline void stretchX(Mult m) {
    mpz_mul_2exp(a.get_mpz_t(), a.get_mpz_t(), (unsigned long int) m);
    mpz_mul_2exp(c.get_mpz_t(), c.get_mpz_t(), (unsigned long int) m);
    mpz_class div;
    mpz_gcd(div.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    if (div!=1) {
      assert(mpz_divisible_p(c.get_mpz_t(), div.get_mpz_t()));
      mpz_divexact(a.get_mpz_t(), a.get_mpz_t(), div.get_mpz_t());
      mpz_divexact(b.get_mpz_t(), b.get_mpz_t(), div.get_mpz_t());
      mpz_divexact(c.get_mpz_t(), c.get_mpz_t(), div.get_mpz_t());
    }
  };

  // Stretch this inequality in y direction.
  inline void stretchY(Mult m) {
    mpz_mul_2exp(b.get_mpz_t(), b.get_mpz_t(), (unsigned long int) m);
    mpz_mul_2exp(c.get_mpz_t(), c.get_mpz_t(), (unsigned long int) m);
    mpz_class div;
    mpz_gcd(div.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    if (div!=1) {
      assert(mpz_divisible_p(c.get_mpz_t(), div.get_mpz_t()));
      mpz_divexact(a.get_mpz_t(), a.get_mpz_t(), div.get_mpz_t());
      mpz_divexact(b.get_mpz_t(), b.get_mpz_t(), div.get_mpz_t());
      mpz_divexact(c.get_mpz_t(), c.get_mpz_t(), div.get_mpz_t());
    }
  };

  // Access the first coefficient of this inequality.
  mpz_class& getA() { return a; }
  mpz_class const& getA() const { return a; }

  // Access the second coefficient of this inequality.
  mpz_class& getB() { return b; }
  mpz_class const& getB() const { return b; }

  // Access the constant of this inequality.
  mpz_class& getC() { return c; }
  mpz_class const& getC() const { return c; }


  // Calculate the intersection with this bound on x.
  mpq_class calculateY(mpq_class x) const {
    // Work around a bug in the assignment operator of the mpq_class
    // that can't handle negative denominators.
    mpq_class y;
    y.get_num()=c-a*x;
    y.get_den()=b;
    y.canonicalize();
    return y;
  };

  // Calculate the intersection with this bound on y.
  mpq_class calculateX(mpq_class y) const {
    mpq_class x;
    x.get_num()=c-b*y;
    x.get_den()=a;
    x.canonicalize();
    return x;
  };

  //  Compare the angle of two inequalities. Returns a negative number
  //  for smaller angle, 0 for equal angle and a positive number for a
  //  larger angle.
  int compareWith(const Inequality& other) const;

  static int compare(const Inequality** eq1, const Inequality** eq2) {
    return (*eq1)->compareWith(**eq2);
  };

  // Calculate an integral cut between an inequality and a bound. The
  // given direction determines if the bound represents the x or y
  // direction and if it is the upper or lower value that is to be
  // used.
  Inequality* calculateCut(const Interval<true>& xBound,
			   const Interval<true>& yBound,
			   Direction dir) const;

  // Calculate an integral cut. This inequality has to have a greater
  // angle than the other. The generated cut will have an angle
  // inbetween the two given. The next cut can be calculated by
  // calling this function again with this inequality and the last
  // cut. The function returns NULL if there are no more cuts.
  Inequality* calculateCut(const Inequality& other) const;

  // Implement the angular relational operators.
  inline bool operator<=(const Inequality& other) const {
    return (compareWith(other)<=0);
  };
  
  inline bool operator<(const Inequality& other) const {
    return (compareWith(other)<0);
  };
  
  inline bool operator>=(const Inequality& other) const {
    return (compareWith(other)>=0);
  };
  
  inline bool operator>(const Inequality& other) const {
    return (compareWith(other)>0);
  };
  
  inline bool operator==(const Inequality& other) const {
    return (compareWith(other)==0);
  };

  std::ostream& output(std::ostream& stream,
		       DomVar varNameX=0,
		       DomVar varNameY=0) const {
    if (a>=0) stream << " ";
    stream << a;
    if (varNameX || varNameY) stream << " x_" << varNameX;
    else stream << "x";
    stream << " ";
    if (b>=0) stream << "+";
    if (varNameX || varNameY) stream << b << " x_" << varNameY;
    else stream << b << " y";
    stream << " <= " << c;
    return stream;
  };
  
  friend std::ostream& operator<<(std::ostream& stream,
				  const Inequality& e) {
    
    return e.output(stream, 0, 0);
  };

  // Determine the direction of this inequality. See Direction in
  // common.hh.
  inline Direction calcDirection() const {
    int sa = sgn(a);
    int sb = sgn(b);
    switch (sa) {
    case 1: return (sb>=0 ? east : south); break;
    case -1: return (sb>0 ? north : west); break;
    default: {
      assert(sb!=0);
      return (sb>0 ? north : south);
    }; break;
    }
  };

  declMemDbg;

};

// Calculate the highest (isUpper=true) integral value of x in the
// intersection space of the two given inequalities. The inequalities
// must be increasing in angle.
template<bool>
mpq_class extremeX(bool isUpper,
		   const Inequality& e1, const Inequality& e2);

// As above for upper y.
template<bool>
mpq_class extremeY(bool isUpper,
		   const Inequality& e1, const Inequality& e2);

// Function to tighten bounds to the closest integral point according
// to the two inequalities.
template<bool isZ>
void refineBoundsToZ(Interval<isZ>& xBound,
		     Interval<isZ>& yBound,
		     const Inequality* currentF,
		     const Inequality* nextF);


}

#endif // __PLANAR_H
