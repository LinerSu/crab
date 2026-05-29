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
// lincomponent.hh: a coefficient-variable pair

#ifndef LINCOMPONENT_H_
#define LINCOMPONENT_H_

#include "common.hh"
#include "interval.hh"
#include <gmpxx.h>
#include <iostream>
#include <vector>

// Define a compound datatype to hold a variable id and a coefficient.
struct LinComponentStruct {
public:
  LinComponentStruct() {};
  LinComponentStruct(mpz_class c, Variable v) : 
    coefficient(c), variable(v) {};
  mpz_class coefficient;
  Variable variable;
  inline void swap(struct LinComponentStruct other) {
    mpz_swap(coefficient.get_mpz_t(), other.coefficient.get_mpz_t());
    Variable tmp=variable;
    variable=other.variable;
    other.variable=tmp;
  };
  inline int operator<(const struct LinComponentStruct other) const {
    return variable<other.variable ||
      (variable==other.variable && coefficient<other.coefficient);
  };
  inline bool operator==(const struct LinComponentStruct other) const {
    return variable==other.variable && coefficient==other.coefficient;
  };
  friend std::ostream& operator<<(std::ostream& s,
				  const std::vector<struct LinComponentStruct>&
				  terms);
};

typedef struct LinComponentStruct LinComponent;

inline std::ostream& operator<<(std::ostream& s,
				const std::vector<struct LinComponentStruct>&
				terms) {
  for(std::vector<LinComponent>::const_iterator termIter=terms.begin();
      termIter!=terms.end(); termIter++) {
    if (termIter!=terms.begin()) {
      s << " ";
      if (termIter->coefficient>=0) s << "+";
    };
    s << termIter->coefficient << " x_" << termIter->variable;
  };
  return s;
};

// Define a compound datatype to hold a TVPI variable, a coefficient and a
// multiplicity
struct TVPIComponentStruct {
public:
  TVPIComponentStruct() : mult(Tvpi::invalidMult) {};
  TVPIComponentStruct(mpz_class c, DomVar v, Mult m=Tvpi::invalidMult) :
    coefficient(c), variable(v), mult(m) {};
  TVPIComponentStruct(const LinComponent& comp) :
    coefficient(comp.coefficient), variable(comp.variable),
    mult(Tvpi::invalidMult) {};
  mpz_class coefficient;
  DomVar variable;
  Mult mult;
  inline int operator<(const struct TVPIComponentStruct other) const {
    return variable<other.variable;
  }
};

typedef struct TVPIComponentStruct TVPIComponent;

#endif // LINCOMPONENT_H_

