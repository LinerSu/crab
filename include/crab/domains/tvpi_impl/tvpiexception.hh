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
// tvpiexception.h
// Definition of exceptions that the TVPI library may throw.

#ifndef __TVPIEXCEPTION_H
#define __TVPIEXCEPTION_H

#include<iostream>

namespace Tvpi {
  class Exception;
  class IllegalArgument;
}; // namespace Tvpi

class Tvpi::Exception {};

//  class Unsatisfiable : TVPIException {
//   public:
//    Unsatisfiable() {};

//    friend ostream& operator<<(ostream& stream, const Unsatisfiable& e) {
//      return stream << "The polyhedron became unsatisfiable.";
//    };
//  };

class Tvpi::IllegalArgument : Tvpi::Exception {
  char* location;
  char* reason;
public:
  IllegalArgument(char* loc, char* rea) : location(loc), reason(rea) {};
  
  friend std::ostream& operator<<(std::ostream& stream,
				  const IllegalArgument& e) {
    return stream << e.location << ": Illegal Argument: " << e.reason;
  };
};
  


#endif // __TVPIEXCEPTION_H
