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
// debug memory allocation

#ifndef __MEMORY_H
#define __MEMORY_H

#undef DEBUG_MEMORY

#ifdef DEBUG_MEMORY

#include <fstream>
#include <stdlib.h>
#include <string.h>

namespace Tvpi {
  extern std::ofstream o;
  extern size_t allocIdx;
}

#define declMemDbg \
  void* operator new(size_t size); \
  void operator delete(void* loc)

#define defMemDbg(templ,class,lowerC,upperC)	\
  templ						\
void* class::operator new(size_t size) {	\
  cerr << #upperC "("; \
  size_t* res= (size_t*) malloc(size+2*sizeof(size));			\
  size_t idx = allocIdx++;						\
  *(res++) = size;							\
  *(res++) = idx;							\
  cerr << (void*) res << "," << idx << "," << size << ")" << std::endl;	\
  return res;								\
};									\
									\
templ									\
void class::operator delete(void* loc_) {				\
  size_t* loc = (size_t*) loc_;						\
  size_t idx =*(--loc);							\
  size_t size =*(--loc);						\
  cerr << #lowerC "(" << (void*) (loc+2) <<				\
    "," << idx << ")" << std::endl;					\
  memset(loc, idx, size+2*sizeof(size));				\
  free(loc);								\
}

#else // !DEBUG_MEMORY

#define declMemDbg

#define defMemDbg(a,b,c,d)

#endif // DEBUG_MEMORY

#endif // __MEMORY_H

