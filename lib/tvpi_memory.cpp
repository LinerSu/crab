/* Derived from the TVPI library (original file: memory.cpp).
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
#include <crab/domains/tvpi_impl/memory.hh>

#ifdef DEBUG_MEMORY

#include <iostream>

namespace Tvpi {
  std::ofstream o("allocations.txt");
  size_t allocIdx = 0;
}


#endif
