#include "../../program_options.hpp"
#include "../../common.hpp"

#include <crab/domains/tvpi_dbm.hpp>

using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace ikos;

using test_domain_t = tvpi_dbm_domain<z_sdbm_domain_t>;

unsigned idx = 1;

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

#if TVPI_DBM_FIXED_COEFFICIENTS
  auto &coeffs = crab_domain_params_man::get().coefficients();
  coeffs.insert(coeffs.end(), {2, 3, 4});
#endif

  variable_factory_t vfac;
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var n(vfac["n"], crab::INT_TYPE, 32);
  z_var o(vfac["o"], crab::INT_TYPE, 32);
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var j(vfac["j"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);

  { // test 1
    crab::outs() << "\n\n---- case "<<idx<<"----\n\n";
    z_number SZ = z_number(3);
    z_number SLICE = z_number(2);

    test_domain_t dom;
    // 0 <= i <= 6
    dom += (i >= z_number(0));
    dom += (i <= z_number(6));
    // 4 <= j <= 10
    dom += (j >= z_number(4));
    dom += (j <= z_number(10));

    // x == 3 * j
    dom += (x == SZ * j);
    // y == 2 * i
    dom += (y == SLICE * i);
    crab::outs() << "dom=" << dom << "\n";

    // assert(y <= x)? ==> 2 * i <= 3 * j?
    bool check = dom.entails(y <= x);
    crab::outs() << "assert(y <= x): " << (check ? "true" : "false") << "\n";
    idx += 1;
  }

  { // test 2
    crab::outs() << "\n\n---- case "<<idx<<"----\n\n";
    z_number SZ = z_number(3);
    z_number SLICE = z_number(2);

    test_domain_t dom;
    // 0 <= i <= 6
    dom += (i >= z_number(0));
    dom += (i <= z_number(6));
    // 0 <= j <= 10
    dom += (j >= z_number(0));
    dom += (j <= z_number(10));

    // i < j => 1 <= j <= 10
    dom += (j > i);

    // x == 3 * j
    dom += (x == SZ * j);
    // y == 2 * i
    dom += (y == SLICE * i);
    crab::outs() << dom << "\n";

    // assert(y <= x)? ==> 2 * i <= 3 * j
    // to prove this it requires coefficients {2, 3} for both i and j.
    bool check = dom.entails(y <= x);
    crab::outs() << "assert(y <= x): " << (check ? "true" : "false") << "\n";
    idx += 1;
  }

  { // test 3
    crab::outs() << "\n\n---- case "<<idx<<"----\n\n";
    z_number SZ = z_number(3);
    z_number OFFSET = z_number(2);
    z_number SLICE = z_number(2);

    test_domain_t dom;
    // 0 <= i <= 4
    dom += (i >= z_number(0));
    dom += (i <= z_number(4));
    // 1 <= j <= 4
    dom += (j >= z_number(1));
    dom += (j <= z_number(4));

    // 0 <= k < i
    dom += (k >= z_number(0));
    dom += (k <= i - 1);

    // x == 3 * i
    dom += (x == SZ * i);

    // o == k + 2
    dom += (o == OFFSET + k);

    // y == 2 * j
    dom += (y == SLICE * j);
    crab::outs() << dom << "\n";

    // assert(o <= x)? k + 2 <= 3 * i?
    // Invariant: 
    bool check = dom.entails(o <= x);
    crab::outs() << "assert(o <= x): " << (check ? "true" : "false") << "\n";

    dom += (z == o + y); // k + 2 + 2 * j < i + 2 + 2 * j

    crab::outs() << dom << "\n";

    // assert(z <= x)? i + 2 + 2 * j <= 4 * i?
    crab::outs() << "assert(z <= x): " << (check ? "true" : "false") << "\n";
    idx += 1;
  }

  { // test 4
    crab::outs() << "\n\n---- case "<<idx<<"----\n\n";
    test_domain_t dom1, dom2;
    // i == 4
    dom1 += (i == z_number(4));
    dom2 += (i == z_number(4));

    // dom1: x == 2 * i
    dom1 += (x == 2 * i);
    // dom2: x == 3 * i
    dom2 += (x == 3 * i);

    crab::outs() << "dom1: "<< dom1  << "\n";
    crab::outs() << "dom2: " << dom2 << "\n";

    auto dom3 = (dom1 | dom2);

    crab::outs() << "dom 1 join dom2: " << dom3 << "\n";

    bool check = dom3.entails(x >= 8);

    crab::outs() << "assert(x >= 2 * 4): " << (check ? "true" : "false") << "\n";

    check = dom3.entails(x <= 12);

    crab::outs() << "assert(x <= 3 * 4): " << (check ? "true" : "false") << "\n";
    idx += 1;
  }

  { // test 5
    crab::outs() << "\n\n---- case "<<idx<<"----\n\n";
    // Check if x = x + 2 join x = x + 3 works
    test_domain_t dom1, dom2, dom3;
    // dom1 : x = 0, i = 0
    dom1 += (x == z_number(0));
    dom1 += (i == z_number(0));

    // dom2:
    dom2 = dom1;
    //  crab_intrinsic(loop_counter,i:int32);
    //  i = i+1;
    dom2.intrinsic("loop_counter", {i}, {});
    dom2.apply(OP_ADDITION, i, i, z_number(1));

    crab::outs() << "Dom2 adds a loop counter " << i << "\n";

    // dom3:
    dom3 = dom2;
    //  crab_intrinsic(loop_counter,i:int32);
    //  i = i+1;
    //  dom2 for x = x+3;
    //  dom3 for x = x+2;

    //  dom2 : x = 3, i = 1, coeff_map = {i: 3}
    //  dom3 : x = 2, i = 1, coeff_map = {i: 2}
    dom2.apply(OP_ADDITION, x, x, z_number(3));
    dom3.apply(OP_ADDITION, x, x, z_number(2));

    // Perform join
    // LIMIT: join will lost coefficients for i
    test_domain_t dom4 = dom2 | dom3;
    crab::outs() << "Dom4 = Dom2 | Dom3 = " << dom4 << "\n";

    // Perform widening
    test_domain_t dom5 = dom1 || (dom1 | dom4);
    crab::outs() << "Dom5 = Dom1 || (Dom1 | Dom4) = " << dom5 << "\n";

    // Check <= order
    bool r1 = dom1 <= dom5;
    crab::outs() << "Dom1 <= Dom5 = " << (r1 ? "true" : "false") << "\n";
    idx++;
  }

  return 0;
}