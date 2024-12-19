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

    z_fixed_tvpi_domain_t dom;
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
    crab::outs() << dom << "\n";

    // assert(y <= x)? ==> 2 * i <= 3 * j?
    bool check = dom.entails(y <= x);
    crab::outs() << "assert(y <= x): " << (check ? "true" : "false") << "\n";
    idx += 1;
  }

  { // test 2
    crab::outs() << "\n\n---- case "<<idx<<"----\n\n";
    z_number SZ = z_number(3);
    z_number SLICE = z_number(2);

    z_fixed_tvpi_domain_t dom;
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

    z_fixed_tvpi_domain_t dom;
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

  {
    crab::outs() << "\n\n---- case "<<idx<<"----\n\n";
    z_fixed_tvpi_domain_t dom1, dom2;
    // i == 4
    dom1 += (i == z_number(4));
    dom2 += (i == z_number(4));

    // dom1: x == 2 * i
    dom1 += (x == 2 * i);
    // dom2: x == 3 * i
    dom2 += (x == 3 * i);

    crab::outs() << "dom1: "<< dom1  << "\n";
    crab::outs() << "dom2: "<< dom1  << "\n";

    auto dom3 = (dom1 | dom2);

    crab::outs() << "dom 1 join dom2: " << dom3 << "\n";

    bool check = dom3.entails(x >= 8);

    crab::outs() << "assert(x >= 2 * 4): " << (check ? "true" : "false") << "\n";

    check = dom3.entails(x <= 12);

    crab::outs() << "assert(x <= 3 * 4): " << (check ? "true" : "false") << "\n";
    idx += 1;
  }

  {
    crab::outs() << "\n\n---- case "<<idx<<"----\n\n";
    z_fixed_tvpi_domain_t dom, dom1, dom2;
    dom += (n >= z_number(1));
    dom += (i == z_number(0));
    dom += (i <= n);
    dom += (x == z_number(0));
    crab::outs() << "init: \n\t" << dom << "\n";


    int count = 0;
    while(count < 3) {
      dom += (i <= n);
      crab::outs() << "entry node: \n\t" << dom << "\n";
      dom1 = dom;
      dom2 = dom;

      dom1.assign(i, i + z_number(1));
      // dom1.apply(OP_ADDITION, i, i, z_number(1));
      dom2.apply(OP_ADDITION, i, i, z_number(1));

      dom1.apply(OP_ADDITION, x, x, z_number(2));
      crab::outs() << "x := x + 2\n\t" << dom1 << "\n";

      dom1.apply(OP_ADDITION, x, x, z_number(3));
      crab::outs() << "x := x + 3\n\t" << dom2 << "\n";
      dom1 |= dom2;
      crab::outs() << "After Join: \n\t" << dom1 << "\n";
      if (count == 2) {
        dom = dom || dom1;
        crab::outs() << "Widen with entry: \n\t" << dom << "\n";
      } else {
        dom |= dom1;
        crab::outs() << "Join with entry: \n\t" << dom << "\n";
      }
      count ++;
    }

    dom += (i > n);
    crab::outs() << "i > n\n\t" << dom << "\n";

    bool check = dom.entails(x >= 2 * n);
    crab::outs() << "assert(x >= 2 * n): \n\t" << (check ? "true" : "false") << "\n";
    check = dom.entails(x <= 3 * n);
    crab::outs() << "assert(x <= 3 * n): \n\t" << (check ? "true" : "false") << "\n";
    idx += 1;
  }
  return 0;
}