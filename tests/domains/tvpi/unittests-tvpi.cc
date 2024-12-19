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

void perfrom_domain_operations(const test_domain_t& dom1, const test_domain_t &dom2) {
  crab::outs() << "\n\n---- case "<<idx<<"----\n\n";
  crab::outs() << "Dom1=" << dom1 << "\n";
  crab::outs() << "Dom2=" << dom2 << "\n";
  bool r1 = dom1 <= dom2;
  crab::outs() << "Dom1 <= Dom2 = " << (r1 ? "true" : "false") << "\n";
  bool r2 = dom2 <= dom1;
  crab::outs() << "Dom2 <= Dom1 = " << (r2 ? "true" : "false") << "\n";
  test_domain_t dom3 = dom1 | dom2;
  crab::outs() << "Dom3 = Dom1 | Dom2 = " << dom3 << "\n";
  bool r3 = dom1 <= dom3;
  bool r4 = dom2 <= dom3;
  crab::outs() << "Dom1 <= Dom3 = " << (r3 ? "true" : "false") << "\n";
  crab::outs() << "Dom2 <= Dom3 = " << (r4 ? "true" : "false") << "\n";
  test_domain_t dom4 = dom1 & dom2;
  crab::outs() << "Dom4 = Dom1 & Dom2 = " << dom4 << "\n";
  bool r5 = dom4 <= dom1;
  bool r6 = dom4 <= dom2;
  crab::outs() << "Dom4 <= Dom1 = " << (r5 ? "true" : "false") << "\n";
  crab::outs() << "Dom4 <= Dom2 = " << (r6 ? "true" : "false") << "\n";
  idx ++;
}

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


  {
    test_domain_t dom1, dom2;
    // dom1 : x = 1, y = 2x, z = 3x
    dom1.assign(x, z_number(1));
    dom1.apply(OP_MULTIPLICATION, y, x, z_number(2));
    dom1.apply(OP_MULTIPLICATION, z, x, z_number(3));

    // dom2 : x = 2, y = 3x, z = 4x
    dom2 += (x == z_number(2));
    dom2 += (y == x * z_number(3) );
    dom2 += (z == x + x + x + x);

    perfrom_domain_operations(dom1, dom2);
    // meet: bot
    // join: x \in 1, 2, y \in 2, 6, z \in 3, 8

    test_domain_t dom5(dom1);
    dom5 -= x;
    // dom 5: {y = 2, z = 3}
    crab::outs() << "After forgetting " << x << " in Dom 1:" << dom5 << "\n";

    test_domain_t dom6(dom1);
    dom6.rename({x}, {i});
    // dom6 : i = 1, y = 2i, z = 3i
    crab::outs() << "After renaming {x} with {i} in Dom1:" << dom6 << "\n";

    test_domain_t dom7(dom1);
    // dom 7: {y = 2, z = 3}
    dom7.project({y, z});
    crab::outs() << "After projecting on y and z in Dom1:" << dom7 << "\n";
  }

  {
    test_domain_t dom1, dom2;
    // dom1 : x >= 1, y = 2x, z = 3x
    dom1 += (x >= z_number(1));
    dom1.apply(OP_MULTIPLICATION, y, x, z_number(2));
    dom1.apply(OP_MULTIPLICATION, z, x, z_number(3));

    // dom2: x <= 20, i = 3x, j = 4x
    dom2 += (x <= z_number(20));
    dom2.assign(i, z_number(3) * x);
    dom2.assign(j, z_number(4) * x);

    perfrom_domain_operations(dom1, dom2);
  }

  {
    crab::outs() << "\n\n---- case "<<idx<<"---- \n\n";
    test_domain_t dom1, dom2;
    // dom1 : x >= 5, y = 4x + 2, z = 2x + 3y
    dom1 += (x >= z_number(5));
    dom1.assign(y, z_number(4) * x + z_number(2));
    dom1.assign(z, z_number(2) * x + z_number(3) * y);
    crab::outs() << "Dom1=" << dom1 << "\n";

    // dom2 : x \in [1, 10], y = 4x + 17, i \in [0, x), z = 4i + 3
    dom2 += (x >= z_number(1));
    dom2 += (x <= z_number(10));
    dom2.assign(y, z_number(4) * x + z_number(17));
    dom2 += (i >= z_number(0));
    dom2 += (i <= x - z_number(1));
    dom2.assign(z, z_number(4) * i + z_number(3));

    crab::outs() << "Dom2=" << dom2 << "\n";
    // To prove:        4i + 3 <= 4x + 17
    // need to prove:   4i <= 4x + 14
    // since            i - x <= -1 => 4i - 4x <= -4
    // -4 <= 14
    bool check = dom2.entails(z <= y);
    crab::outs() << "assert(z <= y): " << (check ? "true" : "false") << "\n";
    idx ++;
  }

  // {
  //   crab::outs() << "\n\n---- case "<<idx<<"---- \n\n";
  //   test_domain_t dom1, dom2;
  //   // dom2 : x \in [1, 10], y = 4x + 17, i \in [0, x), z = 4i + 3
  //   dom2 += (x >= z_number(1));
  //   dom2 += (x <= z_number(10));
  //   dom2.assign(y, z_number(4) * x + z_number(17));
  //   dom2 += (i >= z_number(0));
  //   dom2 += (i <= x - z_number(1));
  //   dom2.assign(z, z_number(4) * i + z_number(3));
  //   idx ++;
  // }

  {
    crab::outs() << "\n\n---- case "<<idx<<"---- \n\n";
    // dom1 : x - y <= 4, 2y - x <= -3, 3x - y <= 5
    test_domain_t dom1;
    dom1 += (x - y <= z_number(4));
    dom1 += (z_number(2) * y - x <= z_number(-3));
    dom1 += (z_number(3) * x - y <= z_number(5));
    crab::outs() << "original dom1=" << dom1 << "\n";
    dom1.reduce();
    crab::outs() << "dom1=" << dom1 << "\n";

    // dom2 :y - x <= 2, y <= 0, x <= 2
    test_domain_t dom2;
    dom2 += (y <= z_number(0));
    dom2 += (x <= z_number(2));
    dom2 += (y - x <= z_number(2));
    dom2.reduce();
    crab::outs() << "dom2=" << dom2 << "\n";

    test_domain_t dom3 = dom1 | dom2;
    dom3.reduce();
    // for tvpi: y<= 0, x <= 2, 2y - x <= 2
    crab::outs() << "dom1 | dom2=" << dom3 << "\n";

    test_domain_t dom4;
    dom4 += (y <= z_number(0));
    dom4 += (x <= z_number(2));
    dom4 += (z_number(2) * y - x <= z_number(2));
    dom4.reduce();
    crab::outs() << "dom4=" << dom4 << "\n";
    bool l1 = dom1 <= dom4;
    bool l2 = dom2 <= dom4;
    crab::outs() << "Dom1 <= Dom4 = " << (l1 ? "true" : "false") << "\n";
    crab::outs() << "Dom2 <= Dom4 = " << (l2 ? "true" : "false") << "\n";
    idx ++;
  }

  {
    crab::outs() << "\n\n---- case "<<idx<<"---- \n\n";
    // dom1 : 2x - y <= 5, 3y - x <= 7, -3x + y <= 4
    test_domain_t dom1;
    dom1 += (z_number(2) * x - y <= z_number(5));
    dom1 += (z_number(3) * y - x <= z_number(7));
    dom1 += (- z_number(3) * x + y <= z_number(4));
    crab::outs() << "original dom1=" << dom1 << "\n";
    dom1.reduce();
    crab::outs() << "dom1=" << dom1 << "\n";

    // dom2 : 3x − 2y <= 8, 4y - x <= 10, -2x - 3y <= -6
    test_domain_t dom2;
    dom2 += (z_number(3) * x - z_number(2) * y <= z_number(8));
    dom2 += (z_number(4) * y - x <= z_number(10));
    dom2 += (z_number(2) * x - z_number(3) * y <= z_number(-6));
    crab::outs() << "original dom2=" << dom2 << "\n";
    dom2.reduce();
    crab::outs() << "dom2=" << dom2 << "\n";

    test_domain_t dom3 = dom1 | dom2;
    crab::outs() << "dom1 | dom2=" << dom3 << "\n";

    test_domain_t dom4;
    dom4 += (y <= z_number(0));
    dom4 += (x <= z_number(2));
    dom4 += (z_number(2) * y - x <= z_number(2));
    dom4.reduce();
    crab::outs() << "dom4=" << dom4 << "\n";
    bool l1 = dom1 <= dom4;
    bool l2 = dom2 <= dom4;
    crab::outs() << "Dom1 <= Dom4 = " << (l1 ? "true" : "false") << "\n";
    crab::outs() << "Dom2 <= Dom4 = " << (l2 ? "true" : "false") << "\n";
    idx ++;
  }

  {
    crab::outs() << "\n\n---- case "<<idx<<"---- \n\n";
    z_var isz(vfac["isz"], crab::INT_TYPE, 32);
    z_var len(vfac["len"], crab::INT_TYPE, 32);
    z_var tsz(vfac["tsz"], crab::INT_TYPE, 32);
    z_var tmp(vfac["tmp"], crab::INT_TYPE, 32);
    z_var tmp2(vfac["tmp2"], crab::INT_TYPE, 32);
    z_var offset(vfac["offset"], crab::INT_TYPE, 32);

    test_domain_t dom1;
    dom1 += (isz == z_number(4));                   // isz = 4
    dom1 += (len >= z_number(1));
    dom1 += (len <= z_number(10));                  // len \in [1, 10]
    dom1.apply(OP_MULTIPLICATION, tmp, isz, len);
    dom1 += (tmp <= tsz);                           // isz * len <= tsz => 4len <= tsz
    dom1 += (i >= z_number(0));
    dom1 += (i <= len - 1);                         // i \in [0, len)   => i <= len - 1
    dom1.apply(OP_MULTIPLICATION, offset, isz, i);  // offset = isz * i => 4i = offset
    crab::outs() << "dom1=" << dom1 << "\n";
    dom1.apply(OP_ADDITION, tmp2, offset, isz);
    bool ret = dom1.entails(tmp2 <= tsz);           // offset + isz => 4i + 4 <= 4(len - 1) + 4 = 4len <= tsz
    crab::outs() << "assert(offset + isz <= tsz): " << (ret ? "true" : "false") << "\n";
    idx ++;
  }

  {
    crab::outs() << "\n\n---- case "<<idx<<"---- \n\n";
    test_domain_t dom1;
    // dom1 : x \in [1, 10], y = 2x, y <= z
    dom1 += (x >= z_number(1));
    dom1 += (x <= z_number(10));
    dom1.apply(OP_MULTIPLICATION, y, x, z_number(2));
    dom1 += (y <= z);

    // i \in [0, x), n = 2i
    dom1 += (i >= z_number(0));
    dom1 += (i <= x - z_number(1));
    dom1.apply(OP_MULTIPLICATION, n, i, z_number(2));

    // n <= z?
    bool ret = dom1.entails(n <= z);
    crab::outs() << "assert(n <= z): " << (ret ? "true" : "false") << "\n";

    idx ++;
  }
  return 0;
}