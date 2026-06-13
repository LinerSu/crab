#include "../../common.hpp"
#include "../../program_options.hpp"

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
  coeffs.insert(coeffs.end(), {10, 255});
#endif

  variable_factory_t vfac;
  z_var c(vfac["c"], crab::INT_TYPE, 32);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var n(vfac["n"], crab::INT_TYPE, 32);
  z_var o(vfac["o"], crab::INT_TYPE, 32);
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var j(vfac["j"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);

  { // tvpi paper example for C string
    // We convert into a simple program with no string but need string property
    /*
        char s[32] = "the string";
        int i = 0;
        while (true) {
            c = s[i]; <------ output
            if (c==0) break;
            i = i + 1;
        };
     */
    // dom1: i \in [0, 9], c \in [1, 255]
    test_domain_t dom1;
    dom1 += (i >= z_number(0));
    dom1 += (i <= z_number(9));     // i \in [0, 10)
    dom1 += (c <= z_number(255));
    dom1 += (c >= z_number(1));      // c \in [1, 255]
    // dom2: i = 10, c = 0
    test_domain_t dom2;
    dom2 += (i == z_number(10));
    dom2 += (c == z_number(0));

    // dom3: i > 10, c \in [0, 255]
    test_domain_t dom3;
    dom3 += (i > z_number(10));
    dom3 += (c <= z_number(255));
    dom3 += (c >= z_number(0));

    // input: i \in [0, 10]
    test_domain_t input;
    input += (i <= z_number(10));
    input += (i >= z_number(0));

    test_domain_t output = input & dom1;
    output |= input & dom2;
    output |= input & dom3;

    crab::outs() << "output: " << output << "\n";

    // TVPI domain captures:
    // i \in [0, 10], c \in [0, 255]
    // 255i + c <= 2550
    // -i - 10c <= -10

    // for us, we can only compute the range if using zones
    // besides, coefficients cannot be extended unless we know extreme points
    // for new convex hull.
    // dom1: i \in [0, 9],       c \in [1, 255]

    // dom2: i = 10, c = 0
    // dom1 join dom2:
    //  i \in [0, 10], c \in [0, 255]
  }
  return 0;
}