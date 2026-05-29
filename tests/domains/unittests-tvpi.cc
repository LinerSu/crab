#include "../common.hpp"
#include "../program_options.hpp"

using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace ikos;

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  variable_factory_t vfac;

  { // top / bottom basics
    z_tvpi_domain_t top;
    z_tvpi_domain_t bot = top.make_bottom();
    crab::outs() << "top=" << top << "\n";
    crab::outs() << "bot=" << bot << "\n";
    crab::outs() << "top.is_top()=" << top.is_top() << "\n";
    crab::outs() << "bot.is_bottom()=" << bot.is_bottom() << "\n";
    crab::outs() << "bot<=top? " << (bot <= top) << "\n";
    crab::outs() << "top<=bot? " << (top <= bot) << "\n";
    crab::outs() << "top<=top? " << (top <= top) << "\n";
    crab::outs() << "bot<=bot? " << (bot <= bot) << "\n";
  }

  { // assign constant and interval query
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_tvpi_domain_t d;
    d.assign(x, z_number(5));
    crab::outs() << "After x:=5: " << d << "\n";
    crab::outs() << "at(x)=" << d.at(x) << "\n";
  }

  { // assign variable copy: x = y
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_tvpi_domain_t d;
    d += (y >= z_number(0));
    d += (y <= z_number(10));
    d.assign(x, z_lin_exp_t(y));
    crab::outs() << "After y=[0,10], x:=y: " << d << "\n";
    crab::outs() << "at(x)=" << d.at(x) << "\n";
  }

  { // affine assign: x = 2*y + 3, y in [1,5]
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_tvpi_domain_t d;
    d += (y >= z_number(1));
    d += (y <= z_number(5));
    d.assign(x, z_number(2) * z_lin_exp_t(y) + z_number(3));
    crab::outs() << "After y=[1,5], x:=2y+3: " << d << "\n";
    crab::outs() << "at(x)=" << d.at(x) << "\n";
  }

  { // apply: x = y + z (both bounded)
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var z(vfac["z"], crab::INT_TYPE, 32);
    z_tvpi_domain_t d;
    d += (y >= z_number(1)); d += (y <= z_number(3));
    d += (z >= z_number(2)); d += (z <= z_number(4));
    d.apply(OP_ADDITION, x, y, z);
    crab::outs() << "After y=[1,3],z=[2,4], x:=y+z: " << d << "\n";
    crab::outs() << "at(x)=" << d.at(x) << "\n";
  }

  { // meet: two domains agree
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_tvpi_domain_t d1, d2;
    d1 += (x >= z_number(0)); d1 += (x <= z_number(10));
    d2 += (x >= z_number(5)); d2 += (x <= z_number(20));
    z_tvpi_domain_t meet = d1 & d2;
    crab::outs() << "d1=" << d1 << "\n";
    crab::outs() << "d2=" << d2 << "\n";
    crab::outs() << "d1 & d2=" << meet << "\n";
    crab::outs() << "at(x) in meet=" << meet.at(x) << "\n";
  }

  { // meet: bottom when contradictory
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_tvpi_domain_t d1, d2;
    d1 += (x >= z_number(10));
    d2 += (x <= z_number(5));
    z_tvpi_domain_t meet = d1 & d2;
    crab::outs() << "Contradictory meet is_bottom=" << meet.is_bottom() << "\n";
  }

  { // join: two intervals for x
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_tvpi_domain_t d1, d2;
    d1 += (x >= z_number(0)); d1 += (x <= z_number(5));
    d2 += (x >= z_number(8)); d2 += (x <= z_number(10));
    z_tvpi_domain_t join = d1 | d2;
    crab::outs() << "d1=" << d1 << "\n";
    crab::outs() << "d2=" << d2 << "\n";
    crab::outs() << "d1 | d2=" << join << "\n";
    crab::outs() << "at(x) in join=" << join.at(x) << "\n";
  }

  { // join: with different variable sets
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_tvpi_domain_t d1, d2;
    d1 += (x == z_number(3));
    d2 += (y == z_number(7));
    z_tvpi_domain_t join = d1 | d2;
    crab::outs() << "Join with different vars: " << join << "\n";
  }

  { // join with bottom
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_tvpi_domain_t d, bot = d.make_bottom();
    d += (x >= z_number(0)); d += (x <= z_number(5));
    z_tvpi_domain_t join = d | bot;
    crab::outs() << "d | bottom=" << join << "\n";
    crab::outs() << "bottom | d=" << (bot | d) << "\n";
  }

  { // operator<=: inclusion check
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_tvpi_domain_t narrow, wide;
    narrow += (x >= z_number(2)); narrow += (x <= z_number(8));
    wide   += (x >= z_number(0)); wide   += (x <= z_number(10));
    crab::outs() << "narrow=[2,8] <= wide=[0,10]? " << (narrow <= wide) << "\n";
    crab::outs() << "wide <= narrow? " << (wide <= narrow) << "\n";
    crab::outs() << "narrow <= narrow? " << (narrow <= narrow) << "\n";
  }

  { // widening: basic loop pattern
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_tvpi_domain_t pre, post;
    pre += (x >= z_number(0)); pre += (x <= z_number(10));
    post += (x >= z_number(0)); post += (x <= z_number(20));
    z_tvpi_domain_t widened = pre || post;
    crab::outs() << "pre=" << pre << " post=" << post << "\n";
    crab::outs() << "pre || post=" << widened << "\n";
    // widening should make x unbounded above
    crab::outs() << "at(x) after widen=" << widened.at(x) << "\n";
  }

  { // forget: remove a variable
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_tvpi_domain_t d;
    d += (x == z_lin_exp_t(y));
    d += (x >= z_number(0)); d += (x <= z_number(5));
    crab::outs() << "Before forget x: " << d << "\n";
    d -= x;
    crab::outs() << "After forget x: " << d << "\n";
  }

  { // relational constraint: x + y <= 10
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_tvpi_domain_t d;
    d += (x >= z_number(0)); d += (y >= z_number(0));
    d += z_lin_cst_t(z_lin_exp_t(x) + z_lin_exp_t(y) <= z_number(10));
    crab::outs() << "x>=0, y>=0, x+y<=10: " << d << "\n";
    crab::outs() << "at(x)=" << d.at(x) << " at(y)=" << d.at(y) << "\n";
  }

  { // equality: x == y
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_tvpi_domain_t d;
    d += (z_lin_exp_t(x) == z_lin_exp_t(y));
    d += (x >= z_number(3)); d += (x <= z_number(7));
    crab::outs() << "x==y, x=[3,7]: " << d << "\n";
    crab::outs() << "at(y)=" << d.at(y) << "\n";
  }

  { // to_linear_constraint_system
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_tvpi_domain_t d;
    d += (x >= z_number(1)); d += (x <= z_number(5));
    d += (y >= z_number(2)); d += (y <= z_number(8));
    d += z_lin_cst_t(z_lin_exp_t(x) + z_lin_exp_t(y) <= z_number(10));
    auto csts = d.to_linear_constraint_system();
    crab::outs() << "Constraints from domain: " << csts << "\n";
  }

  { // domain_name
    z_tvpi_domain_t d;
    crab::outs() << "domain_name=" << d.domain_name() << "\n";
  }

  { // Precision comparison: TVPI join vs DBM join
    //
    // s1: -i<=0 /\ i<=9 /\ -c<=10 /\ c<=-1   (i in [0,9], c in [-10,-1])
    // s2:  i=10 /\ c=0
    //
    // Expected TVPI convex-hull join:
    //   -i<=0 /\ i<=10 /\ -c<=10 /\ c<=0 /\ 10c-i<=-10 /\ 10i-c<=100
    //
    // DBM/zones join loses the two-variable facets, giving only:
    //   -i<=0 /\ i<=10 /\ -c<=10 /\ c<=0
    //
    // The spurious point (i=0, c=0) lies in the DBM box but NOT in the
    // TVPI convex hull (it violates 10c-i<=-10: 10*0-0=0 > -10).
    z_var i(vfac["i"], crab::INT_TYPE, 32);
    z_var c(vfac["c"], crab::INT_TYPE, 32);

    z_tvpi_domain_t s1;
    s1 += (i >= z_number(0));
    s1 += (i <= z_number(9));
    s1 += (c >= z_number(-10));
    s1 += (c <= z_number(-1));

    z_tvpi_domain_t s2;
    s2 += (i == z_number(10));
    s2 += (c == z_number(0));

    crab::outs() << "s1=" << s1 << "\n";
    crab::outs() << "s2=" << s2 << "\n";

    z_tvpi_domain_t tvpi_join = s1 | s2;
    crab::outs() << "TVPI s1|s2=" << tvpi_join << "\n";
    crab::outs() << "TVPI at(i)=" << tvpi_join.at(i)
                 << "  at(c)=" << tvpi_join.at(c) << "\n";

    // DBM join for comparison
    z_sdbm_domain_t s1_dbm, s2_dbm;
    s1_dbm += (i >= z_number(0)); s1_dbm += (i <= z_number(9));
    s1_dbm += (c >= z_number(-10)); s1_dbm += (c <= z_number(-1));
    s2_dbm += (i == z_number(10)); s2_dbm += (c == z_number(0));
    z_sdbm_domain_t dbm_join = s1_dbm | s2_dbm;
    crab::outs() << "DBM  s1|s2=" << dbm_join << "\n";

    // Witness point: (i=5, c=0).
    // - TVPI facet 10c-i<=-10: 10*0-5 = -5, and -5 <= -10 is FALSE
    //   => (5,0) violates the TVPI constraint => TVPI join excludes it.
    // - DBM constraints: i=5 in [0,10], c=0 in [-10,0], c-i=-5<=-1, i-c=5<=19
    //   => (5,0) satisfies all DBM constraints => DBM join includes it.
    z_tvpi_domain_t tvpi_check = tvpi_join;
    tvpi_check += (i == z_number(5));
    tvpi_check += (c == z_number(0));
    crab::outs() << "TVPI join + (i=5,c=0) is_bottom="
                 << tvpi_check.is_bottom()
                 << "  (expected 1: excluded by 10c-i<=-10)\n";

    z_sdbm_domain_t dbm_check = dbm_join;
    dbm_check += (i == z_number(5));
    dbm_check += (c == z_number(0));
    crab::outs() << "DBM  join + (i=5,c=0) is_bottom="
                 << dbm_check.is_bottom()
                 << "  (expected 0: DBM cannot rule it out)\n";
  }

  return 0;
}
