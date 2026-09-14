// Independent review of the taint (tag) analysis in region_domain.
//
// These cases are written from paper/src/4-taint.tex (+ paper/assets/
// simp_trans.tex, dfa_trans.tex), not from the implementation:
//
//   * the abstract propagation is
//        tau#' = tau# \ abs-kills(Tgt)  cup  {l in Tgt | Src cap tau# != {}}
//     with abs-kills(T) = T cap Vars -- a *cell* (region) is never cleared,
//     because it summarises several concrete locations;
//   * a value defined from tainted operands is maybe-tainted.
//
// Each BOOST_CHECK below states the answer the paper's semantics requires.
// A failing case is therefore a finding, not a broken test.
#define BOOST_TEST_MODULE review_taint
#include <boost/test/unit_test.hpp>

#include "../crab_lang.hpp"
#include "../crab_dom.hpp"

#include <crab/analysis/fwd_analyzer.hpp>
#include <crab/analysis/graphs/sccg_bgl.hpp>
#include <crab/analysis/inter/inter_params.hpp>
#include <crab/analysis/inter/top_down_inter_analyzer.hpp>
#include <crab/cg/cg_bgl.hpp>
#include <crab/checkers/assertion.hpp>
#include <crab/checkers/base_property.hpp>
#include <crab/checkers/checker.hpp>
#include <crab/domains/abstract_domain_params.hpp>
#include <crab/domains/region/tags.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/support/os.hpp>

#include <string>
#include <vector>

using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::cg_impl;
using namespace crab::domain_impl;
using namespace crab::domains;

namespace {

struct counts {
  unsigned safe;
  unsigned warn;
  unsigned err;
};

z_var_or_cst_t int32_cst(int n) {
  return z_var_or_cst_t(z_number(n), crab::variable_type(crab::INT_TYPE, 32));
}

void set_params(bool havoc_clean) {
  region_domain_params p(true /*allocation_sites*/, true /*deallocation*/,
                         true /*tag_analysis*/, false /*is_dereferenceable*/,
                         true /*skip_unknown_regions*/,
                         havoc_clean /*tag_havoc_clean*/);
  crab_domain_params_man::get().update_params(p);
}

struct params_fixture {
  params_fixture() { set_params(false); }
};
BOOST_GLOBAL_FIXTURE(params_fixture);

counts intra(z_cfg_t &cfg) {
  using analyzer_t =
      crab::analyzer::intra_fwd_analyzer<z_cfg_ref_t, z_rgn_bool_int_t>;
  using checker_t = crab::checker::intra_checker<analyzer_t>;
  using assert_checker_t = crab::checker::assert_property_checker<analyzer_t>;
  z_rgn_bool_int_t init;
  crab::fixpoint_parameters fixpo_params;
  analyzer_t a(cfg, init.make_top(), nullptr, fixpo_params);
  typename analyzer_t::assumption_map_t assumptions;
  a.run(cfg.entry(), init, assumptions);
  typename checker_t::prop_checker_ptr prop(new assert_checker_t(0));
  checker_t checker(a, {prop});
  checker.run();
  crab::checker::checks_db db = checker.get_all_checks();
  return {db.get_total_safe(), db.get_total_warning(), db.get_total_error()};
}

counts inter(std::vector<z_cfg_ref_t> cfgs) {
  using analyzer_t =
      crab::analyzer::top_down_inter_analyzer<z_cg_t, z_rgn_bool_int_t>;
  z_cg_t cg(cfgs);
  crab::analyzer::inter_analyzer_parameters<z_cg_t> params;
  z_rgn_bool_int_t init;
  analyzer_t a(cg, init, params);
  a.run(init);
  crab::checker::checks_db db = a.get_all_checks();
  return {db.get_total_safe(), db.get_total_warning(), db.get_total_error()};
}

#define EXPECT_COUNTS(r, s, w)                                                 \
  do {                                                                         \
    BOOST_CHECK_EQUAL((r).safe, (unsigned)(s));                                \
    BOOST_CHECK_EQUAL((r).warn, (unsigned)(w));                                \
    BOOST_CHECK_EQUAL((r).err, 0u);                                            \
  } while (0)

} // namespace

//===--------------------------------------------------------------------===//
// (R1) move_tag on a region that summarises several concrete objects.
//
// R has two references (refcount = "one or more"), so R's cell stands for at
// least two concrete locations; one of them is tainted.  A memcpy-style
// propagation from a clean source S into R must not clear R:
// abs-kills(T) = T cap Vars forbids killing a cell.
//
// region_domain.hpp:2847 does m_tag_env.set(rgn2, find_tag_or_not(rgn1)),
// an unconditional strong overwrite with no refcount guard (contrast
// remove_tag, region_domain.hpp:2761-2769, which does guard).
//===--------------------------------------------------------------------===//
BOOST_AUTO_TEST_CASE(move_tag_must_not_clear_a_summary_region) {
  variable_factory_t vfac;
  crab::tag_manager as;
  z_var R(vfac["R"], crab::REG_INT_TYPE, 32), S(vfac["S"], crab::REG_INT_TYPE, 32);
  z_var p(vfac["p"], crab::REF_TYPE), p2(vfac["p2"], crab::REF_TYPE),
      s(vfac["s"], crab::REF_TYPE);
  z_var b(vfac["b"], crab::BOOL_TYPE);

  z_cfg_t cfg("entry", "ret");
  auto &entry = cfg.insert("entry");
  auto &ret = cfg.insert("ret");
  entry >> ret;
  entry.region_init(R);
  entry.make_ref(p, R, int32_cst(4), as.mk_tag());
  entry.make_ref(p2, R, int32_cst(4), as.mk_tag()); // R now summarises 2 objects
  entry.intrinsic("add_tag", {}, {R, p, int32_cst(1)});
  entry.region_init(S);
  entry.make_ref(s, S, int32_cst(4), as.mk_tag());  // S is clean
  entry.intrinsic("move_tag", {}, {S, s, R, p});    // clean S --> tainted R
  ret.intrinsic("check_does_not_have_tag", {b}, {R, p, int32_cst(1)});
  ret.bool_assert(b);
  // Paper: the tag on R's cell survives.
  EXPECT_COUNTS(intra(cfg), 0, 1);
}

//===--------------------------------------------------------------------===//
// (R2) A tainted *reference* passed as an actual parameter.
//
// inter_abstract_operations_impl.hpp:44-48 unifies a reference formal with
// its actual by `inv -= lhs; ref_assume(lhs == rhs)`.  Neither statement
// propagates tags: operator-= (region_domain.hpp:992-1007) makes lhs
// unknown, or *clean* when region.tag_havoc_clean is on -- which is the
// configuration the benchmarks use (yaml/clam_base.yaml:43).  ref_assume
// does not touch m_tag_env at all.
//
// main builds a tainted reference tp (int_to_ref of a tainted integer: the
// tag rule at region_domain.hpp:1869) and passes it to foo, which converts
// it back to an integer and returns it.
//===--------------------------------------------------------------------===//
static counts reference_argument_case(bool havoc_clean) {
  set_params(havoc_clean);
  variable_factory_t vfac;
  crab::tag_manager as;
  z_var f(vfac["f"], crab::REF_TYPE), RF(vfac["RF"], crab::REG_INT_TYPE, 32),
      o(vfac["o"], crab::INT_TYPE, 32);
  z_var RI(vfac["RI"], crab::REG_INT_TYPE, 32), pi(vfac["pi"], crab::REF_TYPE),
      iv(vfac["iv"], crab::INT_TYPE, 32);
  z_var RT(vfac["RT"], crab::REG_INT_TYPE, 32), tp(vfac["tp"], crab::REF_TYPE);
  z_var RO(vfac["RO"], crab::REG_INT_TYPE, 32), po(vfac["po"], crab::REF_TYPE);
  z_var res(vfac["res"], crab::INT_TYPE, 32), b(vfac["b"], crab::BOOL_TYPE);

  // int foo(ref f, region RF) { return ref_to_int(RF, f); }
  function_decl<z_number, varname_t> dfoo("foo", {f, RF}, {o});
  z_cfg_t foo("entry", "exit", dfoo);
  auto &fe = foo.insert("entry");
  auto &fx = foo.insert("exit");
  fe >> fx;
  fe.ref_to_int(RF, f, o);

  function_decl<z_number, varname_t> dmain("main", {}, {});
  z_cfg_t m("entry", "exit", dmain);
  auto &me = m.insert("entry");
  auto &mx = m.insert("exit");
  me >> mx;
  // iv is tainted
  me.region_init(RI);
  me.make_ref(pi, RI, int32_cst(4), as.mk_tag());
  me.intrinsic("add_tag", {}, {RI, pi, int32_cst(1)});
  me.load_from_ref(iv, pi, RI);
  // tp := (ref) iv  -- a tainted reference
  me.region_init(RT);
  me.int_to_ref(iv, RT, tp);
  me.region_init(RO);
  me.make_ref(po, RO, int32_cst(4), as.mk_tag());
  mx.callsite("foo", {res}, {tp, RT});
  mx.store_to_ref(po, RO, res);
  mx.intrinsic("check_does_not_have_tag", {b}, {RO, po, int32_cst(1)});
  mx.bool_assert(b);
  counts c = inter({foo, m});
  set_params(false);
  return c;
}

BOOST_AUTO_TEST_CASE(reference_actual_keeps_its_tags_at_a_call) {
  // Paper: the value of tp is tainted, so is foo's result.
  EXPECT_COUNTS(reference_argument_case(false), 0, 1);
  EXPECT_COUNTS(reference_argument_case(true), 0, 1);
}

//===--------------------------------------------------------------------===//
// (R3) region_copy overwrites the destination's tags -- documenting the
// contract, not a finding.
//
// Round 1 flagged region_copy and region_cast (region_domain.hpp:1099,
// :1175) as a kill on a cell.  Round 2 made region_cast into a *typed*
// region a union (region_cast_into_typed_region_is_a_union in
// region_taint.cc) and left region_copy a definition: crab's contract is
// "lhs_rgn := rhs_rgn", clam uses it only to rename function input
// parameters into fresh names (CfgBuilder.cc, initializeRegions puts the
// lhs in mustNotInitVars), and its numeric and allocation-site components
// overwrite as well.  I could not find a clam program that reaches
// region_copy with a non-fresh lhs.  This case therefore records the
// precondition: a non-fresh lhs loses its tags, exactly as it loses its
// numeric content.
//===--------------------------------------------------------------------===//
BOOST_AUTO_TEST_CASE(region_copy_is_a_definition_by_contract) {
  variable_factory_t vfac;
  crab::tag_manager as;
  z_var R(vfac["R"], crab::REG_INT_TYPE, 32), S(vfac["S"], crab::REG_INT_TYPE, 32);
  z_var p(vfac["p"], crab::REF_TYPE), s(vfac["s"], crab::REF_TYPE);
  z_var b(vfac["b"], crab::BOOL_TYPE);
  z_cfg_t cfg("entry", "ret");
  auto &entry = cfg.insert("entry");
  auto &ret = cfg.insert("ret");
  entry >> ret;
  entry.region_init(R);
  entry.make_ref(p, R, int32_cst(4), as.mk_tag());
  entry.intrinsic("add_tag", {}, {R, p, int32_cst(1)}); // R tainted
  entry.region_init(S);
  entry.make_ref(s, S, int32_cst(4), as.mk_tag());      // S clean
  entry.region_copy(R, S);                              // R := S  (kills R)
  ret.intrinsic("check_does_not_have_tag", {b}, {R, p, int32_cst(1)});
  ret.bool_assert(b);
  // Contract: region_copy defines lhs_rgn, so R's own tags are gone.
  EXPECT_COUNTS(intra(cfg), 1, 0);
}

//===--------------------------------------------------------------------===//
// (R4) region.tag_havoc_clean turns "no rule" into "clean".
//
// This is the policy of plan.md Decision 1.  It is recorded here as
// behaviour, not as a paper rule: 4-taint.tex has no havoc statement and no
// call/external-call rule, so nothing in the paper licenses it.  The test
// passes under the current code; it exists so that a later change of the
// default is noticed.
//===--------------------------------------------------------------------===//
BOOST_AUTO_TEST_CASE(havoc_clean_is_a_policy) {
  for (int clean = 0; clean < 2; ++clean) {
    set_params(clean);
    variable_factory_t vfac;
    crab::tag_manager as;
    z_var R(vfac["R"], crab::REG_INT_TYPE, 32), p(vfac["p"], crab::REF_TYPE),
        x(vfac["x"], crab::INT_TYPE, 32), y(vfac["y"], crab::INT_TYPE, 32),
        b(vfac["b"], crab::BOOL_TYPE);
    z_cfg_t cfg("entry", "ret");
    auto &entry = cfg.insert("entry");
    auto &ret = cfg.insert("ret");
    entry >> ret;
    entry.region_init(R);
    entry.make_ref(p, R, int32_cst(4), as.mk_tag());
    entry.intrinsic("add_tag", {}, {R, p, int32_cst(1)});
    entry.load_from_ref(x, p, R);   // x tainted
    entry.havoc(y);                 // "y := ext(x)" with no rule for ext
    entry.store_to_ref(p, R, y);    // strong store: kills R's tag
    ret.intrinsic("check_does_not_have_tag", {b}, {R, p, int32_cst(1)});
    ret.bool_assert(b);
    if (clean) {
      EXPECT_COUNTS(intra(cfg), 1, 0); // taint disappears (closed-world policy)
    } else {
      EXPECT_COUNTS(intra(cfg), 0, 1);
    }
  }
  set_params(false);
}

//===--------------------------------------------------------------------===//
// (R5) The map lattice, checked against gamma independently of
// tests/unit/tag_env_soundness.cc: absent = top, an explicit empty set =
// clean, order/join/meet pointwise.  The cases below are the ones the
// patricia_trees.hpp::compare fix (:1172-1178) is about: two single-leaf
// trees with different keys, and a leaf against a node.
//===--------------------------------------------------------------------===//
BOOST_AUTO_TEST_CASE(separate_discrete_domain_order) {
  using tag_t = crab::domains::region_domain_impl::tag<ikos::z_number>;
  using env_t = separate_discrete_domain<z_var, tag_t>;
  using set_t = typename env_t::mapped_type;

  variable_factory_t vfac;
  z_var x(vfac["x"], crab::INT_TYPE, 32), y(vfac["y"], crab::INT_TYPE, 32),
      z(vfac["z"], crab::INT_TYPE, 32);
  tag_t t1(ikos::z_number(1)), t2(ikos::z_number(2));

  auto empty = set_t::bottom();
  auto s1 = set_t(t1);
  auto s2 = set_t(t2);

  env_t top = env_t::top();
  BOOST_CHECK(top.is_top());

  // A key constrained only on the right must make leq fail: gamma(A) has
  // states with y tainted, gamma(B) has none.
  env_t A, B;
  A.set(x, empty);
  B.set(z, empty);
  BOOST_CHECK(!(A <= B));  // the compare() fix
  BOOST_CHECK(!(B <= A));

  env_t C;
  C.set(x, empty);
  C.set(y, empty);
  env_t D;
  D.set(x, empty);
  BOOST_CHECK(C <= D);     // more keys = more constraints = lower
  BOOST_CHECK(!(D <= C));

  // leaf vs node, both directions.
  env_t E;
  E.set(x, s1);
  E.set(y, s2);
  env_t F;
  F.set(x, s1);
  BOOST_CHECK(E <= F);
  BOOST_CHECK(!(F <= E));
  env_t G;
  G.set(z, s1);
  BOOST_CHECK(!(E <= G));
  BOOST_CHECK(!(G <= E));

  // Join drops one-sided keys (they become unknown); meet keeps them.
  env_t J = A | B;
  BOOST_CHECK(J.is_top());
  env_t M = A & B;
  BOOST_CHECK(M.at(x).is_bottom());
  BOOST_CHECK(M.at(z).is_bottom());

  // top/bottom of the per-key lattice.
  env_t H;
  H.set(x, s1);
  BOOST_CHECK(H.at(y).is_top());       // absent = any tag
  BOOST_CHECK(!H.at(x).is_top());
  H.set(x, set_t::top());
  BOOST_CHECK(H.is_top());             // storing top removes the key

  // forget and project move *up*.
  env_t K;
  K.set(x, s1);
  K.set(y, empty);
  env_t K1(K);
  K1 -= y;
  BOOST_CHECK(K <= K1);
  BOOST_CHECK(K1.at(y).is_top());
  env_t K2(K);
  K2.project({x});
  BOOST_CHECK(K <= K2);
  BOOST_CHECK(K2.at(y).is_top());
}

//===--------------------------------------------------------------------===//
// (R6) Completeness of the patricia_trees.hpp::compare fix.
//
// tests/unit/tag_env_soundness.cc is exhaustive over three keys, so the
// trees it builds have at most three leaves and the deep branch-vs-branch
// cases of compare() are never reached.  Here leq is checked against its
// definition (pointwise, absent = top) on random maps over eight keys, which
// produces multi-level patricia trees.
//===--------------------------------------------------------------------===//
BOOST_AUTO_TEST_CASE(leq_matches_pointwise_definition_on_deep_trees) {
  using tag_t = crab::domains::region_domain_impl::tag<ikos::z_number>;
  using env_t = separate_discrete_domain<z_var, tag_t>;
  using set_t = typename env_t::mapped_type;

  variable_factory_t vfac;
  std::vector<z_var> vars;
  for (unsigned i = 0; i < 8; ++i) {
    vars.push_back(z_var(vfac["v" + std::to_string(i)], crab::INT_TYPE, 32));
  }
  tag_t t1(ikos::z_number(1)), t2(ikos::z_number(2));
  // The four representable non-top values plus "absent" (= top).
  std::vector<set_t> vals = {set_t::bottom(), set_t(t1), set_t(t2),
                             set_t(t1) | set_t(t2)};

  auto build = [&](unsigned seed) {
    env_t e;
    unsigned s = seed;
    for (unsigned i = 0; i < vars.size(); ++i) {
      unsigned c = s % 5;
      s /= 5;
      if (c < 4) {
        e.set(vars[i], vals[c]);
      } // c == 4: leave the key absent
    }
    return e;
  };
  auto pointwise_leq = [&](const env_t &a, const env_t &b) {
    for (auto const &v : vars) {
      if (!(a.at(v) <= b.at(v))) {
        return false;
      }
    }
    return true;
  };

  unsigned x = 12345;
  unsigned mismatches = 0;
  for (unsigned n = 0; n < 20000; ++n) {
    x = x * 1103515245u + 12345u;
    env_t a = build((x >> 8) % 390625u);
    x = x * 1103515245u + 12345u;
    env_t b = build((x >> 8) % 390625u);
    if ((a <= b) != pointwise_leq(a, b)) {
      if (mismatches < 5) {
        crab::crab_string_os os_a, os_b;
        a.write(os_a);
        b.write(os_b);
        BOOST_TEST_MESSAGE("leq mismatch: " << os_a.str() << " <= " << os_b.str()
                                            << " gave " << (a <= b)
                                            << ", expected "
                                            << pointwise_leq(a, b));
      }
      ++mismatches;
    }
  }
  BOOST_CHECK_EQUAL(mismatches, 0u);

  // Join and meet must be pointwise too (absent = top), and must be the lub
  // resp. glb for the order above.
  unsigned jm_mismatches = 0;
  x = 987654321u;
  for (unsigned n = 0; n < 20000; ++n) {
    x = x * 1103515245u + 12345u;
    env_t a = build((x >> 8) % 390625u);
    x = x * 1103515245u + 12345u;
    env_t b = build((x >> 8) % 390625u);
    env_t j = a | b;
    env_t m = a & b;
    // NOTE: discrete_domain::operator== (discrete_domains.hpp:140-142) says
    // top == bottom (both carry an empty m_set), so compare is_top()
    // separately instead of relying on it.
    auto same = [](const set_t &u, const set_t &w) {
      return u.is_top() == w.is_top() && u.is_bottom() == w.is_bottom() &&
             u <= w && w <= u;
    };
    for (auto const &v : vars) {
      if (!same(j.at(v), a.at(v) | b.at(v))) {
        ++jm_mismatches;
      }
      if (!same(m.at(v), a.at(v) & b.at(v))) {
        ++jm_mismatches;
      }
    }
    if (!(a <= j) || !(b <= j) || !(m <= a) || !(m <= b)) {
      ++jm_mismatches;
    }
  }
  BOOST_CHECK_EQUAL(jm_mismatches, 0u);
}

//===--------------------------------------------------------------------===//
// (R7) The sink verdict is produced by asserting a boolean that
// check_does_not_have_tag either sets to true (region_domain.hpp:2785) or
// havocs (:2802).  So an assert on an *undefined* boolean must warn,
// otherwise a sink on an unknown region could be reported safe.
//
// This is worth pinning because crab/tests/domains/region/region-8.cc uses
// the intrinsic name "does_not_have_tag", which the region domain does not
// match (:2771); the fallback at :2916 forwards it to the base domain, and
// flat_boolean_domain::intrinsic (flat_boolean_domain.hpp:440-444) is a
// no-op -- so region-8's five asserts are on a boolean no statement defines,
// and it still reports 4 safe / 1 warning.
//===--------------------------------------------------------------------===//
BOOST_AUTO_TEST_CASE(assert_on_an_undefined_bool_must_warn) {
  variable_factory_t vfac;
  crab::tag_manager as;
  z_var R(vfac["R"], crab::REG_INT_TYPE, 32), p(vfac["p"], crab::REF_TYPE),
      b(vfac["b"], crab::BOOL_TYPE);
  z_cfg_t cfg("entry", "ret");
  auto &entry = cfg.insert("entry");
  auto &ret = cfg.insert("ret");
  entry >> ret;
  entry.region_init(R);
  entry.make_ref(p, R, int32_cst(4), as.mk_tag());
  // b is never defined.
  ret.bool_assert(b);
  EXPECT_COUNTS(intra(cfg), 0, 1);
}

// (R7b) ... and the reason region-8 still shows 4 safe is that `assert(b)`
// also *assumes* b afterwards, so the second and later asserts on the same
// undefined boolean are trivially safe.  region-8 therefore contains exactly
// one real check, and it warns because no statement defines b1.
BOOST_AUTO_TEST_CASE(repeated_assert_on_the_same_bool_is_assumed) {
  variable_factory_t vfac;
  crab::tag_manager as;
  z_var R(vfac["R"], crab::REG_INT_TYPE, 32), p(vfac["p"], crab::REF_TYPE),
      b(vfac["b"], crab::BOOL_TYPE);
  z_cfg_t cfg("entry", "ret");
  auto &entry = cfg.insert("entry");
  auto &ret = cfg.insert("ret");
  entry >> ret;
  entry.region_init(R);
  entry.make_ref(p, R, int32_cst(4), as.mk_tag());
  ret.bool_assert(b);
  ret.bool_assert(b);
  ret.bool_assert(b);
  EXPECT_COUNTS(intra(cfg), 2, 1);
}

//===--------------------------------------------------------------------===//
// (R8) Round 2: ref_assume now *reduces* the tags of both sides of a
// reference equality to their intersection (region_domain.hpp, ref_assume,
// "p == q: both references denote the same value").
//
// That is the right thing for the renaming in
// inter_transformers_impl::unify (forget(lhs) makes tags(lhs) = top, so the
// intersection is tags(rhs)).  But ref_assume is also the transfer function
// of the program statement `assume(p == q)`, and there the refinement does
// not hold: in the paper's concrete semantics taint is attached to
// *locations*, not to values (4-taint.tex: "tau in P(Locs) records tainted
// locations"), and it is written by the statement that defines the variable.
// Two variables can hold the same address with different taint, e.g.
//
//     t = *pp;        // prop(cell(pp), t):  t is tainted if *pp's cell is
//     u = make_ref(); // clean
//     assume(t == u); // both denote the same address; sigma(t) = {1} still
//
// so intersecting loses the tag of t.  Below, RP (a region of references)
// is tainted and holds u; t is loaded from it, then compared with u.
//===--------------------------------------------------------------------===//
BOOST_AUTO_TEST_CASE(ref_assume_equality_must_not_intersect_tags) {
  variable_factory_t vfac;
  crab::tag_manager as;
  z_var RP(vfac["RP"], crab::REG_REF_TYPE), pp(vfac["pp"], crab::REF_TYPE);
  z_var RT(vfac["RT"], crab::REG_INT_TYPE, 32), u(vfac["u"], crab::REF_TYPE),
      t(vfac["t"], crab::REF_TYPE);
  z_var RO(vfac["RO"], crab::REG_INT_TYPE, 32), po(vfac["po"], crab::REF_TYPE),
      o(vfac["o"], crab::INT_TYPE, 32), b(vfac["b"], crab::BOOL_TYPE);

  z_cfg_t cfg("entry", "ret");
  auto &entry = cfg.insert("entry");
  auto &ret = cfg.insert("ret");
  entry >> ret;
  entry.region_init(RT);
  entry.make_ref(u, RT, int32_cst(4), as.mk_tag());      // u: clean reference
  entry.region_init(RP);
  entry.make_ref(pp, RP, int32_cst(8), as.mk_tag());
  entry.store_to_ref(pp, RP, u);                        // *pp = u
  entry.intrinsic("add_tag", {}, {RP, pp, int32_cst(1)}); // the cell is tainted
  entry.load_from_ref(t, pp, RP);                       // t = *pp  -> tainted
  entry.assume_ref(z_ref_cst_t::mk_eq(t, u, z_number(0)));
  entry.ref_to_int(RT, t, o);                           // o = (int) t
  entry.region_init(RO);
  entry.make_ref(po, RO, int32_cst(4), as.mk_tag());
  entry.store_to_ref(po, RO, o);
  ret.intrinsic("check_does_not_have_tag", {b}, {RO, po, int32_cst(1)});
  ret.bool_assert(b);
  // Paper: t is still tainted on the branch where t == u.
  EXPECT_COUNTS(intra(cfg), 0, 1);
}

//===--------------------------------------------------------------------===//
// (R9) Round 2, completeness of the F4 fix: region_domain::caller_continuation
// saves and unions back the caller's tags only for output-only regions whose
// type is *known* ("unknown-typed outputs are fresh views").  That holds for
// the aux regions clam creates when the callsite and the callee disagree on
// the type (CfgBuilder.cc, unifyRgnType), but a callsite whose actual and
// formal are *both* unknown-typed gets no aux: the caller's own
// unknown-typed region is then an output-only parameter and is not saved.
//===--------------------------------------------------------------------===//
BOOST_AUTO_TEST_CASE(unknown_typed_new_region_accumulates_across_calls) {
  variable_factory_t vfac;
  crab::tag_manager as;
  z_var U(vfac["U"], crab::REG_UNKNOWN_TYPE);
  z_var r(vfac["r"], crab::REF_TYPE), q(vfac["q"], crab::REF_TYPE);
  z_var b1(vfac["b1"], crab::BOOL_TYPE);

  function_decl<z_number, varname_t> dfoo("foo", {}, {U});
  z_cfg_t foo("entry", "exit", dfoo);
  auto &fe = foo.insert("entry");
  auto &fx = foo.insert("exit");
  fe >> fx;
  fe.region_init(U);
  fe.make_ref(r, U, int32_cst(4), as.mk_tag());

  function_decl<z_number, varname_t> dmain("main", {}, {});
  z_cfg_t m("entry", "exit", dmain);
  auto &me = m.insert("entry");
  auto &mx = m.insert("exit");
  me >> mx;
  me.region_init(U);
  me.make_ref(q, U, int32_cst(4), as.mk_tag());
  me.intrinsic("add_tag", {}, {U, q, int32_cst(2)}); // the caller tags an object
  mx.callsite("foo", {U}, {});                       // output-only, unknown type
  mx.intrinsic("check_does_not_have_tag", {b1}, {U, q, int32_cst(2)});
  mx.bool_assert(b1);
  // The caller's object still carries tag 2 after the call.
  EXPECT_COUNTS(inter({foo, m}), 0, 1);
}
