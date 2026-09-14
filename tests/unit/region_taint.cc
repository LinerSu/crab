// Region-level checks of the taint (tag) analysis: the transfer
// functions of region_domain against the intended semantics.
//
//   clean  = explicit empty tag set, written by a rule
//   T      = maybe tainted, only by tags in T
//   unknown= not in the tag environment (any tag); a sink check on an
//            unknown region warns.
//
// Each test builds a small CFG, runs the intra (or top-down inter)
// analyzer with the assertion checker and compares the number of
// safe / warning checks with what the semantics predicts.
#define BOOST_TEST_MODULE region_taint
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

// Tag analysis on; havoc'ed values unknown unless a test says otherwise.
struct params_fixture {
  params_fixture() { set_havoc_clean(false); }
  static void set_havoc_clean(bool b) {
    region_domain_params p(true /*allocation_sites*/, true /*deallocation*/,
                           true /*tag_analysis*/, false /*is_dereferenceable*/,
                           true /*skip_unknown_regions*/, b /*tag_havoc_clean*/);
    crab_domain_params_man::get().update_params(p);
  }
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

// Common variables for the intra tests.
struct vars {
  variable_factory_t vfac;
  crab::tag_manager as;
  z_var c, i, x, b, b2, p, q, R, S;
  vars()
      : c(vfac["c"], crab::INT_TYPE, 32), i(vfac["i"], crab::INT_TYPE, 32),
        x(vfac["x"], crab::INT_TYPE, 32), b(vfac["b"], crab::BOOL_TYPE),
        b2(vfac["b2"], crab::BOOL_TYPE), p(vfac["p"], crab::REF_TYPE),
        q(vfac["q"], crab::REF_TYPE), R(vfac["R"], crab::REG_INT_TYPE, 32),
        S(vfac["S"], crab::REG_INT_TYPE, 32) {}
};

} // namespace

// (i) The witness of docs/taint_dfa/issue.md: a region initialised
// and tainted on one path only, checked after the join. The other
// path says nothing about R, so after the join R is unknown and the
// check must warn.
BOOST_AUTO_TEST_CASE(one_sided_join_warns) {
  vars v;
  z_cfg_t cfg("entry", "ret");
  auto &entry = cfg.insert("entry");
  auto &bt = cfg.insert("bb_t");
  auto &bf = cfg.insert("bb_f");
  auto &ret = cfg.insert("ret");
  entry >> bt;
  entry >> bf;
  bt >> ret;
  bf >> ret;
  entry.havoc(v.c);
  bt.assume(v.c >= 1);
  bt.region_init(v.R);
  bt.make_ref(v.p, v.R, int32_cst(4), v.as.mk_tag());
  bt.intrinsic("add_tag", {}, {v.R, v.p, int32_cst(1)});
  bf.assume(v.c <= 0);
  ret.intrinsic("check_does_not_have_tag", {v.b}, {v.R, v.p, int32_cst(1)});
  ret.bool_assert(v.b);
  EXPECT_COUNTS(intra(cfg), 0, 1);
}

// (ii) Both paths define R (the SMALL path_3 shape): tainted on one,
// overwritten with a constant on the other -> maybe tainted; a
// constant on both -> clean.
BOOST_AUTO_TEST_CASE(both_sided_join) {
  for (int tainted_branch = 0; tainted_branch < 2; ++tainted_branch) {
    vars v;
    z_cfg_t cfg("entry", "ret");
    auto &entry = cfg.insert("entry");
    auto &bt = cfg.insert("bb_t");
    auto &bf = cfg.insert("bb_f");
    auto &ret = cfg.insert("ret");
    entry >> bt;
    entry >> bf;
    bt >> ret;
    bf >> ret;
    entry.region_init(v.R);
    entry.make_ref(v.p, v.R, int32_cst(4), v.as.mk_tag());
    entry.havoc(v.c);
    bt.assume(v.c >= 1);
    if (tainted_branch) {
      bt.intrinsic("add_tag", {}, {v.R, v.p, int32_cst(1)});
    } else {
      bt.store_to_ref(v.p, v.R, int32_cst(1));
    }
    bf.assume(v.c <= 0);
    bf.store_to_ref(v.p, v.R, int32_cst(0));
    ret.intrinsic("check_does_not_have_tag", {v.b}, {v.R, v.p, int32_cst(1)});
    ret.bool_assert(v.b);
    if (tainted_branch) {
      EXPECT_COUNTS(intra(cfg), 0, 1);
    } else {
      EXPECT_COUNTS(intra(cfg), 1, 0);
    }
  }
}

// (iii) Loop: a source in the body taints R for the check after the
// loop; a body that only stores constants leaves R clean.
BOOST_AUTO_TEST_CASE(loop_join) {
  for (int tainted_body = 0; tainted_body < 2; ++tainted_body) {
    vars v;
    z_cfg_t cfg("entry", "ret");
    auto &entry = cfg.insert("entry");
    auto &head = cfg.insert("head");
    auto &body = cfg.insert("body");
    auto &ret = cfg.insert("ret");
    entry >> head;
    head >> body;
    head >> ret;
    body >> head;
    entry.region_init(v.R);
    entry.make_ref(v.p, v.R, int32_cst(4), v.as.mk_tag());
    entry.assign(v.i, 0);
    body.assume(v.i <= 9);
    if (tainted_body) {
      body.intrinsic("add_tag", {}, {v.R, v.p, int32_cst(1)});
    } else {
      body.store_to_ref(v.p, v.R, int32_cst(0));
    }
    body.add(v.i, v.i, 1);
    ret.assume(v.i >= 10);
    ret.intrinsic("check_does_not_have_tag", {v.b}, {v.R, v.p, int32_cst(1)});
    ret.bool_assert(v.b);
    if (tainted_body) {
      EXPECT_COUNTS(intra(cfg), 0, 1);
    } else {
      EXPECT_COUNTS(intra(cfg), 1, 0);
    }
  }
}

// (iv) Call boundary (top-down inter analysis, generic call
// transformer: forget outputs, meet). foo taints its output region
// R; bar only reads its input region Q; main has its own tainted
// region S. After the calls: S keeps tag 2 (the callee said nothing
// about it), R carries tag 1 (the caller had forgotten it), S does
// not carry tag 1, and Q stays clean.
BOOST_AUTO_TEST_CASE(call_return_keeps_both_sides) {
  variable_factory_t vfac;
  crab::tag_manager as;
  z_var x(vfac["x"], crab::INT_TYPE, 32), t(vfac["t"], crab::INT_TYPE, 32);
  z_var R(vfac["R"], crab::REG_INT_TYPE, 32), S(vfac["S"], crab::REG_INT_TYPE, 32),
      Q(vfac["Q"], crab::REG_INT_TYPE, 32);
  z_var r(vfac["r"], crab::REF_TYPE), s(vfac["s"], crab::REF_TYPE),
      q(vfac["q"], crab::REF_TYPE), q2(vfac["q2"], crab::REF_TYPE);
  z_var b1(vfac["b1"], crab::BOOL_TYPE), b2(vfac["b2"], crab::BOOL_TYPE),
      b3(vfac["b3"], crab::BOOL_TYPE), b4(vfac["b4"], crab::BOOL_TYPE);

  function_decl<z_number, varname_t> dfoo("foo", {x}, {R});
  z_cfg_t foo("entry", "exit", dfoo);
  auto &fe = foo.insert("entry");
  auto &fx = foo.insert("exit");
  fe >> fx;
  fe.region_init(R);
  fe.make_ref(r, R, int32_cst(4), as.mk_tag());
  fe.intrinsic("add_tag", {}, {R, r, int32_cst(1)});

  function_decl<z_number, varname_t> dbar("bar", {Q}, {});
  z_cfg_t bar("entry", "exit", dbar);
  auto &be = bar.insert("entry");
  auto &bx = bar.insert("exit");
  be >> bx;
  be.make_ref(q2, Q, int32_cst(4), as.mk_tag());
  be.load_from_ref(t, q2, Q);

  function_decl<z_number, varname_t> dmain("main", {}, {});
  z_cfg_t m("entry", "exit", dmain);
  auto &me = m.insert("entry");
  auto &mx = m.insert("exit");
  me >> mx;
  me.region_init(S);
  me.make_ref(s, S, int32_cst(4), as.mk_tag());
  me.intrinsic("add_tag", {}, {S, s, int32_cst(2)});
  me.region_init(Q);
  me.make_ref(q, Q, int32_cst(4), as.mk_tag());
  me.store_to_ref(q, Q, int32_cst(0));
  me.havoc(x);
  mx.callsite("foo", {R}, {x});
  mx.callsite("bar", {}, {Q});
  mx.intrinsic("check_does_not_have_tag", {b1}, {S, s, int32_cst(2)});
  mx.bool_assert(b1); // warning: main's own taint survives the call
  mx.intrinsic("check_does_not_have_tag", {b2}, {R, s, int32_cst(1)});
  mx.bool_assert(b2); // warning: the callee's taint reaches the caller
  mx.intrinsic("check_does_not_have_tag", {b3}, {S, s, int32_cst(1)});
  mx.bool_assert(b3); // safe: tag 1 never reached S
  mx.intrinsic("check_does_not_have_tag", {b4}, {Q, q, int32_cst(1)});
  mx.bool_assert(b4); // safe: bar only read Q

  EXPECT_COUNTS(inter({foo, bar, m}), 2, 2);
}

// A callee that creates objects in a region it returns ("new" region,
// passed as an output only) is called twice; between the calls the
// caller tags the object it got from the first call. The second call
// must not erase that tag: the region collects the objects of both
// calls (the caller initialises it at entry; the second call's clean
// objects join it).
BOOST_AUTO_TEST_CASE(new_region_accumulates_across_calls) {
  variable_factory_t vfac;
  crab::tag_manager as;
  z_var R(vfac["R"], crab::REG_INT_TYPE, 32);
  z_var r(vfac["r"], crab::REF_TYPE), q(vfac["q"], crab::REF_TYPE);
  z_var b1(vfac["b1"], crab::BOOL_TYPE), b2(vfac["b2"], crab::BOOL_TYPE);

  function_decl<z_number, varname_t> dfoo("foo", {}, {R});
  z_cfg_t foo("entry", "exit", dfoo);
  auto &fe = foo.insert("entry");
  auto &fx = foo.insert("exit");
  fe >> fx;
  fe.region_init(R);
  fe.make_ref(r, R, int32_cst(4), as.mk_tag());
  fe.store_to_ref(r, R, int32_cst(0));

  function_decl<z_number, varname_t> dmain("main", {}, {});
  z_cfg_t m("entry", "exit", dmain);
  auto &me = m.insert("entry");
  auto &mx = m.insert("exit");
  me >> mx;
  me.region_init(R);                 // empty at entry (clam does this)
  me.callsite("foo", {R}, {});
  me.make_ref(q, R, int32_cst(4), as.mk_tag());
  me.intrinsic("add_tag", {}, {R, q, int32_cst(2)}); // the caller tags an object
  mx.callsite("foo", {R}, {});       // clean objects join R
  mx.intrinsic("check_does_not_have_tag", {b1}, {R, q, int32_cst(2)});
  mx.bool_assert(b1);                // warning: tag 2 survives
  mx.intrinsic("check_does_not_have_tag", {b2}, {R, q, int32_cst(1)});
  mx.bool_assert(b2);                // safe: nobody used tag 1
  EXPECT_COUNTS(inter({foo, m}), 1, 1);
}

// A reference returned by a callee carries the callee's tags for it
// (int_to_ref of a tainted integer inside foo); storing it into a
// region of references in the caller taints that region.
BOOST_AUTO_TEST_CASE(reference_output_carries_callee_tags) {
  variable_factory_t vfac;
  crab::tag_manager as;
  z_var RI(vfac["RI"], crab::REG_INT_TYPE, 32), pi(vfac["pi"], crab::REF_TYPE),
      iv(vfac["iv"], crab::INT_TYPE, 32), RT(vfac["RT"], crab::REG_INT_TYPE, 32),
      ret(vfac["ret"], crab::REF_TYPE), RP(vfac["RP"], crab::REG_REF_TYPE),
      pp(vfac["pp"], crab::REF_TYPE), q(vfac["q"], crab::REF_TYPE),
      b(vfac["b"], crab::BOOL_TYPE);

  function_decl<z_number, varname_t> dfoo("foo", {}, {ret});
  z_cfg_t foo("entry", "exit", dfoo);
  auto &fe = foo.insert("entry");
  auto &fx = foo.insert("exit");
  fe >> fx;
  fe.region_init(RI);
  fe.make_ref(pi, RI, int32_cst(4), as.mk_tag());
  fe.intrinsic("add_tag", {}, {RI, pi, int32_cst(1)});
  fe.load_from_ref(iv, pi, RI);
  fe.region_init(RT);
  fe.int_to_ref(iv, RT, ret);            // a tainted reference

  function_decl<z_number, varname_t> dmain("main", {}, {});
  z_cfg_t m("entry", "exit", dmain);
  auto &me = m.insert("entry");
  auto &mx = m.insert("exit");
  me >> mx;
  me.region_init(RP);
  me.make_ref(pp, RP, int32_cst(8), as.mk_tag());
  me.callsite("foo", {q}, {});
  mx.store_to_ref(pp, RP, q);
  mx.intrinsic("check_does_not_have_tag", {b}, {RP, pp, int32_cst(1)});
  mx.bool_assert(b);                     // warning
  EXPECT_COUNTS(inter({foo, m}), 0, 1);
}

// (v) havoc: the statement "x := nondet" defines x. Without the
// closed-world policy x is unknown and a value stored from it makes
// the region unknown; with region.tag_havoc_clean it is clean. A
// havoc of a region (its contents) is unknown under both.
BOOST_AUTO_TEST_CASE(havoc_policy) {
  for (int clean = 0; clean < 2; ++clean) {
    params_fixture::set_havoc_clean(clean);
    {
      vars v;
      z_cfg_t cfg("entry", "ret");
      auto &entry = cfg.insert("entry");
      auto &ret = cfg.insert("ret");
      entry >> ret;
      entry.region_init(v.R);
      entry.make_ref(v.p, v.R, int32_cst(4), v.as.mk_tag());
      entry.havoc(v.x);
      entry.store_to_ref(v.p, v.R, v.x);
      ret.intrinsic("check_does_not_have_tag", {v.b}, {v.R, v.p, int32_cst(1)});
      ret.bool_assert(v.b);
      if (clean) {
        EXPECT_COUNTS(intra(cfg), 1, 0);
      } else {
        EXPECT_COUNTS(intra(cfg), 0, 1);
      }
    }
    {
      vars v;
      z_cfg_t cfg("entry", "ret");
      auto &entry = cfg.insert("entry");
      auto &ret = cfg.insert("ret");
      entry >> ret;
      entry.region_init(v.R);
      entry.make_ref(v.p, v.R, int32_cst(4), v.as.mk_tag());
      entry.havoc(v.R);
      ret.intrinsic("check_does_not_have_tag", {v.b}, {v.R, v.p, int32_cst(1)});
      ret.bool_assert(v.b);
      EXPECT_COUNTS(intra(cfg), 0, 1);
    }
  }
  params_fixture::set_havoc_clean(false);
}

// (vi) Sanitiser: remove_tag on a single object clears the tag; on an
// unknown region there is nothing to remove from and it stays unknown.
BOOST_AUTO_TEST_CASE(sanitiser) {
  {
    vars v;
    z_cfg_t cfg("entry", "ret");
    auto &entry = cfg.insert("entry");
    auto &ret = cfg.insert("ret");
    entry >> ret;
    entry.region_init(v.R);
    entry.make_ref(v.p, v.R, int32_cst(4), v.as.mk_tag());
    entry.intrinsic("add_tag", {}, {v.R, v.p, int32_cst(1)});
    entry.intrinsic("remove_tag", {}, {v.R, v.p, int32_cst(1)});
    ret.intrinsic("check_does_not_have_tag", {v.b}, {v.R, v.p, int32_cst(1)});
    ret.bool_assert(v.b);
    EXPECT_COUNTS(intra(cfg), 1, 0);
  }
  {
    vars v;
    z_cfg_t cfg("entry", "ret");
    auto &entry = cfg.insert("entry");
    auto &ret = cfg.insert("ret");
    entry >> ret;
    entry.region_init(v.R);
    entry.make_ref(v.p, v.R, int32_cst(4), v.as.mk_tag());
    entry.havoc(v.R);
    entry.intrinsic("remove_tag", {}, {v.R, v.p, int32_cst(1)});
    ret.intrinsic("check_does_not_have_tag", {v.b}, {v.R, v.p, int32_cst(1)});
    ret.bool_assert(v.b);
    EXPECT_COUNTS(intra(cfg), 0, 1);
  }
}

// (vii) Statements that define a value must write its tags explicitly
// (with "absent = unknown" nothing is clean by default).

// A fresh reference is clean: storing it into a region of references
// leaves that region clean.
BOOST_AUTO_TEST_CASE(make_ref_is_clean) {
  vars v;
  z_var RP(v.vfac["RP"], crab::REG_REF_TYPE), pp(v.vfac["pp"], crab::REF_TYPE);
  z_cfg_t cfg("entry", "ret");
  auto &entry = cfg.insert("entry");
  auto &ret = cfg.insert("ret");
  entry >> ret;
  entry.region_init(v.R);
  entry.make_ref(v.p, v.R, int32_cst(4), v.as.mk_tag());
  entry.region_init(RP);
  entry.make_ref(pp, RP, int32_cst(8), v.as.mk_tag());
  entry.store_to_ref(pp, RP, v.p);
  ret.intrinsic("check_does_not_have_tag", {v.b}, {RP, pp, int32_cst(1)});
  ret.bool_assert(v.b);
  EXPECT_COUNTS(intra(cfg), 1, 0);
}

// gep defines its result with exactly the tags of its base pointer:
// clean base -> clean, tainted base (via int_to_ref of a tainted
// integer) -> tainted.
BOOST_AUTO_TEST_CASE(gep_carries_base_tags) {
  for (int tainted = 0; tainted < 2; ++tainted) {
    vars v;
    z_var RP(v.vfac["RP"], crab::REG_REF_TYPE), pp(v.vfac["pp"], crab::REF_TYPE),
        p2(v.vfac["p2"], crab::REF_TYPE), RI(v.vfac["RI"], crab::REG_INT_TYPE, 32),
        pi(v.vfac["pi"], crab::REF_TYPE);
    z_cfg_t cfg("entry", "ret");
    auto &entry = cfg.insert("entry");
    auto &ret = cfg.insert("ret");
    entry >> ret;
    entry.region_init(v.R);
    entry.region_init(RP);
    entry.make_ref(pp, RP, int32_cst(8), v.as.mk_tag());
    if (tainted) {
      entry.region_init(RI);
      entry.make_ref(pi, RI, int32_cst(4), v.as.mk_tag());
      entry.intrinsic("add_tag", {}, {RI, pi, int32_cst(1)});
      entry.load_from_ref(v.i, pi, RI);
      entry.int_to_ref(v.i, v.R, p2);
    } else {
      entry.make_ref(p2, v.R, int32_cst(8), v.as.mk_tag());
    }
    entry.gep_ref(v.q, v.R, p2, v.R, 4);
    entry.store_to_ref(pp, RP, v.q);
    ret.intrinsic("check_does_not_have_tag", {v.b}, {RP, pp, int32_cst(1)});
    ret.bool_assert(v.b);
    if (tainted) {
      EXPECT_COUNTS(intra(cfg), 0, 1);
    } else {
      EXPECT_COUNTS(intra(cfg), 1, 0);
    }
  }
}

// region_copy / region_cast define their destination with exactly the
// tags of the source: a clean source gives a clean destination (not
// unknown), a tainted source a tainted one.
BOOST_AUTO_TEST_CASE(region_copy_and_cast_are_definitions) {
  for (int tainted = 0; tainted < 2; ++tainted) {
    {
      vars v;
      z_cfg_t cfg("entry", "ret");
      auto &entry = cfg.insert("entry");
      auto &ret = cfg.insert("ret");
      entry >> ret;
      entry.region_init(v.R);
      entry.make_ref(v.p, v.R, int32_cst(4), v.as.mk_tag());
      if (tainted) {
        entry.intrinsic("add_tag", {}, {v.R, v.p, int32_cst(1)});
      }
      entry.region_copy(v.S, v.R);
      ret.intrinsic("check_does_not_have_tag", {v.b}, {v.S, v.p, int32_cst(1)});
      ret.bool_assert(v.b);
      if (tainted) {
        EXPECT_COUNTS(intra(cfg), 0, 1);
      } else {
        EXPECT_COUNTS(intra(cfg), 1, 0);
      }
    }
    {
      vars v;
      z_var U(v.vfac["U"], crab::REG_UNKNOWN_TYPE);
      z_cfg_t cfg("entry", "ret");
      auto &entry = cfg.insert("entry");
      auto &ret = cfg.insert("ret");
      entry >> ret;
      entry.region_init(v.R);
      entry.make_ref(v.p, v.R, int32_cst(4), v.as.mk_tag());
      if (tainted) {
        entry.intrinsic("add_tag", {}, {v.R, v.p, int32_cst(1)});
      }
      entry.region_cast(v.R, U);
      ret.intrinsic("check_does_not_have_tag", {v.b}, {U, v.p, int32_cst(1)});
      ret.bool_assert(v.b);
      if (tainted) {
        EXPECT_COUNTS(intra(cfg), 0, 1);
      } else {
        EXPECT_COUNTS(intra(cfg), 1, 0);
      }
    }
  }
}

// Casting into a typed region merges the incoming objects' tags with
// those the region already holds (objects created by a callee join the
// caller's region at every call site): tag 2 held by Z survives a cast
// of a region carrying tag 1, and a clean incoming region does not
// clean Z.
BOOST_AUTO_TEST_CASE(region_cast_into_typed_region_is_a_union) {
  for (int incoming_tainted = 0; incoming_tainted < 2; ++incoming_tainted) {
    vars v;
    z_var U(v.vfac["U"], crab::REG_UNKNOWN_TYPE), Z(v.vfac["Z"], crab::REG_INT_TYPE, 32),
        z(v.vfac["z"], crab::REF_TYPE);
    z_cfg_t cfg("entry", "ret");
    auto &entry = cfg.insert("entry");
    auto &ret = cfg.insert("ret");
    entry >> ret;
    entry.region_init(Z);
    entry.make_ref(z, Z, int32_cst(4), v.as.mk_tag());
    entry.intrinsic("add_tag", {}, {Z, z, int32_cst(2)});
    entry.region_init(v.R);
    entry.make_ref(v.p, v.R, int32_cst(4), v.as.mk_tag());
    if (incoming_tainted) {
      entry.intrinsic("add_tag", {}, {v.R, v.p, int32_cst(1)});
    }
    entry.region_cast(v.R, U);   // view of R
    entry.region_cast(U, Z);     // R's objects join Z
    ret.intrinsic("check_does_not_have_tag", {v.b}, {Z, z, int32_cst(2)});
    ret.bool_assert(v.b);        // warning: Z still holds tag 2
    ret.intrinsic("check_does_not_have_tag", {v.b2}, {Z, z, int32_cst(1)});
    ret.bool_assert(v.b2);       // tag 1 only if it came in
    if (incoming_tainted) {
      EXPECT_COUNTS(intra(cfg), 0, 2);
    } else {
      EXPECT_COUNTS(intra(cfg), 1, 1);
    }
  }
}

// A boolean defined by a reference comparison whose operands have
// disjoint allocation sites (the region domain's shortcut) is still
// tagged from its operands: two clean references give a clean boolean.
BOOST_AUTO_TEST_CASE(bool_from_ref_comparison_is_tagged) {
  vars v;
  z_var RB(v.vfac["RB"], crab::REG_BOOL_TYPE), pb(v.vfac["pb"], crab::REF_TYPE);
  z_cfg_t cfg("entry", "ret");
  auto &entry = cfg.insert("entry");
  auto &ret = cfg.insert("ret");
  entry >> ret;
  entry.region_init(v.R);
  entry.make_ref(v.p, v.R, int32_cst(4), v.as.mk_tag());
  entry.make_ref(v.q, v.R, int32_cst(4), v.as.mk_tag());
  entry.bool_assign(v.b2, z_ref_cst_t::mk_eq(v.p, v.q));
  entry.region_init(RB);
  entry.make_ref(pb, RB, int32_cst(1), v.as.mk_tag());
  entry.store_to_ref(pb, RB, v.b2);
  ret.intrinsic("check_does_not_have_tag", {v.b}, {RB, pb, int32_cst(1)});
  ret.bool_assert(v.b);
  EXPECT_COUNTS(intra(cfg), 1, 0);
}
