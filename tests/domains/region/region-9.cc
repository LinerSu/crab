#include "../../common.hpp"
#include "../../program_options.hpp"

#include <crab/analysis/fwd_analyzer.hpp>
#include <crab/checkers/assertion.hpp>
#include <crab/checkers/base_property.hpp>
#include <crab/checkers/checker.hpp>

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/*
 * Tests for the tag (taint) analysis: interaction between the tag
 * intrinsics (add_tag, move_tag, remove_tag) and the strong/weak
 * update policy of the region domain.
 *
 * Invariant: the tags attached to a region are part of its content.
 * Therefore:
 *
 * 1. Any intrinsic that modifies the tags of a region marks the
 *    region as possibly written, exactly like a store does, so that
 *    a later store to a region that stands for more than one
 *    concrete cell is a weak update that keeps the existing tags.
 *
 * (sequence regions are a separate follow-up: mark_sequence_region is not part of this test)
 *    for an unbounded number of cells and is never strongly updated,
 *    regardless of its reference count.
 *
 * 3. A genuine singleton region (one reference, not a sequence) is
 *    still strongly updated: overwriting it with an untainted value
 *    legitimately drops its tags.
 *
 * Unlike most tests, this one also verifies programmatically the
 * number of safe/warning checks and returns a non-zero exit code if
 * they differ from the expected ones.
 */

static z_var_or_cst_t int32_cst(int n) {
  return z_var_or_cst_t(z_number(n), crab::variable_type(crab::INT_TYPE, 32));
}

// Two references (p and q) to the same region R: R stands for two
// cells. After add_tag(R, p, TAG_1), the store through q must be a
// weak update because only one of the two cells is overwritten.
z_cfg_t *cfg1(variable_factory_t &vfac) {
  /*
    int *p = malloc(4);  // region R
    int *q = malloc(4);  // region R
    add_tag(p, TAG_1);
    *q = 0;
    assert(check_does_not_have_tag(p, TAG_1)); // EXPECTED: FAIL
  */
  z_var p(vfac["p"], crab::REF_TYPE);
  z_var q(vfac["q"], crab::REF_TYPE);
  z_var b1(vfac["b1"], crab::BOOL_TYPE);
  z_var R(vfac["R"], crab::REG_INT_TYPE, 32);
  crab::tag_manager as_man;
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &ret = cfg->insert("ret");
  entry.add_succ(ret);

  entry.region_init(R);
  entry.make_ref(p, R, int32_cst(4), as_man.mk_tag());
  entry.make_ref(q, R, int32_cst(4), as_man.mk_tag());
  entry.intrinsic("add_tag", {}, {R, p, int32_cst(1)});
  entry.store_to_ref(q, R, int32_cst(0));
  ret.intrinsic("check_does_not_have_tag", {b1}, {R, p, int32_cst(1)});
  // EXPECTED: FAIL (the cell pointed by p is still tagged)
  ret.bool_assert(b1);
  return cfg;
}

// The memcpy scenario: the tags of a tagged region S are moved
// (memcpy) into a region D that has two references (d1 and d1+4).
// The store through the second reference must be a weak update so
// that the tags moved into D survive.
z_cfg_t *cfg2(variable_factory_t &vfac) {
  /*
    char *s = malloc(4);        // region S
    add_tag(s, TAG_1);          // e.g., untrusted input
    char *d1 = malloc(8);       // region D
    char *d2 = d1 + 4;          // region D
    memcpy(d1, s, 4);           // move_tag(S, s, D, d1)
    *d2 = 0;                    // e.g., null-terminate the buffer
    assert(check_does_not_have_tag(d1, TAG_1)); // EXPECTED: FAIL
  */
  z_var s(vfac["s"], crab::REF_TYPE);
  z_var d1(vfac["d1"], crab::REF_TYPE);
  z_var d2(vfac["d2"], crab::REF_TYPE);
  z_var b1(vfac["b1"], crab::BOOL_TYPE);
  z_var S(vfac["S"], crab::REG_INT_TYPE, 32);
  z_var D(vfac["D"], crab::REG_INT_TYPE, 32);
  crab::tag_manager as_man;
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &ret = cfg->insert("ret");
  entry.add_succ(ret);

  entry.region_init(S);
  entry.region_init(D);
  entry.make_ref(s, S, int32_cst(4), as_man.mk_tag());
  entry.intrinsic("add_tag", {}, {S, s, int32_cst(1)});
  entry.make_ref(d1, D, int32_cst(8), as_man.mk_tag());
  entry.gep_ref(d2, D, d1, D, 4);
  entry.intrinsic("move_tag", {}, {S, s, D, d1});
  entry.store_to_ref(d2, D, int32_cst(0));
  ret.intrinsic("check_does_not_have_tag", {b1}, {D, d1, int32_cst(1)});
  // EXPECTED: FAIL (the cells pointed by d1 are still tagged)
  ret.bool_assert(b1);
  return cfg;
}

// Genuine singleton region: one reference, not a sequence. The store
// of an untainted constant is a strong update that legitimately
// drops the tag: after the store, the only cell of R contains 0.
z_cfg_t *cfg3(variable_factory_t &vfac) {
  /*
    int *p = malloc(4);  // region R
    add_tag(p, TAG_1);
    *p = 0;
    assert(check_does_not_have_tag(p, TAG_1)); // EXPECTED: OK
  */
  z_var p(vfac["p"], crab::REF_TYPE);
  z_var b1(vfac["b1"], crab::BOOL_TYPE);
  z_var R(vfac["R"], crab::REG_INT_TYPE, 32);
  crab::tag_manager as_man;
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &ret = cfg->insert("ret");
  entry.add_succ(ret);

  entry.region_init(R);
  entry.make_ref(p, R, int32_cst(4), as_man.mk_tag());
  entry.intrinsic("add_tag", {}, {R, p, int32_cst(1)});
  entry.store_to_ref(p, R, int32_cst(0));
  ret.intrinsic("check_does_not_have_tag", {b1}, {R, p, int32_cst(1)});
  // EXPECTED: OK (strong update of a singleton with an untainted value)
  ret.bool_assert(b1);
  return cfg;
}

// Sequence region with a single reference: even if the reference
// count is one, R stands for an unbounded number of cells (e.g., an
// array) so the store through p modifies only one of them and must
// be a weak update.

// Run the analysis and the assertion checker and compare the number
// of safe and warning checks with the expected ones.
template <typename Dom>
static bool verify(z_cfg_t *cfg, Dom init, unsigned expected_safe,
                   unsigned expected_warnings) {
  using analyzer_t = crab::analyzer::intra_fwd_analyzer<z_cfg_ref_t, Dom>;
  using checker_t = crab::checker::intra_checker<analyzer_t>;
  using assert_checker_t = crab::checker::assert_property_checker<analyzer_t>;

  crab::fixpoint_parameters fixpo_params;
  analyzer_t a(*cfg, init.make_top(), nullptr, fixpo_params);
  typename analyzer_t::assumption_map_t assumptions;
  a.run(cfg->entry(), init, assumptions);
  typename checker_t::prop_checker_ptr prop(new assert_checker_t(0));
  checker_t checker(a, {prop});
  checker.run();
  crab::checker::checks_db db = checker.get_all_checks();
  return (db.get_total_safe() == expected_safe &&
          db.get_total_warning() == expected_warnings &&
          db.get_total_error() == 0);
}

template <typename Dom>
static bool run_test(const char *name, z_cfg_t *cfg, Dom init,
                     unsigned expected_safe, unsigned expected_warnings,
                     bool stats_enabled) {
  crab::outs() << *cfg << "\n";
  run_and_check(cfg, cfg->entry(), init, false, 2, 2, 20, stats_enabled);
#ifdef USE_GENERIC_WRAPPER
  z_abs_domain_t init_wrapper(init);
  bool ok = verify(cfg, init_wrapper, expected_safe, expected_warnings);
#else
  bool ok = verify(cfg, init, expected_safe, expected_warnings);
#endif
  crab::outs() << name << ": expected " << expected_safe << " safe and "
               << expected_warnings
               << " warning checks: " << (ok ? "OK" : "FAILED") << "\n\n";
  return ok;
}

int main(int argc, char **argv) {
  region_domain_params p(true /*allocation_sites*/, true /*deallocation*/,
                         true /*tag_analysis*/, false /*is_dereferenceable*/,
                         true /*skip_unknown_regions*/);
  crab_domain_params_man::get().update_params(p);

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  variable_factory_t vfac;
  z_rgn_bool_int_t init;
  bool ok = true;

  {
    z_cfg_t *cfg = cfg1(vfac);
    ok = run_test("cfg1 (add_tag, store through another reference)", cfg,
                  init, 0, 1, stats_enabled) && ok;
    delete cfg;
  }
  {
    z_cfg_t *cfg = cfg2(vfac);
    ok = run_test("cfg2 (move_tag, store through another reference)", cfg,
                  init, 0, 1, stats_enabled) && ok;
    delete cfg;
  }
  {
    z_cfg_t *cfg = cfg3(vfac);
    ok = run_test("cfg3 (singleton: strong update drops the tag)", cfg, init,
                  1, 0, stats_enabled) && ok;
    delete cfg;
  }

  if (!ok) {
    crab::outs() << "region-9: some checks produced unexpected results\n";
    return 1;
  }
  return 0;
}
