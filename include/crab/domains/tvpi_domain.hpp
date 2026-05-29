/*******************************************************************************
 * Two-Variables-Per-Inequality (TVPI) abstract domain for Crab.
 *
 * Algorithm and original C++ implementation:
 *   Axel Simon, Andy King, and Jacob M. Howe.
 *   "Two Variables per Linear Inequality as an Abstract Domain."
 *   Higher-Order and Symbolic Computation, 2010.
 *   https://doi.org/10.1007/s10990-010-9062-8
 *
 *   Original TVPI library authored by Axel Simon <A.Simon@kent.ac.uk>.
 *   Source: https://github.com/axelsimon/tvpi  (GPL v2)
 *
 * Crab abstract domain adapter:
 *   Yusen Su <yusen.su@uwaterloo.ca>
 *   (with assistance from Claude Sonnet 4.6, Anthropic)
 *   2026
 *
 * License note:
 *   This file (tvpi_domain.hpp) is the Crab adapter and is part of the Crab
 *   project (NASA Open Source Agreement v1.3).  The underlying TVPI engine
 *   bundled in tvpi_impl/ and lib/tvpi_*.cpp is derived from the TVPI library
 *   by Axel Simon, which is distributed under the GNU General Public License
 *   version 2 (GPL v2).  Redistribution of the combined work must comply with
 *   the GPL v2; see tvpi_impl/COPYING for the full GPL v2 text.
 *
 * Supported operations:
 *   - All standard lattice operations (top/bottom, join, meet, widening,
 *     narrowing, inclusion check)
 *   - assign(x, c)                     constant assignment
 *   - assign(x, y)                     variable copy
 *   - assign(x, a*y + c)               affine assignment (one variable)
 *   - assign(x, a*y + b*z + c)         two-variable affine assignment
 *   - operator+=(linear constraints)   assume constraints (>=, <=, ==)
 *   - operator-=(v)                    forget variable
 *   - at(v) / operator[](v)            query interval of variable
 *   - to_linear_constraint_system()    extract constraints
 *   - forget, project, rename, expand
 *   - apply (arithmetic: +, -, *)      TVPI rewrite; /, % conservative forget
 *   - apply (bitwise)                  conservative forget
 *   - backward_assign / backward_apply not implemented
 *
 * Number type: integers only (DenseTvpi<true>, isZ=true).
 ******************************************************************************/

#pragma once

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/interval.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

// TVPI internal headers bundled under tvpi_impl/
#include <crab/domains/tvpi_impl/tvpi.hh>
#include <crab/domains/tvpi_impl/planar.hh>
#include <crab/domains/tvpi_impl/polyhedron.hh>
#include <crab/domains/tvpi_impl/interval.hh>

#include <boost/optional.hpp>
#include <algorithm>
#include <cassert>
#include <string>
#include <unordered_map>
#include <vector>

namespace crab {
namespace domains {

// ---------------------------------------------------------------------------
// Helper utilities for converting between Crab and TVPI number types
// ---------------------------------------------------------------------------
namespace tvpi_detail {

// z_number -> mpz_class
inline mpz_class to_mpz(const ikos::z_number &n) {
  mpz_class result;
  mpz_set(result.get_mpz_t(),
          const_cast<ikos::z_number &>(n).get_mpz_t());
  return result;
}

// mpz_class -> z_number (via string representation)
inline ikos::z_number from_mpz(const mpz_class &m) {
  return ikos::z_number(m.get_str());
}

// Tvpi::Interval<true> -> ikos::interval<z_number>
inline ikos::interval<ikos::z_number>
tvpi_to_crab_interval(const Tvpi::Interval<true> &ti) {
  using bound_t = ikos::bound<ikos::z_number>;
  using interval_t = ikos::interval<ikos::z_number>;

  bound_t lb = bound_t::minus_infinity();
  bound_t ub = bound_t::plus_infinity();

  if (ti.lowerIsFinite()) {
    mpq_class q = ti.getLower();
    mpz_class val;
    // Ceiling division: smallest integer >= q
    mpz_cdiv_q(val.get_mpz_t(), q.get_num().get_mpz_t(),
               q.get_den().get_mpz_t());
    lb = bound_t(from_mpz(val));
  }
  if (ti.upperIsFinite()) {
    mpq_class q = ti.getUpper();
    mpz_class val;
    // Floor division: largest integer <= q
    mpz_fdiv_q(val.get_mpz_t(), q.get_num().get_mpz_t(),
               q.get_den().get_mpz_t());
    ub = bound_t(from_mpz(val));
  }
  return interval_t(lb, ub);
}

// ikos::interval<z_number> -> Tvpi::Interval<true>
inline Tvpi::Interval<true>
crab_to_tvpi_interval(const ikos::interval<ikos::z_number> &ci) {
  Tvpi::Interval<true> result;
  if (ci.lb().is_finite()) {
    auto n = ci.lb().number();
    assert(n);
    mpq_class q(to_mpz(*n));
    result.updateLower(q);
  }
  if (ci.ub().is_finite()) {
    auto n = ci.ub().number();
    assert(n);
    mpq_class q(to_mpz(*n));
    result.updateUpper(q);
  }
  return result;
}

} // namespace tvpi_detail

// ---------------------------------------------------------------------------
// tvpi_domain: wrapper around DenseTvpi<true>
// ---------------------------------------------------------------------------

template <typename Number, typename VariableName>
class tvpi_domain final
    : public abstract_domain_api<tvpi_domain<Number, VariableName>> {
public:
  using tvpi_domain_t = tvpi_domain<Number, VariableName>;
  using abstract_domain_t = abstract_domain_api<tvpi_domain_t>;
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::interval_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_or_constant_t;
  using typename abstract_domain_t::variable_or_constant_vector_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using typename abstract_domain_t::varname_t;
  using number_t = Number;

private:
  using var_map_t = std::unordered_map<variable_t, TvpiVar>;
  using rev_map_t = std::vector<boost::optional<variable_t>>; // indexed by TvpiVar

  Tvpi::DenseTvpi<true> m_tvpi;  // underlying TVPI domain
  var_map_t m_var_map;            // variable_t -> TvpiVar index
  rev_map_t m_rev_map;            // TvpiVar index -> variable_t
  bool m_is_bottom;

  // -------------------------------------------------------------------
  // Variable management helpers
  // -------------------------------------------------------------------

  // Return TvpiVar for v (read-only, does not create)
  boost::optional<TvpiVar> get_var(const variable_t &v) const {
    auto it = m_var_map.find(v);
    if (it != m_var_map.end())
      return it->second;
    return boost::none;
  }

  // Return TvpiVar for v, creating an unbounded TVPI variable if absent
  TvpiVar get_or_create_var(const variable_t &v) {
    auto it = m_var_map.find(v);
    if (it != m_var_map.end())
      return it->second;

    TvpiVar idx = m_tvpi.createVariable();
    m_var_map[v] = idx;
    if (idx >= (TvpiVar)m_rev_map.size())
      m_rev_map.resize(idx + 1, boost::none);
    m_rev_map[idx] = v;
    return idx;
  }

  // -------------------------------------------------------------------
  // Number conversion helpers
  // -------------------------------------------------------------------

  static mpz_class to_mpz(const number_t &n) {
    return tvpi_detail::to_mpz(n);
  }
  static number_t from_mpz(const mpz_class &m) {
    return tvpi_detail::from_mpz(m);
  }
  static Tvpi::Interval<true> to_tvpi_iv(const interval_t &ci) {
    return tvpi_detail::crab_to_tvpi_interval(ci);
  }
  static interval_t from_tvpi_iv(const Tvpi::Interval<true> &ti) {
    return tvpi_detail::tvpi_to_crab_interval(ti);
  }

  // -------------------------------------------------------------------
  // Internal: add a single linear constraint to m_tvpi
  // -------------------------------------------------------------------

  // Add  sum_i c_i * v_i <= rhs  to the domain.
  // Returns false if unsatisfiable.
  bool add_inequality(const linear_expression_t &lhs_vars, const number_t &rhs) {
    std::vector<LinComponent> comps;
    for (auto it = lhs_vars.begin(), et = lhs_vars.end(); it != et; ++it) {
      if (it->first == number_t(0))
        continue;
      // DenseTvpi::approximateInequality uses LinComponent::variable as the
      // direct TvpiVar index (0-based) into the bounds[] array.
      TvpiVar idx = get_or_create_var(it->second);
      comps.push_back(LinComponent(to_mpz(it->first),
                                         static_cast<Variable>(idx)));
    }
    if (comps.empty())
      return true; // tautology (0 <= rhs assumed)

    std::sort(comps.begin(), comps.end());
    // approximateInequality takes: sum(comps) + constant <= 0
    // We want: sum(comps) <= rhs, i.e., sum(comps) + (-rhs) <= 0
    mpz_class neg_rhs = -to_mpz(rhs);
    Result res = m_tvpi.approximateInequality(comps, neg_rhs, false);
    return (res != resUnsatisfiable);
  }

  void add_constraint(const linear_constraint_t &cst) {
    if (cst.is_disequation())
      return; // TVPI has no disequality support; skip conservatively

    // cst.expression() = sum c_i v_i + k   with constraint kind  e <= 0 / == 0
    // Rewrite as: sum c_i v_i <= -k
    const linear_expression_t &e = cst.expression();
    number_t minus_const = -e.constant();
    // Strip constant from expression
    linear_expression_t vars_only = e - e.constant();

    if (cst.is_inequality()) {
      if (!add_inequality(vars_only, minus_const))
        m_is_bottom = true;
    } else if (cst.is_equality()) {
      // sum c_i v_i <= -k   AND   sum -c_i v_i <= k
      if (!add_inequality(vars_only, minus_const)) {
        m_is_bottom = true;
        return;
      }
      linear_expression_t neg_vars = vars_only * number_t(-1);
      if (!add_inequality(neg_vars, -minus_const))
        m_is_bottom = true;
    }
  }

  // -------------------------------------------------------------------
  // Internal: assignment helpers
  // -------------------------------------------------------------------

  // x := a*y + c (single variable linear expression)
  void assign_affine_one_var(const variable_t &x, number_t a,
                              const variable_t &y, number_t c) {
    auto opt_y = get_var(y);
    if (!opt_y) {
      // y not known; forget x conservatively
      *this -= x;
      return;
    }
    TvpiVar x_idx = get_or_create_var(x);
    // updateVariable(a1, x1, a2, x2, cMin, cMax): a1*x1 = a2*x2 + [cMin..cMax]
    // Here: 1*x = a*y + c  =>  a1=1, x1=x, a2=a, x2=y, c=c
    mpz_class a_mpz = to_mpz(a);
    mpz_class c_mpz = to_mpz(c);
    m_tvpi.updateVariable(mpz_class(1), x_idx, a_mpz, *opt_y, c_mpz, c_mpz);
  }

  // Compute interval of linear expression e using known variable intervals
  interval_t eval_expr_interval(const linear_expression_t &e) const {
    interval_t result(e.constant());
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      result = result + interval_t(it->first) * at(it->second);
    }
    return result;
  }

  // Apply interval bounds for variable x
  void apply_interval(const variable_t &x, const interval_t &xi) {
    if (xi.is_bottom()) {
      m_is_bottom = true;
      return;
    }
    TvpiVar idx = get_or_create_var(x);
    Tvpi::Interval<true> tvpi_i = to_tvpi_iv(xi);
    bool sat = m_tvpi.intersectBound(idx, tvpi_i);
    if (!sat)
      m_is_bottom = true;
  }

public:
  // -------------------------------------------------------------------
  // Constructors / destructor / assignment
  // -------------------------------------------------------------------

  tvpi_domain() : m_tvpi(), m_is_bottom(false) {}

  tvpi_domain(const tvpi_domain_t &o) = default;
  tvpi_domain(tvpi_domain_t &&o) = default;
  tvpi_domain_t &operator=(const tvpi_domain_t &o) = default;
  tvpi_domain_t &operator=(tvpi_domain_t &&o) = default;

  // -------------------------------------------------------------------
  // Lattice: top / bottom
  // -------------------------------------------------------------------

  tvpi_domain_t make_top() const override { return tvpi_domain_t(); }

  tvpi_domain_t make_bottom() const override {
    tvpi_domain_t res;
    res.m_is_bottom = true;
    return res;
  }

  void set_to_top() override {
    m_tvpi = Tvpi::DenseTvpi<true>();
    m_var_map.clear();
    m_rev_map.clear();
    m_is_bottom = false;
  }

  void set_to_bottom() override { m_is_bottom = true; }

  bool is_bottom() const override { return m_is_bottom; }

  bool is_top() const override {
    return !m_is_bottom && m_var_map.empty();
  }

  // -------------------------------------------------------------------
  // Internal: build a restricted copy with only the common variables
  // shared with `other`, returning the copy and two perm arrays.
  //
  // The TVPI join/includes/widen functions require that every variable
  // in *this has a matching variable in other (no Tvpi::invalidTvpiVar allowed
  // as a hole). This helper projects both sides to the common variables.
  //
  // Returns false if either side has no variables after intersection.
  // -------------------------------------------------------------------

  struct CommonVarInfo {
    tvpi_domain_t left;   // restricted copy of *this
    tvpi_domain_t right;  // restricted copy of other
    // perm[i] = right's TvpiVar for left's TvpiVar i
    std::vector<TvpiVar> perm;
    std::vector<Mult> mults;
  };

  CommonVarInfo restrict_to_common(const tvpi_domain_t &other) const {
    CommonVarInfo info;

    // Find common variables
    std::vector<variable_t> common_vars;
    for (auto &kv : m_var_map) {
      if (other.m_var_map.count(kv.first))
        common_vars.push_back(kv.first);
    }

    if (common_vars.empty()) {
      // No common variables: both become top after projection
      info.left = make_top();
      info.right = make_top();
      return info;
    }

    // Build left restricted
    info.left = *this;
    {
      std::vector<TvpiVar> seq;
      for (auto &v : common_vars) {
        auto opt = info.left.get_var(v);
        if (opt) seq.push_back(*opt);
      }
      std::sort(seq.begin(), seq.end());
      info.left.m_tvpi.projectOnto(seq.size(), seq.data());
      // Rebuild left maps
      var_map_t nm; rev_map_t nr(seq.size(), boost::none);
      for (size_t i = 0; i < seq.size(); i++) {
        TvpiVar old_i = seq[i];
        assert(info.left.m_rev_map[old_i]);
        const variable_t &v = *info.left.m_rev_map[old_i];
        nm[v] = static_cast<TvpiVar>(i);
        nr[i] = v;
      }
      info.left.m_var_map = std::move(nm);
      info.left.m_rev_map = std::move(nr);
    }

    // Build right restricted
    info.right = other;
    {
      std::vector<TvpiVar> seq;
      for (auto &v : common_vars) {
        auto opt = info.right.get_var(v);
        if (opt) seq.push_back(*opt);
      }
      std::sort(seq.begin(), seq.end());
      info.right.m_tvpi.projectOnto(seq.size(), seq.data());
      var_map_t nm; rev_map_t nr(seq.size(), boost::none);
      for (size_t i = 0; i < seq.size(); i++) {
        TvpiVar old_i = seq[i];
        assert(info.right.m_rev_map[old_i]);
        const variable_t &v = *info.right.m_rev_map[old_i];
        nm[v] = static_cast<TvpiVar>(i);
        nr[i] = v;
      }
      info.right.m_var_map = std::move(nm);
      info.right.m_rev_map = std::move(nr);
    }

    // Build perm: left TvpiVar -> right TvpiVar
    TvpiVar n = static_cast<TvpiVar>(info.left.m_tvpi.size());
    info.perm.assign(n, Tvpi::invalidTvpiVar);
    info.mults.assign(
        info.right.m_tvpi.size() > 0 ? info.right.m_tvpi.size() : 1, 0);
    for (auto &kv : info.left.m_var_map) {
      auto opt = info.right.get_var(kv.first);
      if (opt)
        info.perm[kv.second] = *opt;
    }
    return info;
  }

  // -------------------------------------------------------------------
  // Lattice: inclusion (<=)
  // -------------------------------------------------------------------

  bool operator<=(const tvpi_domain_t &other) const override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");

    if (is_bottom() || other.is_top())
      return true;
    if (is_top() || other.is_bottom())
      return false;

    // *this <= other means: concrete_states(*this) ⊆ concrete_states(other).
    // In polyhedra terms: other.includes(*this) should be true.
    //
    // After restrict_to_common: info.perm maps left (= *this) indices to
    // right (= other) indices. The TVPI includes() takes perm mapping
    // "this" indices to "other" indices (perm[i] = index in 'other' of
    // variable i in 'this'). So we call:
    //   info.right.includes(info.left, inv_perm)
    // where inv_perm maps right->left, or equivalently:
    //   info.right.includes(info.left, perm_right_to_left)
    //
    // Since our common-variable projection ensures perm is identity,
    // we can simply call info.right.m_tvpi.includes(info.left.m_tvpi, perm)
    // with the identity perm (which is the same either way).
    CommonVarInfo info = restrict_to_common(other);
    if (info.left.is_top() || info.right.is_top())
      return true; // no common variables, so trivially included

    TvpiVar n = static_cast<TvpiVar>(info.right.m_tvpi.size());
    if (n == 0) return true;

    // info.perm[i] = i (identity due to common-var ordering), so passing
    // it to right.includes(left, perm) checks right ⊇ left, i.e., *this ⊆ other.
    return const_cast<Tvpi::DenseTvpi<true> &>(info.right.m_tvpi).includes(
        const_cast<Tvpi::DenseTvpi<true> &>(info.left.m_tvpi),
        info.perm.data());
  }

  // -------------------------------------------------------------------
  // Lattice: join (|)
  // -------------------------------------------------------------------

  void operator|=(const tvpi_domain_t &other) override {
    *this = *this | other;
  }

  tvpi_domain_t operator|(const tvpi_domain_t &other) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    if (is_bottom())
      return other;
    if (other.is_bottom())
      return *this;
    if (is_top() || other.is_top())
      return make_top();

    // Join is defined on the common variables only.
    // Variables exclusive to one side become unconstrained in the result
    // (since the other side knows nothing about them = top).
    CommonVarInfo info = restrict_to_common(other);
    if (info.left.is_top() || info.right.is_top())
      return make_top();

    TvpiVar n = static_cast<TvpiVar>(info.left.m_tvpi.size());
    if (n == 0) return make_top();

    tvpi_domain_t result(info.left);
    result.m_tvpi.join(
        const_cast<Tvpi::DenseTvpi<true> &>(info.right.m_tvpi),
        info.perm.data(), info.mults.data());
    return result;
  }

  // -------------------------------------------------------------------
  // Lattice: meet (&)
  // -------------------------------------------------------------------

  void operator&=(const tvpi_domain_t &other) override {
    *this = *this & other;
  }

  tvpi_domain_t operator&(const tvpi_domain_t &other) const override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    if (is_bottom() || other.is_bottom())
      return make_bottom();
    if (is_top())
      return other;
    if (other.is_top())
      return *this;

    // Meet: copy this, then add all constraints from other
    tvpi_domain_t result(*this);
    linear_constraint_system_t csts = other.to_linear_constraint_system();
    result += csts;
    return result;
  }

  // -------------------------------------------------------------------
  // Lattice: widening (||)
  // -------------------------------------------------------------------

  tvpi_domain_t operator||(const tvpi_domain_t &other) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");

    if (is_bottom())
      return other;
    if (other.is_bottom())
      return *this;

    // Widening is defined on the common variables (same reasoning as join).
    CommonVarInfo info = restrict_to_common(other);
    if (info.left.is_top() || info.right.is_top())
      return make_top();

    TvpiVar n = static_cast<TvpiVar>(info.left.m_tvpi.size());
    if (n == 0) return make_top();

    tvpi_domain_t result(info.left);
    // extrapolate=0 => standard widening
    result.m_tvpi.widen(
        const_cast<Tvpi::DenseTvpi<true> &>(info.right.m_tvpi),
        info.perm.data(), mpz_class(0));
    return result;
  }

  tvpi_domain_t widening_thresholds(const tvpi_domain_t &other,
                                    const thresholds<number_t> &) const override {
    return *this || other;
  }

  // -------------------------------------------------------------------
  // Lattice: narrowing (&&)
  // -------------------------------------------------------------------

  tvpi_domain_t operator&&(const tvpi_domain_t &other) const override {
    // Narrowing: meet is a sound (conservative) narrowing
    return *this & other;
  }

  // -------------------------------------------------------------------
  // Variable forget / interval query
  // -------------------------------------------------------------------

  void operator-=(const variable_t &v) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");

    if (is_bottom())
      return;

    auto opt = get_var(v);
    if (!opt)
      return;

    TvpiVar idx = *opt;

    // Build projection sequence: all variables except idx, sorted
    std::vector<TvpiVar> seq;
    seq.reserve(m_tvpi.size() > 1 ? m_tvpi.size() - 1 : 0);
    for (TvpiVar i = 0; i < (TvpiVar)m_tvpi.size(); i++) {
      if (i != idx)
        seq.push_back(i);
    }

    if (seq.empty()) {
      // Only variable: reset to top
      set_to_top();
      return;
    }

    m_tvpi.projectOnto(seq.size(), seq.data());

    // Rebuild maps: indices > idx shift down by one
    var_map_t new_map;
    rev_map_t new_rev(seq.size(), boost::none);
    for (auto &kv : m_var_map) {
      if (kv.first == v)
        continue;
      TvpiVar old_i = kv.second;
      TvpiVar new_i = (old_i > idx) ? old_i - 1 : old_i;
      new_map[kv.first] = new_i;
      new_rev[new_i] = kv.first;
    }
    m_var_map = std::move(new_map);
    m_rev_map = std::move(new_rev);
  }

  interval_t operator[](const variable_t &v) override { return at(v); }

  interval_t at(const variable_t &v) const override {
    if (is_bottom())
      return interval_t::bottom();
    auto opt = get_var(v);
    if (!opt)
      return interval_t::top();
    const Tvpi::Interval<true> &ti = m_tvpi.getInterval(*opt);
    return from_tvpi_iv(ti);
  }

  // -------------------------------------------------------------------
  // Assume constraints
  // -------------------------------------------------------------------

  void operator+=(const linear_constraint_system_t &csts) override {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");

    if (is_bottom())
      return;

    for (auto const &cst : csts) {
      if (cst.is_contradiction()) {
        set_to_bottom();
        return;
      }
      if (cst.is_tautology())
        continue;
      add_constraint(cst);
      if (is_bottom())
        return;
    }
  }

  bool entails(const linear_constraint_t &cst) const override {
    if (is_bottom())
      return true;
    if (cst.is_tautology())
      return true;
    if (cst.is_contradiction())
      return false;
    // Entailment: NOT(cst) leads to bottom?
    tvpi_domain_t tmp(*this);
    tmp += cst.negate();
    return tmp.is_bottom();
  }

  // -------------------------------------------------------------------
  // Assignments
  // -------------------------------------------------------------------

  void assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");
    if (is_bottom())
      return;
    do_assign(x, e);
  }

  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.weak_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".weak_assign");
    if (is_bottom())
      return;
    tvpi_domain_t other(*this);
    other.do_assign(x, e);
    *this = *this | other;
  }

  // -------------------------------------------------------------------
  // Arithmetic apply
  // -------------------------------------------------------------------

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");
    if (is_bottom())
      return;

    switch (op) {
    case OP_ADDITION:
      do_assign(x, linear_expression_t(y) + k);
      break;
    case OP_SUBTRACTION:
      do_assign(x, linear_expression_t(y) - k);
      break;
    case OP_MULTIPLICATION:
      do_assign(x, k * linear_expression_t(y));
      break;
    default:
      // Division / remainder: conservative
      *this -= x;
      break;
    }
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");
    if (is_bottom())
      return;

    switch (op) {
    case OP_ADDITION:
      do_assign(x, linear_expression_t(y) + linear_expression_t(z));
      break;
    case OP_SUBTRACTION:
      do_assign(x, linear_expression_t(y) - linear_expression_t(z));
      break;
    case OP_MULTIPLICATION: {
      // If one operand is a singleton, reduce to scalar multiplication
      interval_t yi = at(y), zi = at(z);
      if (auto ky = yi.singleton()) {
        apply(op, x, z, *ky);
        return;
      }
      if (auto kz = zi.singleton()) {
        apply(op, x, y, *kz);
        return;
      }
      *this -= x;
      apply_interval(x, yi * zi);
      break;
    }
    default:
      *this -= x;
      break;
    }
  }

  // -------------------------------------------------------------------
  // Integer conversion / bitwise
  // -------------------------------------------------------------------

  void apply(int_conv_operation_t, const variable_t &dst,
             const variable_t &src) override {
    if (is_bottom())
      return;
    // Treat as copy (ignore bit-width differences)
    do_assign(dst, linear_expression_t(src));
  }

  void apply(bitwise_operation_t, const variable_t &x, const variable_t &,
             const variable_t &) override {
    if (!is_bottom())
      *this -= x;
  }

  void apply(bitwise_operation_t, const variable_t &x, const variable_t &,
             number_t) override {
    if (!is_bottom())
      *this -= x;
  }

  // -------------------------------------------------------------------
  // Backward operations (not implemented)
  // -------------------------------------------------------------------

  void backward_assign(const variable_t &, const linear_expression_t &,
                       const tvpi_domain_t &) override {
    CRAB_WARN(domain_name(), "::backward_assign not implemented");
  }

  void backward_apply(arith_operation_t, const variable_t &, const variable_t &,
                      number_t, const tvpi_domain_t &) override {
    CRAB_WARN(domain_name(), "::backward_apply not implemented");
  }

  void backward_apply(arith_operation_t, const variable_t &, const variable_t &,
                      const variable_t &, const tvpi_domain_t &) override {
    CRAB_WARN(domain_name(), "::backward_apply not implemented");
  }

  DEFAULT_SELECT(tvpi_domain_t)
  BOOL_OPERATIONS_NOT_IMPLEMENTED(tvpi_domain_t)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(tvpi_domain_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(tvpi_domain_t)

  // -------------------------------------------------------------------
  // Intrinsics
  // -------------------------------------------------------------------

  void intrinsic(std::string name,
                 const variable_or_constant_vector_t &,
                 const variable_vector_t &) override {
    CRAB_WARN(domain_name(), "::intrinsic for ", name, " not implemented");
  }

  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &,
                          const variable_vector_t &,
                          const tvpi_domain_t &) override {
    CRAB_WARN(domain_name(), "::backward_intrinsic for ", name,
              " not implemented");
  }

  // -------------------------------------------------------------------
  // Variable management: forget, project, rename, expand
  // -------------------------------------------------------------------

  void forget(const variable_vector_t &variables) override {
    if (is_bottom() || is_top())
      return;
    for (auto const &v : variables)
      *this -= v;
  }

  void project(const variable_vector_t &variables) override {
    crab::CrabStats::count(domain_name() + ".count.project");
    crab::ScopedCrabStats __st__(domain_name() + ".project");
    if (is_bottom() || is_top())
      return;

    // Identify which TvpiVar indices to keep
    std::vector<bool> keep(m_tvpi.size(), false);
    for (auto &v : variables) {
      auto opt = get_var(v);
      if (opt && *opt < (TvpiVar)keep.size())
        keep[*opt] = true;
    }

    // Build sorted projection sequence
    std::vector<TvpiVar> seq;
    for (TvpiVar i = 0; i < (TvpiVar)m_tvpi.size(); i++) {
      if (keep[i])
        seq.push_back(i);
    }

    if (seq.size() == m_tvpi.size())
      return; // Nothing to remove

    if (seq.empty()) {
      set_to_top();
      return;
    }

    m_tvpi.projectOnto(seq.size(), seq.data());

    // Rebuild maps
    var_map_t new_map;
    rev_map_t new_rev(seq.size(), boost::none);
    for (size_t i = 0; i < seq.size(); i++) {
      assert(m_rev_map[seq[i]]);
      const variable_t &v = *m_rev_map[seq[i]];
      new_map[v] = static_cast<TvpiVar>(i);
      new_rev[i] = v;
    }
    m_var_map = std::move(new_map);
    m_rev_map = std::move(new_rev);
  }

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    crab::CrabStats::count(domain_name() + ".count.rename");
    crab::ScopedCrabStats __st__(domain_name() + ".rename");
    if (is_bottom() || is_top())
      return;
    assert(from.size() == to.size());
    for (size_t i = 0; i < from.size(); i++) {
      auto it = m_var_map.find(from[i]);
      if (it == m_var_map.end())
        continue;
      TvpiVar idx = it->second;
      m_var_map.erase(it);
      m_var_map[to[i]] = idx;
      if (idx < (TvpiVar)m_rev_map.size())
        m_rev_map[idx] = to[i];
    }
  }

  void expand(const variable_t &x, const variable_t &new_x) override {
    crab::CrabStats::count(domain_name() + ".count.expand");
    crab::ScopedCrabStats __st__(domain_name() + ".expand");
    if (is_bottom() || is_top())
      return;
    // Conservative: set new_x to the interval of x
    interval_t xi = at(x);
    TvpiVar new_idx = get_or_create_var(new_x);
    Tvpi::Interval<true> ti = to_tvpi_iv(xi);
    m_tvpi.intersectBound(new_idx, ti);
  }

  void normalize() override {}
  void minimize() override {}

  // -------------------------------------------------------------------
  // Output
  // -------------------------------------------------------------------

  void write(crab_os &o) const override {
    if (is_bottom()) {
      o << "_|_";
      return;
    }
    if (is_top()) {
      o << "top";
      return;
    }
    linear_constraint_system_t csts = to_linear_constraint_system();
    o << csts;
  }

  friend crab_os &operator<<(crab_os &o, const tvpi_domain_t &val) {
    val.write(o);
    return o;
  }

  std::string domain_name() const override { return "TvpiDomain"; }

  // -------------------------------------------------------------------
  // Extract linear constraints
  // -------------------------------------------------------------------

  linear_constraint_system_t to_linear_constraint_system() const override {
    crab::CrabStats::count(domain_name() +
                           ".count.to_linear_constraint_system");
    crab::ScopedCrabStats __st__(domain_name() +
                                 ".to_linear_constraint_system");

    linear_constraint_system_t csts;
    if (is_bottom()) {
      csts += linear_constraint_t::get_false();
      return csts;
    }
    if (is_top())
      return csts;

    // 1. Interval bounds per variable
    for (auto &kv : m_var_map) {
      const variable_t &v = kv.first;
      TvpiVar idx = kv.second;
      interval_t ci = from_tvpi_iv(m_tvpi.getInterval(idx));
      if (auto lb = ci.lb().number())
        csts += linear_constraint_t(v >= *lb);
      if (auto ub = ci.ub().number())
        csts += linear_constraint_t(v <= *ub);
    }

    // 2. Relational constraints between pairs (x, y)
    TvpiVar n = static_cast<TvpiVar>(m_tvpi.size());
    for (TvpiVar yi = 0; yi < n; yi++) {
      for (TvpiVar xi = 0; xi < yi; xi++) {
        if (xi >= (TvpiVar)m_rev_map.size() ||
            yi >= (TvpiVar)m_rev_map.size())
          continue;
        if (!m_rev_map[xi] || !m_rev_map[yi])
          continue;

        const variable_t &vx = *m_rev_map[xi];
        const variable_t &vy = *m_rev_map[yi];
        if (m_var_map.find(vx) == m_var_map.end()) continue;
        if (m_var_map.find(vy) == m_var_map.end()) continue;

        const Tvpi::Polyhedron<true> &poly = m_tvpi.getProjection(xi, yi);
        size_t nineqs = poly.getNoOfInequalities();
        for (size_t k = 0; k < nineqs; k++) {
          const Tvpi::Inequality *ineq = poly[k];
          if (!ineq)
            continue;
          number_t a = from_mpz(ineq->getA());
          number_t b = from_mpz(ineq->getB());
          number_t c = from_mpz(ineq->getC());
          // a*vx + b*vy <= c
          linear_expression_t lhs =
              a * linear_expression_t(vx) + b * linear_expression_t(vy);
          csts += linear_constraint_t(lhs <= c);
        }
      }
    }
    return csts;
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    auto lin_csts = to_linear_constraint_system();
    if (lin_csts.is_false())
      return disjunctive_linear_constraint_system_t(true /*is_false*/);
    if (lin_csts.is_true())
      return disjunctive_linear_constraint_system_t(false /*is_false*/);
    return disjunctive_linear_constraint_system_t(lin_csts);
  }

private:
  // -------------------------------------------------------------------
  // Core assign implementation
  // -------------------------------------------------------------------

  void do_assign(const variable_t &x, const linear_expression_t &e) {
    // Case 1: e is a constant
    if (e.is_constant()) {
      number_t val = e.constant();
      *this -= x;
      TvpiVar idx = get_or_create_var(x);
      mpz_class v = to_mpz(val);
      Tvpi::Interval<true> i(v);
      m_tvpi.intersectBound(idx, i);
      return;
    }

    // Case 2: e = variable (possibly with coefficient 1)
    if (boost::optional<variable_t> yv = e.get_variable()) {
      assign_affine_one_var(x, number_t(1), *yv, number_t(0));
      return;
    }

    // Count non-zero terms
    int nterms = 0;
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      if (it->first != number_t(0))
        nterms++;
    }

    // Case 3: single variable with coefficient: a*y + c
    if (nterms == 1) {
      auto it = e.begin();
      while (it != e.end() && it->first == number_t(0))
        ++it;
      if (it == e.end()) {
        // Constant (already handled above, but just in case)
        do_assign(x, linear_expression_t(e.constant()));
        return;
      }
      assign_affine_one_var(x, it->first, it->second, e.constant());
      return;
    }

    // Case 4: two variables: a*y + b*z + c
    if (nterms == 2) {
      // Collect the two terms
      std::vector<std::pair<number_t, variable_t>> terms;
      for (auto it = e.begin(), et = e.end(); it != et; ++it) {
        if (it->first != number_t(0))
          terms.push_back({it->first, it->second});
      }
      assert(terms.size() == 2);

      number_t a = terms[0].first, b = terms[1].first;
      const variable_t &y = terms[0].second;
      const variable_t &z = terms[1].second;
      number_t c = e.constant();

      // Ensure y and z are tracked before forgetting x
      get_or_create_var(y);
      get_or_create_var(z);

      // Forget x: indices of y and z may shift if they were after x
      *this -= x;

      // Re-fetch y and z indices after the deletion
      auto opt_y = get_var(y);
      auto opt_z = get_var(z);
      assert(opt_y && opt_z);

      TvpiVar x_idx = get_or_create_var(x);

      // DenseTvpi::approximateInequality uses LinComponent::variable as the
      // direct TvpiVar index.
      // Add: x - a*y - b*z <= c   and   -x + a*y + b*z <= -c
      // Together these encode  x == a*y + b*z + c
      auto add2 = [&](number_t cx, number_t cy, number_t cz,
                      number_t rhs) -> bool {
        std::vector<LinComponent> comps;
        if (cx != number_t(0))
          comps.push_back(LinComponent(
              to_mpz(cx), static_cast<Variable>(x_idx)));
        if (cy != number_t(0))
          comps.push_back(LinComponent(
              to_mpz(cy), static_cast<Variable>(*opt_y)));
        if (cz != number_t(0))
          comps.push_back(LinComponent(
              to_mpz(cz), static_cast<Variable>(*opt_z)));
        if (comps.empty())
          return true;
        std::sort(comps.begin(), comps.end());
        Result res =
            m_tvpi.approximateInequality(comps, to_mpz(rhs), false);
        return (res != resUnsatisfiable);
      };

      if (!add2(number_t(1), -a, -b, c)) {
        m_is_bottom = true;
        return;
      }
      if (!add2(number_t(-1), a, b, -c))
        m_is_bottom = true;
      return;
    }

    // General case (>2 terms): forget x, use interval arithmetic
    *this -= x;
    interval_t xi = eval_expr_interval(e);
    if (!xi.is_top())
      apply_interval(x, xi);
    else
      get_or_create_var(x); // just add x unbounded
  }
};

// ---------------------------------------------------------------------------
// abstract_domain_traits
// ---------------------------------------------------------------------------

template <typename Number, typename VariableName>
struct abstract_domain_traits<tvpi_domain<Number, VariableName>> {
  using number_t = Number;
  using varname_t = VariableName;
};

} // namespace domains
} // namespace crab
