#pragma once

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>
#include <string>

namespace crab {
namespace domains {

/**
 ** This domain when instantiated with the Octagon domain simulates
 ** the Two Variables Per Inequalities (tvpi)
 ** (http://www2.in.tum.de/bib/files/simon02two.pdf) domain by Simon,
 ** King, and Howe but it fixes a priori the values of the non-unit
 ** coefficients to a small predefined set which can be chosen by
 ** users.
 **
 ** Algorithmically, this implementation has nothing to do with
 ** Simon's. This is much simpler (but less expressive) based on
 ** pretty simple rewrites. For each variable v we add a ghost
 ** variable GHOST(v,N) where N is a fixed coefficient and the ghost
 ** variable represents v/N. We then rewrite operands in assignments
 ** and arithmetic operations, and linear constraints to use
 ** GHOST(v,N).
 **
 **
 ** For instance, robust_tvpi_domain cannot prove this program while tvpi does:
 **    int i,x;
 **    int N = nd_int();
 **   __CRAB_assume(N > 0);
 **   i = 0;
 **   x = 0;
 **   while (i < N) {
 **     i++;
 **     if (*) {
 **      x = x+2;
 **     } else {
 **      x = x+3;
 **     }
 **   }
 **   __CRAB_assert(x >= 2*N);
 **   __CRAB_assert(x <= 3*N);
 **
 ** The reason is that we assign two different ghost variables for x/2
 ** and x/3 so at the join point of the if-then-else we lose track of
 ** them.
 **/
template <typename OctLikeDomain>
class robust_tvpi_domain
    : public abstract_domain_api<robust_tvpi_domain<OctLikeDomain>> {
public:
  using robust_tvpi_domain_t = robust_tvpi_domain<OctLikeDomain>;
  using abstract_domain_api_t = abstract_domain_api<robust_tvpi_domain_t>;
  using typename abstract_domain_api_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_api_t::interval_t;
  using typename abstract_domain_api_t::linear_constraint_system_t;
  using typename abstract_domain_api_t::linear_constraint_t;
  using typename abstract_domain_api_t::linear_expression_t;
  using typename abstract_domain_api_t::number_t;
  using typename abstract_domain_api_t::reference_constraint_t;
  using typename abstract_domain_api_t::variable_or_constant_t;
  using typename abstract_domain_api_t::variable_or_constant_vector_t;
  using typename abstract_domain_api_t::variable_t;
  using typename abstract_domain_api_t::variable_vector_t;
  using typename abstract_domain_api_t::varname_t;

private:
  using base_domain_t = OctLikeDomain;
  using coeff_range_t = typename base_domain_t::interval_t;
  using var_map_t = ikos::separate_domain<variable_t, coeff_range_t>;

  base_domain_t m_base_absval;
  var_map_t m_var_coeff_map;

  static variable_t get_ghost_var(const variable_t &v, number_t coefficient) {
    if (coefficient <= number_t(1)) {
      CRAB_ERROR("Coefficient must be greater than 1");
    }
    auto &vfac = const_cast<varname_t *>(&(v.name()))->get_var_factory();
    std::string ghost_name = ".tvpi.ghost_var(" + coefficient.get_str() + ")";
    return variable_t(vfac.get(v.name(), ghost_name), v.get_type());
  }

  // (TODO): needed for pretty-printing
  // static std::pair<variable_t, unsigned> get_rev_ghost_var(const variable_t
  // &gv) {
  //
  //}

  boost::optional<number_t>
  get_range_coefficient_lb(const variable_t &x) const {
    if (const coeff_range_t *ival = m_var_coeff_map.find(x)) {
      if ((*ival).is_bottom()) {
        return boost::none;
      } else {
        return (*ival).lb().number();
      }
    }
    return boost::none;
  }

  boost::optional<number_t>
  get_range_coefficient_ub(const variable_t &x) const {
    if (const coeff_range_t *ival = m_var_coeff_map.find(x)) {
      if ((*ival).is_bottom()) {
        return boost::none;
      } else {
        return (*ival).ub().number();
      }
    }
    return boost::none;
  }

  bool is_coefficient_in_bounds(boost::optional<number_t> lb_opt,
                                boost::optional<number_t> ub_opt) const {

    if (lb_opt == boost::none || ub_opt == boost::none) {
      return false;
    }

    return true;
  }

  void add_coeff_range(const variable_t &x, number_t coefficient) {
    if (coefficient <= number_t(1)) {
      return;
    }

    coeff_range_t coeff_range(coefficient);
    if (const coeff_range_t *ival = m_var_coeff_map.find(x)) {
      if (coeff_range <= *ival) {
        return;
      }
      coeff_range = coeff_range | *ival;
    }
    m_var_coeff_map.set(x, coeff_range);
  }

  void add_coeff_range(const variable_t &x, coeff_range_t coeff_range) {
    if (coeff_range.is_bottom()) {
      return;
    }

    if (const coeff_range_t *ival = m_var_coeff_map.find(x)) {
      if (coeff_range <= *ival) {
        return;
      }
      coeff_range = coeff_range | *ival;
    }
    m_var_coeff_map.set(x, coeff_range);
  }

  coeff_range_t
  compute_linear_expression_coeff_range(const linear_expression_t &e) const {
    coeff_range_t ret_ival = coeff_range_t::bottom();
    for (auto it = e.variables_begin(); it != e.variables_end(); it++) {
      if (const coeff_range_t *ival = m_var_coeff_map.find(*it)) {
        ret_ival = ret_ival | *ival;
      }
    }
    return ret_ival;
  }

  coeff_range_t
  compute_linear_constraint_coeff_range(const linear_constraint_t &cst) const {
    return compute_linear_expression_coeff_range(cst.expression());
  }

  coeff_range_t
  infer_linear_expression_coeff(const linear_expression_t &e) const {
    coeff_range_t ret_ival = coeff_range_t::bottom();
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const number_t &coeff = (*it).first;
      if (coeff > number_t(1)) {
        ret_ival = ret_ival | coeff_range_t(coeff);
      } else if (coeff < number_t(-1)) {
        ret_ival = ret_ival | coeff_range_t(-coeff);
      }
    }
    return ret_ival;
  }

  void assign_ghost_var(const variable_t &var, number_t coefficient) {
    if (coefficient <= number_t(1)) {
      return;
    }
    variable_t ghost_var = get_ghost_var(var, coefficient);
    // m_base_absval.assign(ghost_var, var * (number_t(1) / coefficient));
    m_base_absval += linear_constraint_t(var == ghost_var * coefficient);
  }

  void propagate_coeff(const variable_t &dst, const variable_t &src) {
    if (const coeff_range_t *ival = m_var_coeff_map.find(src)) {
      m_var_coeff_map.set(dst, *ival);
    }
  }

  void insert_coeff(const variable_t &var, number_t z) {
    if (z <= number_t(1)) {
      return;
    }

    coeff_range_t coeff_range(z);
    const coeff_range_t *ival = m_var_coeff_map.find(var);
    if (ival) { // found, check range
      if (coeff_range <= *ival) {
        return;
      } else {
        assign_ghost_var(var, z);
      }
      coeff_range = coeff_range | *ival;
    } else {
      assign_ghost_var(var, z);
    }
    m_var_coeff_map.set(var, coeff_range);
  }

  bool can_rewrite_linear_expression(const linear_expression_t &e,
                                     number_t coefficient) const {
    assert(coefficient > number_t(1));

    number_t tracked_coeff = coefficient;
    if (!((e.constant() % tracked_coeff) == 0)) {
      return false;
    }

    return true;

    // for (auto it = e.begin(), et = e.end(); it != et; ++it) {
    //   const number_t &coefficient = (*it).first;
    //   if (coefficient == tracked_coeff || coefficient == -tracked_coeff) {
    //     return true;
    //   }
    // }
    // return false;
  }

  // clang-format off
  /**
   *
   * Given c1*x1 + c2*x2 +... + k, rewrite each ci*xi into
   *
   * ui*xi         if ci is one of tracked coefficients then ui = +1 (-1) if ci > 0 (c1 < 0)
   * ci*ghost(xi)  otherwise
   **/
  // clang-format on

  linear_expression_t rewrite_linear_expression(const linear_expression_t &e,
                                                number_t coefficient) const {
    assert(can_rewrite_linear_expression(e, coefficient));

    number_t tracked_coeff = coefficient;
    // This is okay because can_rewrite_linear_expression guarantees
    // that e.constant() is divisible by coefficient.
    linear_expression_t res(e.constant() / tracked_coeff);
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const variable_t &v = (*it).second;
      const number_t &coeff = (*it).first;
      if (coeff == number_t(0)) {
        continue;
      } else if (coeff == tracked_coeff) {
        res = res + v;
      } else if (coeff == -tracked_coeff) {
        res = res - v;
      } else {
        variable_t ghost_v = get_ghost_var(v, coefficient);
        res = res + coeff * ghost_v;
      }
    }
    return res;
  }

  bool can_rewrite_linear_constraint(const linear_constraint_t &cst,
                                     number_t coefficient) const {
    return can_rewrite_linear_expression(cst.expression(), coefficient);
  }

  linear_constraint_t rewrite_linear_constraint(const linear_constraint_t &cst,
                                                number_t coefficient) const {
    if (cst.is_equality() || cst.is_disequation()) {
      return linear_constraint_t(
          rewrite_linear_expression(cst.expression(), coefficient), cst.kind());
    } else {
      return linear_constraint_t(
          rewrite_linear_expression(cst.expression(), coefficient), cst.kind(),
          cst.is_signed());
    }
  }

  void rewrite_assign(const variable_t &x, const linear_expression_t &e,
                      number_t coefficient) {
    assert(coefficient > number_t(1));

    insert_coeff(x, coefficient);

    variable_t ghost_x = get_ghost_var(x, coefficient);
    auto x_intv = m_base_absval.at(x);
    number_t tracked_coefficient = coefficient;
    if (boost::optional<number_t> n = x_intv.singleton()) {
      if ((*n) % tracked_coefficient == 0) {
        // rewrite("x := n") = "x/COEF := n/COEF"
        m_base_absval.assign(ghost_x, (*n) / tracked_coefficient);
      } else {
        m_base_absval -= ghost_x;
      }
    } else if (boost::optional<variable_t> y = e.get_variable()) {
      // rewrite("x := y") = "x/COEF := y/COEF"
      insert_coeff(*y, coefficient);
      variable_t ghost_y = get_ghost_var(*y, coefficient);
      m_base_absval.assign(ghost_x, ghost_y);
    } else {
      if (can_rewrite_linear_expression(e, coefficient)) {
        m_base_absval.assign(ghost_x,
                             rewrite_linear_expression(e, coefficient));
      } else {
        m_base_absval -= ghost_x;
      }
    }
  }

  void rewrite_apply(arith_operation_t op, const variable_t &x,
                     const variable_t &y, number_t z, number_t coefficient) {
    assert(coefficient > number_t(1));

    insert_coeff(y, coefficient);
    insert_coeff(x, coefficient);

    variable_t ghost_x = get_ghost_var(x, coefficient);
    number_t tracked_coefficient = coefficient;
    if (op == OP_ADDITION || op == OP_SUBTRACTION) {
      if ((z % tracked_coefficient) == 0) {
        // rewrite("x := y + COEF") = "x/COEF := y/COEF + 1"
        // rewrite("x := y - COEF") = "x/COEF := y/COEF - 1"
        variable_t ghost_y = get_ghost_var(y, coefficient);
        m_base_absval.apply(op, ghost_x, ghost_y, z / tracked_coefficient);
        return;
      }
    } else if (op == OP_MULTIPLICATION) {
      if (z == number_t(1)) {
        // rewrite("x := y") = "x/COEF := y/COEF"
        variable_t ghost_y = get_ghost_var(y, coefficient);
        m_base_absval.assign(ghost_x, ghost_y);
        return;
      } else if ((z % tracked_coefficient) == 0) {
        // rewrite("x := COEF * y") = "x/COEF := y"
        m_base_absval.assign(ghost_x, y);
        return;
      }
    } else if (op == OP_SDIV) {
      if (z == number_t(1)) {
        // rewrite("x := y") = "x/COEF := y/COEF"
        variable_t ghost_y = get_ghost_var(y, coefficient);
        m_base_absval.assign(ghost_x, ghost_y);
        return;
      } else if ((z % tracked_coefficient) == 0) {
        // rewrite("x := y/COEF") =  "x := y/COEF"
        variable_t ghost_y = get_ghost_var(y, coefficient);
        m_base_absval.assign(x, ghost_y);
        return;
      }
    }

    // default case: forget x
    m_base_absval -= ghost_x;
  }

  void rewrite_apply_var(arith_operation_t op, const variable_t &x,
                         const variable_t &y, const variable_t &z,
                         number_t coefficient) {
    assert(coefficient > number_t(1));
    // TODO: ignored case if op is a subtraction or division
    // TODO: ignored cases if y or z are not singleton

    insert_coeff(y, coefficient);
    insert_coeff(x, coefficient);
    insert_coeff(z, coefficient);

    if (boost::optional<number_t> n = m_base_absval.at(z).singleton()) {
      rewrite_apply(op, x, y, *n, coefficient);
      return;
    }

    if (op == OP_ADDITION || op == OP_MULTIPLICATION) {
      if (boost::optional<number_t> n = m_base_absval.at(y).singleton()) {
        rewrite_apply(op, x, z, *n, coefficient);
        return;
      }
    }

    // default case: forget x
    variable_t ghost_x = get_ghost_var(x, coefficient);
    m_base_absval -= ghost_x;
  }

  void insert_coeff_by_var(const variable_t &x, const variable_t &y) {

    if (boost::optional<number_t> n = m_base_absval.at(y).singleton()) {
      if (*n > number_t(1)) {
        insert_coeff(x, *n);
        insert_coeff(y, *n);
      }
    }
  }

  robust_tvpi_domain(base_domain_t &&num, var_map_t &&var_coeff_map)
      : m_base_absval(std::move(num)),
        m_var_coeff_map(std::move(var_coeff_map)) {}

public:
  robust_tvpi_domain() {}

  robust_tvpi_domain(const robust_tvpi_domain_t &o) = default;
  robust_tvpi_domain(robust_tvpi_domain_t &&o) = default;
  robust_tvpi_domain_t &operator=(const robust_tvpi_domain_t &o) {
    if (this != &o) {
      m_base_absval = o.m_base_absval;
      m_var_coeff_map = o.m_var_coeff_map;
    }
    return *this;
  }
  robust_tvpi_domain_t &operator=(robust_tvpi_domain_t &&o) {
    if (this != &o) {
      m_base_absval = std::move(o.m_base_absval);
      m_var_coeff_map = std::move(o.m_var_coeff_map);
    }
    return *this;
  }

  void set_to_top() override { m_base_absval.set_to_top(); }

  void set_to_bottom() override { m_base_absval.set_to_bottom(); }

  robust_tvpi_domain_t make_bottom() const override {
    robust_tvpi_domain_t res;
    res.set_to_bottom();
    return res;
  }

  robust_tvpi_domain_t make_top() const override {
    robust_tvpi_domain_t res;
    return res;
  }

  bool is_bottom() const override { return m_base_absval.is_bottom(); }

  bool is_top() const override { return m_base_absval.is_top(); }

  bool operator<=(const robust_tvpi_domain_t &other) const override {
    if (is_bottom() || other.is_top()) {
      return true;
    } else if (is_top() || other.is_bottom()) {
      return false;
    } else {
      bool res = m_base_absval <= other.m_base_absval;
      res &= (m_var_coeff_map <= other.m_var_coeff_map);
      return res;
    }
  }

  void operator|=(const robust_tvpi_domain_t &other) override {
    if (is_bottom() || other.is_top()) {
      *this = other;
    } else if (other.is_bottom() || is_top()) {
      // do nothing
    } else {
      m_base_absval |= other.m_base_absval;
      m_var_coeff_map = m_var_coeff_map | other.m_var_coeff_map;
    }
  }

  robust_tvpi_domain_t
  operator|(const robust_tvpi_domain_t &other) const override {
    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    } else {
      base_domain_t out_base_absval = m_base_absval | other.m_base_absval;
      var_map_t out_var_coeff_map = m_var_coeff_map | other.m_var_coeff_map;
      return robust_tvpi_domain_t(std::move(out_base_absval),
                                  std::move(out_var_coeff_map));
    }
  }

  void operator&=(const robust_tvpi_domain_t &other) override {
    if (is_bottom() || other.is_top()) {
      // do nothing
    } else if (other.is_bottom() || is_top()) {
      *this = other;
    } else {
      m_base_absval &= other.m_base_absval;
      m_var_coeff_map = m_var_coeff_map & other.m_var_coeff_map;
    }
  }

  robust_tvpi_domain_t
  operator&(const robust_tvpi_domain_t &other) const override {
    if (is_bottom() || other.is_top()) {
      return *this;
    } else if (other.is_bottom() || is_top()) {
      return other;
    } else {
      base_domain_t out_base_absval = m_base_absval & other.m_base_absval;
      var_map_t out_var_coeff_map = m_var_coeff_map & other.m_var_coeff_map;
      return robust_tvpi_domain_t(std::move(out_base_absval),
                                  std::move(out_var_coeff_map));
    }
  }

  robust_tvpi_domain_t
  operator||(const robust_tvpi_domain_t &other) const override {
    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    } else {
      base_domain_t out_base_absval = m_base_absval || other.m_base_absval;
      var_map_t out_var_coeff_map = m_var_coeff_map || other.m_var_coeff_map;
      return robust_tvpi_domain_t(std::move(out_base_absval),
                                  std::move(out_var_coeff_map));
    }
  }

  robust_tvpi_domain_t
  widening_thresholds(const robust_tvpi_domain_t &other,
                      const thresholds<number_t> &ts) const override {
    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    } else {
      base_domain_t out_base_absval =
          m_base_absval.widening_thresholds(other.m_base_absval, ts);
      var_map_t out_var_coeff_map =
          m_var_coeff_map.widening_thresholds(other.m_var_coeff_map, ts);
      return robust_tvpi_domain_t(std::move(out_base_absval),
                                  std::move(out_var_coeff_map));
    }
  }

  robust_tvpi_domain_t
  operator&&(const robust_tvpi_domain_t &other) const override {
    if (is_bottom() || other.is_top()) {
      return *this;
    } else if (other.is_bottom() || is_top()) {
      return other;
    } else {
      base_domain_t out_base_absval = m_base_absval && other.m_base_absval;
      var_map_t out_var_coeff_map = m_var_coeff_map && other.m_var_coeff_map;
      return robust_tvpi_domain_t(std::move(out_base_absval),
                                  std::move(out_var_coeff_map));
    }
  }

  void operator-=(const variable_t &var) override {
    if (!(is_bottom() || is_top())) {
      m_base_absval -= var;

      boost::optional<number_t> lb_opt = get_range_coefficient_lb(var);
      boost::optional<number_t> ub_opt = get_range_coefficient_ub(var);

      if (is_coefficient_in_bounds(lb_opt, ub_opt)) {
        for (number_t coefficient = *lb_opt; coefficient <= *ub_opt;
             coefficient++) {
          variable_t ghost_var = get_ghost_var(var, coefficient);
          m_base_absval -= ghost_var;
        }
      }
    }
  }

  interval_t operator[](const variable_t &v) override {
    if (is_bottom()) {
      return interval_t::bottom();
    }
    // REVISIT: not sure we need to do some rewriting here
    return m_base_absval[v];
  }

  interval_t at(const variable_t &v) const override {
    if (is_bottom()) {
      return interval_t::bottom();
    }
    // REVISIT: not sure we need to do some rewriting here
    return m_base_absval.at(v);
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    CRAB_LOG("fixed-tvpi",
             crab::outs() << "Before assume(" << csts << ")=" << *this << "\n");
    if (!is_bottom()) {
      for (auto const &cst : csts) {
        if (cst.is_contradiction()) {
          set_to_bottom();
          break;
        }
        if (cst.is_tautology()) {
          continue;
        }

        CRAB_LOG("fixed-tvpi", crab::outs() << "robust_tvpi_domain processing "
                                            << cst << "\n");
        m_base_absval += cst;
        if (m_base_absval.is_bottom()) {
          break;
        }

        interval_t ival = compute_linear_constraint_coeff_range(cst);

        boost::optional<number_t> lb_opt = ival.lb().number();
        boost::optional<number_t> ub_opt = ival.ub().number();

        if (is_coefficient_in_bounds(lb_opt, ub_opt)) {
          for (number_t coefficient = *lb_opt; coefficient <= *ub_opt;
               coefficient++) {
            if (can_rewrite_linear_constraint(cst, coefficient)) {
              auto gcst = rewrite_linear_constraint(cst, coefficient);
              CRAB_LOG("fixed-tvpi", crab::outs() << "\tRewritten " << cst
                                                  << " into " << gcst << "\n");
              m_base_absval += gcst;
              if (m_base_absval.is_bottom()) {
                break;
              }
            } else {
              CRAB_LOG("fixed-tvpi", crab::outs()
                                         << "\tCannot rewrite " << cst << "\n");
            }
          }
        }
      } // end for
    }

    CRAB_LOG("fixed-tvpi",
             crab::outs() << "After assume(" << csts << ")=" << *this << "\n");
  }

  bool entails(const linear_constraint_t &cst) const override {
    if (is_bottom()) {
      return true;
    } else if (cst.is_tautology()) {
      return true;
    } else if (cst.is_contradiction()) {
      return false;
    }

    if (m_base_absval.entails(cst)) {
      return true;
    }

    interval_t ival = compute_linear_constraint_coeff_range(cst);

    boost::optional<number_t> lb_opt = ival.lb().number();
    boost::optional<number_t> ub_opt = ival.ub().number();

    if (is_coefficient_in_bounds(lb_opt, ub_opt)) {
      for (number_t coefficient = *lb_opt; coefficient <= *ub_opt;
           coefficient++) {
        if (can_rewrite_linear_constraint(cst, coefficient)) {
          linear_constraint_t gcst =
              rewrite_linear_constraint(cst, coefficient);
          if (m_base_absval.entails(gcst)) {
            return true;
          }
        }
      }
    }

    return false;
  }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    if (!is_bottom()) {
      m_base_absval.assign(x, e);

      coeff_range_t ival = infer_linear_expression_coeff(e);
      add_coeff_range(x, ival);
      if (boost::optional<variable_t> y = e.get_variable()) {
        propagate_coeff(x, *y);
      }

      boost::optional<number_t> lb_opt = get_range_coefficient_lb(x);
      boost::optional<number_t> ub_opt = get_range_coefficient_ub(x);

      if (is_coefficient_in_bounds(lb_opt, ub_opt)) {
        for (number_t coefficient = *lb_opt; coefficient <= *ub_opt;
             coefficient++) {
          rewrite_assign(x, e, coefficient);
        }
      }
    }
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t z) override {
    if (!is_bottom()) {
      m_base_absval.apply(op, x, y, z);

      insert_coeff(y, z);
      insert_coeff(x, z);

      boost::optional<number_t> lb_opt = get_range_coefficient_lb(y);
      boost::optional<number_t> ub_opt = get_range_coefficient_ub(y);

      if (is_coefficient_in_bounds(lb_opt, ub_opt)) {
        for (number_t coefficient = *lb_opt; coefficient <= *ub_opt;
             coefficient++) {
          rewrite_apply(op, x, y, z, coefficient);
        }
      }
    }
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    if (!is_bottom()) {
      m_base_absval.apply(op, x, y, z);

      insert_coeff_by_var(x, y);
      insert_coeff_by_var(x, z);

      boost::optional<number_t> lb_opt = get_range_coefficient_lb(y);
      boost::optional<number_t> ub_opt = get_range_coefficient_ub(y);

      if (is_coefficient_in_bounds(lb_opt, ub_opt)) {
        for (number_t coefficient = *lb_opt; coefficient <= *ub_opt;
             coefficient++) {
          rewrite_apply_var(op, x, y, z, coefficient);
        }
      }

      lb_opt = get_range_coefficient_lb(z);
      ub_opt = get_range_coefficient_ub(z);

      if (is_coefficient_in_bounds(lb_opt, ub_opt)) {
        for (number_t coefficient = *lb_opt; coefficient <= *ub_opt;
             coefficient++) {
          rewrite_apply_var(op, x, y, z, coefficient);
        }
      }
    }
  }

  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    if (!is_bottom()) {
      m_base_absval.apply(op, dst, src);

      boost::optional<number_t> lb_opt = get_range_coefficient_lb(src);
      boost::optional<number_t> ub_opt = get_range_coefficient_ub(src);

      if (is_coefficient_in_bounds(lb_opt, ub_opt)) {
        for (number_t coefficient = *lb_opt; coefficient <= *ub_opt;
             coefficient++) {
          // WARNING: we ignore the actual conversion because we assume
          // mathematical integers.
          variable_t ghost_dst = get_ghost_var(dst, coefficient);
          variable_t ghost_src = get_ghost_var(src, coefficient);
          m_base_absval.assign(ghost_dst, ghost_src);
        }
      }
    }
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    if (!is_bottom()) {
      m_base_absval.apply(op, x, y, z);

      boost::optional<number_t> lb_opt = get_range_coefficient_lb(x);
      boost::optional<number_t> ub_opt = get_range_coefficient_ub(x);

      if (is_coefficient_in_bounds(lb_opt, ub_opt)) {
        for (number_t coefficient = *lb_opt; coefficient <= *ub_opt;
             coefficient++) {
          variable_t ghost_x = get_ghost_var(x, coefficient);
          m_base_absval -= ghost_x;
        }
      }
    }
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    if (!is_bottom()) {
      m_base_absval.apply(op, x, y, k);

      boost::optional<number_t> lb_opt = get_range_coefficient_lb(x);
      boost::optional<number_t> ub_opt = get_range_coefficient_ub(x);

      if (is_coefficient_in_bounds(lb_opt, ub_opt)) {
        for (number_t coefficient = *lb_opt; coefficient <= *ub_opt;
             coefficient++) {
          variable_t ghost_x = get_ghost_var(x, coefficient);
          m_base_absval -= ghost_x;
        }
      }
    }
  }

  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const robust_tvpi_domain_t &invariant) override {
    CRAB_WARN(domain_name(), "::backward_assign not implemented");
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const robust_tvpi_domain_t &invariant) override {
    CRAB_WARN(domain_name(), "::backward_apply not implemented");
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const robust_tvpi_domain_t &invariant) override {
    CRAB_WARN(domain_name(), "::backward_apply not implemented");
  }

  DEFAULT_SELECT(robust_tvpi_domain_t)
  BOOL_OPERATIONS_NOT_IMPLEMENTED(robust_tvpi_domain_t)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(robust_tvpi_domain_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(robust_tvpi_domain_t)

  linear_constraint_system_t to_linear_constraint_system() const override {

    if (is_bottom()) {
      return linear_constraint_system_t(linear_constraint_t::get_false());
    }

    if (is_top()) {
      return linear_constraint_system_t(linear_constraint_t::get_true());
    }

    // TODO: eliminate ghost variables
    return m_base_absval.to_linear_constraint_system();
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    CRAB_WARN(domain_name(),
              "::to_disjunctive_linear_constraint_system not implemented");
    disjunctive_linear_constraint_system_t res;
    return res;
  }

  void forget(const variable_vector_t &variables) override {
    if (!(is_bottom() || is_top())) {
      variable_vector_t allvars(variables);

      for (auto const &v : variables) {
        boost::optional<number_t> lb_opt = get_range_coefficient_lb(v);
        boost::optional<number_t> ub_opt = get_range_coefficient_ub(v);

        if (is_coefficient_in_bounds(lb_opt, ub_opt)) {
          for (number_t coefficient = *lb_opt; coefficient <= *ub_opt;
               coefficient++) {
            variable_t gv = get_ghost_var(v, coefficient);
            allvars.push_back(gv);
          }
        }
      }

      m_base_absval.forget(allvars);
    }
  }

  void project(const variable_vector_t &variables) override {
    if (!is_bottom()) {
      variable_vector_t allvars(variables);

      for (auto const &v : variables) {
        boost::optional<number_t> lb_opt = get_range_coefficient_lb(v);
        boost::optional<number_t> ub_opt = get_range_coefficient_ub(v);

        if (is_coefficient_in_bounds(lb_opt, ub_opt)) {
          for (number_t coefficient = *lb_opt; coefficient <= *ub_opt;
               coefficient++) {
            variable_t gv = get_ghost_var(v, coefficient);
            allvars.push_back(gv);
          }
        }
      }

      m_base_absval.project(allvars);
    }
  }

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    CRAB_WARN(domain_name(), "::rename not implemented");
  }

  void expand(const variable_t &var, const variable_t &new_var) override {
    if (is_bottom() || is_top()) {
      return;
    }

    m_base_absval.expand(var, new_var);

    boost::optional<number_t> lb_opt = get_range_coefficient_lb(var);
    boost::optional<number_t> ub_opt = get_range_coefficient_ub(var);

    if (is_coefficient_in_bounds(lb_opt, ub_opt)) {
      for (number_t coefficient = *lb_opt; coefficient <= *ub_opt;
           coefficient++) {
        variable_t gv = get_ghost_var(var, coefficient);
        variable_t gnv = get_ghost_var(new_var, coefficient);
        m_base_absval.expand(gv, gnv);
      }
    }
  }

  void normalize() override {}
  void minimize() override {}

  void intrinsic(std::string name, const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    CRAB_WARN(domain_name(), "::instrinsic for ", name, " not implemented");
  }

  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const robust_tvpi_domain_t &invariant) override {
    CRAB_WARN(domain_name(), "::backward_intrinsic for ", name,
              " not implemented");
  }

  void write(crab_os &o) const override {
    if (is_bottom()) {
      o << "_|_\n";
      m_var_coeff_map.write(o);
    } else if (is_top()) {
      o << "top\n";
      m_var_coeff_map.write(o);
    } else {
      // TODO: eliminate ghost variables
      m_base_absval.write(o);
      o << "\n";
      m_var_coeff_map.write(o);
    }
  }

  friend crab_os &operator<<(crab_os &o, const robust_tvpi_domain_t &val) {
    val.write(o);
    return o;
  }

  std::string domain_name() const override {
    base_domain_t dom;
    return "RobustTVPIDomain(" + dom.domain_name() + ")";
  }
};

template <typename OctLikeDomain>
struct abstract_domain_traits<robust_tvpi_domain<OctLikeDomain>> {
  using number_t = typename OctLikeDomain::number_t;
  using varname_t = typename OctLikeDomain::varname_t;
};

} // end namespace domains
} // end namespace crab
