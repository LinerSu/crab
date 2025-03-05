#pragma once

#include <algorithm>
#include <boost/optional.hpp>
#include <chrono>
#include <string>

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_params.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/inter_abstract_operations.hpp>
#include <crab/domains/tvpi/coefficient_map.hpp>
#include <crab/numbers/bignums.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>

namespace crab {
namespace domains {

#define tvpi_dbm_domain_SCOPED_STATS(NAME)                                     \
  CRAB_DOMAIN_SCOPED_STATS(this, NAME, 1)
#define tvpi_dbm_domain_SCOPED_STATS_ASSIGN_CTOR(NAME)                         \
  CRAB_DOMAIN_SCOPED_STATS(&o, NAME, 0)

class TVPIDBMDefaultParams {
public:
  enum { implement_inter_transformers = 0 };
};

template <typename OctLikeDomain, typename Params = TVPIDBMDefaultParams>
class tvpi_dbm_domain
    : public abstract_domain_api<tvpi_dbm_domain<OctLikeDomain, Params>> {
public:
  using tvpi_dbm_domain_t = tvpi_dbm_domain<OctLikeDomain, Params>;
  using abstract_domain_api_t = abstract_domain_api<tvpi_dbm_domain_t>;
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
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
  using coefficient_map_t =
      typename tvpi_utils::coefficient_map<variable_t,
                                           unsigned>::coefficient_map_t;
  using coefficient_set_t = typename coefficient_map_t::coefficient_set_t;
#else
  using coefficient_set_t = std::vector<unsigned>;
#endif

  static_assert(std::is_same<typename abstract_domain_api_t::number_t,
                             ikos::z_number>::value,
                "abstract_domain_api_t::number_t must be the ikos::z_number");

  // This domain is a special handling of TVPI constraints based on Difference
  // Bound Matrices (DBMs).
  // The TVPI constraint we are interested in is of the form:
  //          ax - by <= c | +/-x <= c
  // where a, b are positive integers and c is a number.
  // The domain keeps two DBMs:
  // - Classical DBM (m_base): Supports constraints of the form ±x <= c and x -
  // y <= c.
  // - Extended DBM (m_ext): Supports constraints of the form ax - by <= c.

  // For coefficients which is greater than 1, we introduce ghost variables
  //  ax as a * x.
  // The domain will keep a map from each variable to a set of coefficients.
  // This set of coefficients will be imprecise interms of join and widening.
  // The reason we ignore merge coefficients is that
  // loop may introduce new coefficient but later those may not be needed.

  // The way to infer coefficients is based on the following pattern:
  //  y := x * c | x / c
  //  Or semantically knows:
  //  y := x * z where z is a constant value

  // The domain will only keep a TVPI constraint:
  //          ax - by <= c | ax <= c
  // in a representative form:
  //          a/dx - b/dy <= c/d | x <= c/a
  //   where d = gcd(a, b).
  // For integer constant c, we approximate by \floor(c/d) and \floor(c/a).
  // Diophantine? In code, we called normalize.

  // The majority reduction is based on linear arithmetic:
  // Combining two inequalities and remove one intermidiate variable.
  // In general, it follows Fourier-Motzkin elimination method and is similar to
  // TVPI resultant operation. Since DBM has closure operation, we only need to
  // perform resultant for inequalities between m_base and m_ext.
  // Basically,
  //  ax - by <= c && y - z <= f => a'x - b'z <= c'
  // The coefficients and constant are normalized.
  // Besides, we skip to keep inequalities if the coefficients are not tracked.

private:
  using base_domain_t = OctLikeDomain;
  using coeff_map_t = std::unordered_map<variable_t, std::set<unsigned>>;
  using bound_t = ikos::bound<number_t>;

  base_domain_t m_base_absval;
  base_domain_t m_ext_absval;
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
  coefficient_map_t m_coeff_map;
#endif
  boost::optional<variable_t> counter;

  // ============================================================
  // Ghost variable APIs
  // ============================================================

  variable_t get_ghost_var(const variable_t &v, unsigned coefficient) {
    if (coefficient == 0) {
      CRAB_ERROR("Coefficient must be > 0");
    } else if (coefficient == 1) {
      return v;
    }

#if TVPI_DBM_FIXED_COEFFICIENTS == 0
    auto it = m_coeff_map.find(v);
    if (it != m_coeff_map.end()) {
      it->second.insert(coefficient);
    } else {
      m_coeff_map.insert({v, coefficient});
    }
#endif
    auto &vfac = const_cast<varname_t *>(&(v.name()))->get_var_factory();

    variable_t coeff_v(vfac.get_or_insert_varname(v.name(), coefficient),
                       v.get_type());
    return coeff_v;
  }

  variable_t get_original_var(const variable_t &v) {
    auto &vfac = const_cast<varname_t *>(&(v.name()))->get_var_factory();
    auto ret = vfac.find_original_varname(v.name());
    return variable_t(ret.first, v.get_type());
  }

  unsigned convert(const number_t &coefficient) const {
    if (!coefficient.fits_int64()) {
      CRAB_ERROR("Coefficient must be an 64 bit integer");
    }
    int64_t coeff = static_cast<int64_t>(coefficient);
    if (coeff < 1) {
      CRAB_ERROR("Coefficient must be positive");
    }
    return static_cast<unsigned>(coeff);
  }

  variable_t get_ghost_var(const variable_t &v, const number_t &coefficient) {
    return get_ghost_var(v, convert(coefficient));
  }

  boost::optional<variable_t> find_ghost_var(const variable_t &v,
                                             unsigned coefficient) const {
    if (coefficient == 0) {
      return boost::none;
    } else if (coefficient == 1) {
      return v;
    }
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
    else if (m_coeff_map.empty()) {
      return boost::none;
    } else {
      auto it = m_coeff_map.find(v);
      if (it != m_coeff_map.end()) {
        auto &coeff_set = it->second;
        if (coeff_set.find(coefficient) != coeff_set.end()) {
          auto &vfac = const_cast<varname_t *>(&(v.name()))->get_var_factory();
          variable_t coeff_v(vfac.get_or_insert_varname(v.name(), coefficient),
                             v.get_type());
          return coeff_v;
        }
      }
      return boost::none;
    }
#else
    else if (tvpi_utils::find(crab_domain_params_man::get().coefficients(),
                              coefficient) != boost::none) {
      auto &vfac = const_cast<varname_t *>(&(v.name()))->get_var_factory();
      variable_t coeff_v(vfac.get_or_insert_varname(v.name(), coefficient),
                         v.get_type());
      return coeff_v;
    } else {
      return boost::none;
    }
#endif
  }

  boost::optional<variable_t>
  find_ghost_var(const variable_t &v, const number_t &coefficient) const {
    return find_ghost_var(v, convert(coefficient));
  }

  // ============================================================
  // Linear Constraints APIs
  // ============================================================

  linear_expression_t rewrite_linear_expression(const linear_expression_t &e) {
    /**
     *
     * Given c1*x1 + c2*x2 +... + k, rewrite into
     *       c1x1 + c2x2 +... + kc
     **/
    linear_expression_t res;
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const variable_t &v = (*it).second;
      const number_t &coeff = (*it).first;
      if (coeff == 0) {
        continue;
      } else if (coeff > 0) {
#if TVPI_DBM_FIXED_COEFFICIENTS
        if (find_ghost_var(v, coeff) == boost::none) { // give up
          return e;
        }
#endif
        res = res + get_ghost_var(v, coeff);
      } else if (coeff < 0) {
#if TVPI_DBM_FIXED_COEFFICIENTS
        if (find_ghost_var(v, -coeff) == boost::none) { // give up
          return e;
        }
#endif
        res = res - get_ghost_var(v, -coeff);
      } else { // give up
        res = res + coeff * v;
      }
    }
    res = res + e.constant();
    return res;
  }

  linear_constraint_t
  rewrite_linear_constraint(const linear_constraint_t &cst) {
    /**
     *
     * Given c1*x1 + c2*x2 +... <= k rewrite into
     *       c1'x1 + c2'x2 + ... <= k', if all coefficients have common divisor,
     *  c1' = c1 / d, c2' = c2 / d, ..., k' = k / d
     **/
    linear_expression_t res;
    const linear_expression_t e = cst.expression();
    number_t d(0);

    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const number_t &coeff = (*it).first;
      number_t abs_coeff = coeff < 0 ? -coeff : coeff;
      d = d == number_t(0) ? abs_coeff : tvpi_utils::gcd(d, abs_coeff);
      if (d == number_t(1)) {
        break;
      }
    }
    const number_t &k = e.constant();
    number_t abs_k = k < 0 ? -k : k;
    d = d == number_t(0) ? abs_k : tvpi_utils::gcd(d, abs_k);
    if (d != number_t(1) && d != number_t(0)) {
      for (auto it = e.begin(), et = e.end(); it != et; ++it) {
        const variable_t &v = (*it).second;
        const number_t &coeff = (*it).first;
        number_t abs_coeff = coeff < 0 ? -coeff : coeff;
        bool neg = coeff < 0;
        number_t new_coeff = abs_coeff / d;
#if TVPI_DBM_FIXED_COEFFICIENTS
        if (find_ghost_var(v, new_coeff) == boost::none) { // give up
          return cst;
        }
#endif
        auto gv = get_ghost_var(v, new_coeff);
        res = neg ? res - gv : res + gv;
      }
      number_t new_k = k / d;
      res = res + new_k;
    } else {
      // rewrite to same form but using ghost variable
      for (auto it = e.begin(), et = e.end(); it != et; ++it) {
        const variable_t &v = (*it).second;
        const number_t &coeff = (*it).first;
        number_t abs_coeff = coeff < 0 ? -coeff : coeff;
        bool neg = coeff < 0;
#if TVPI_DBM_FIXED_COEFFICIENTS
        if (find_ghost_var(v, abs_coeff) == boost::none) { // give up
          return cst;
        }
#endif
        auto gv = get_ghost_var(v, abs_coeff);
        res = neg ? res - gv : res + gv;
      }
      res = res + k;
    }
    return linear_constraint_t(res, cst.kind());
  }

  boost::optional<linear_constraint_t>
  try_rewrite_linear_constraint(const linear_constraint_t &cst) const {

    linear_expression_t res;
    const linear_expression_t e = cst.expression();
    number_t d(0);

    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const number_t &coeff = (*it).first;
      number_t abs_coeff = coeff < 0 ? -coeff : coeff;
      d = d == number_t(0) ? abs_coeff : tvpi_utils::gcd(d, abs_coeff);
      if (d == number_t(1)) {
        break;
      }
    }
    const number_t k = e.constant();
    number_t abs_k = k < 0 ? -k : k;
    d = d == number_t(0) ? abs_k : tvpi_utils::gcd(d, abs_k);
    if (d != number_t(1) && d != number_t(0)) {
      for (auto it = e.begin(), et = e.end(); it != et; ++it) {
        const variable_t &v = (*it).second;
        const number_t &coeff = (*it).first;
        number_t abs_coeff = coeff < 0 ? -coeff : coeff;
        number_t new_coeff = abs_coeff / d;
        bool neg = coeff < 0;
        auto gv = find_ghost_var(v, new_coeff);
        if (!gv)
          return boost::none;
        res = neg ? res - (*gv) : res + (*gv);
      }
      number_t new_k = k / d;
      res = res + new_k;
    } else {
      // rewrite to same form but using ghost variable
      for (auto it = e.begin(), et = e.end(); it != et; ++it) {
        const variable_t &v = (*it).second;
        const number_t &coeff = (*it).first;
        number_t abs_coeff = coeff < 0 ? -coeff : coeff;
        bool neg = coeff < 0;
        auto gv = find_ghost_var(v, abs_coeff);
        if (!gv)
          return boost::none;
        res = neg ? res - (*gv) : res + (*gv);
      }
      res = res + k;
    }
    return linear_constraint_t(res, cst.kind());
  }

  void rewrite_assign(const variable_t &x, const linear_expression_t &e,
                      bool weak) {
    if (e.is_constant()) {
      // (1) e == constant
      CRAB_LOG("tvpi-dbm-assign",
               crab::outs() << "cannot rewrite: " << x << " := " << e << "\n");
      return;
    } else if (e.size() == 1) {
      // (2) e == b*y +/- c
      auto it = e.begin();
      const number_t &b = (*it).first;
      const variable_t &y = (*it).second;
      number_t abs_b = b > number_t(0) ? b : -b;
      bool neg_b = (b < number_t(0));
      const number_t c = e.constant();
      number_t abs_c = c > number_t(0) ? c : -c;
      if (abs_b != number_t(0) && abs_b != number_t(1)) {
#if TVPI_DBM_FIXED_COEFFICIENTS
        if (find_ghost_var(y, abs_b) == boost::none) { // give up
          CRAB_LOG("tvpi-dbm-assign", crab::outs() << "cannot rewrite: " << x
                                                   << " := " << e << "\n");
          return;
        }
#endif
        auto by = get_ghost_var(y, abs_b);
        if (c == number_t(0)) {
          // rewrite("x := b*y") = "x := by"
          linear_expression_t e1 =
              neg_b ? linear_expression_t(-by) : linear_expression_t(by);
          CRAB_LOG("tvpi-dbm-assign", crab::outs()
                                          << "processing rewritten " << x
                                          << " := " << e1 << "\n");
          if (!weak) {
            m_ext_absval.assign(x, e1);
          } else {
            m_ext_absval.weak_assign(x, e1);
          }
        } else {
          linear_expression_t e1 = neg_b ? linear_expression_t(-by + c)
                                         : linear_expression_t(by + c);
          CRAB_LOG("tvpi-dbm-assign", crab::outs()
                                          << "processing rewritten " << x
                                          << " := " << e1 << "\n");
          if (!weak) {
            m_ext_absval.assign(x, e1);
          } else {
            m_ext_absval.weak_assign(x, e1);
          }
          // rewrite("x := b*y +/- c") = "x / d := (b * y  +/- c) / d "
          // compute gcd(b, c) = d
          number_t d = tvpi_utils::gcd(abs_b, abs_c);
          // if d > 1, let b' = b / d, c' = c / d
          if (d > number_t(1)) {
            number_t new_abs_b = abs_b / d;
            number_t new_c = c / d;
            // => find b' * y + / - c' => assign x by (b' * y + / - c') * d
            if (new_abs_b > number_t(1)) {
#if TVPI_DBM_FIXED_COEFFICIENTS
              if (find_ghost_var(y, new_abs_b) == boost::none) { // give up
                CRAB_LOG("tvpi-dbm-assign", crab::outs()
                                                << "cannot rewrite: " << x
                                                << " := " << e << "\n");
                return;
              }
#endif
              auto by = get_ghost_var(y, new_abs_b);
              linear_expression_t e2 = neg_b ? linear_expression_t(-by + new_c)
                                             : linear_expression_t(by + new_c);
              CRAB_LOG("tvpi-dbm-assign", crab::outs()
                                              << "processing rewritten " << x
                                              << " := " << e2 << "\n");
              if (!weak) {
                m_ext_absval.assign(x, e2);
              } else {
                m_ext_absval.weak_assign(x, e2);
              }
              m_ext_absval.apply(OP_MULTIPLICATION, x, x, d);
            } else { // new_abs_b == 1
              linear_expression_t e2 = neg_b ? linear_expression_t(-y + new_c)
                                             : linear_expression_t(y + new_c);
              CRAB_LOG("tvpi-dbm-assign", crab::outs()
                                              << "processing rewritten " << x
                                              << " := " << e2 << "\n");
              if (!weak) {
                m_ext_absval.assign(x, e2);
              } else {
                m_ext_absval.weak_assign(x, e2);
              }
              m_base_absval.apply(OP_MULTIPLICATION, x, x, d);
            }
          }
        }
      }
    } else if (e.size() == 2) {
      // (3) e == b*y +/- c*z +/- d
      auto it = e.begin();
      const number_t &b = (*it).first;
      number_t abs_b = b > number_t(0) ? b : -b;
      bool neg_b = (b < number_t(0));
      const variable_t &y = (*it).second;
      it = boost::next(it);
      const number_t &c = (*it).first;
      number_t abs_c = c > number_t(0) ? c : -c;
      bool neg_c = (c < number_t(0));
      const variable_t &z = (*it).second;
      number_t d = e.constant();
      number_t abs_d = d > number_t(0) ? d : -d;
      if (b == number_t(0)) {
        linear_expression_t e1 = linear_expression_t(c * z + d);
        if (!weak) {
          assign(x, e1);
        } else {
          weak_assign(x, e1);
        }
      } else if (c == number_t(0)) {
        linear_expression_t e1 = linear_expression_t(b * y + d);
        if (!weak) {
          assign(x, e1);
        } else {
          weak_assign(x, e1);
        }
      } else if (abs_b != number_t(1) && abs_c != number_t(1)) {
#if TVPI_DBM_FIXED_COEFFICIENTS
        if (find_ghost_var(y, abs_b) == boost::none ||
            find_ghost_var(z, abs_c) == boost::none) { // give up
          CRAB_LOG("tvpi-dbm-assign", crab::outs() << "cannot rewrite: " << x
                                                   << " := " << e << "\n");
          return;
        }
#endif
        auto by = get_ghost_var(y, abs_b);
        auto cz = get_ghost_var(z, abs_c);
        linear_expression_t e1 = linear_expression_t(neg_b ? -by : by) +
                                 linear_expression_t(neg_c ? -cz : cz) + d;
        CRAB_LOG("tvpi-dbm-assign", crab::outs() << "processing rewritten " << x
                                                 << " := " << e1 << "\n");
        if (!weak) {
          m_ext_absval.assign(x, e1);
        } else {
          m_ext_absval.weak_assign(x, e1);
        }
        // rewrite("x := b*y +/- c*z +/- d") = "x / d := (b * y  +/- c * z +/-
        // d) / d "
        // TODO: ignore this for now.
        // number_t d = tvpi_utils::gcd3(abs_b, abs_c, abs_d);
      }
    } else {
      // (4) e == general form
      auto e1 = rewrite_linear_expression(e);
      if (!e1.equal(e)) {
        if (!weak) {
          m_ext_absval.assign(x, e1);
        } else {
          m_ext_absval.weak_assign(x, e1);
        }
      }
    }
  }

  void rewrite_apply(arith_operation_t op, const variable_t &x,
                     const variable_t &y, number_t z, unsigned coefficient) {
    assert(coefficient > 1);

    if (find_ghost_var(x, coefficient) == boost::none ||
        find_ghost_var(y, coefficient) == boost::none) {
      return;
    }
    variable_t ghost_x = get_ghost_var(x, coefficient);
    variable_t ghost_y = get_ghost_var(y, coefficient);
    number_t tracked_coefficient(coefficient);
    switch (op) {
    case OP_MULTIPLICATION:
    case OP_SDIV:
    case OP_UDIV: // x := y * z or x := y / z
      // rewrite to x*COEF := y*COEF */div z
      m_ext_absval.apply(op, ghost_x, ghost_y, z);
      break;
    case OP_ADDITION:
    case OP_SUBTRACTION: // x := y + z or x := y - z
      // rewrite to x*COEF := y*COEF +/- z*COEF
      m_ext_absval.apply(op, ghost_x, ghost_y, z * tracked_coefficient);
      break;
    default:
      break;
    }
  }

  void rewrite_apply(arith_operation_t op, const variable_t &x,
                     const variable_t &y, const variable_t &z,
                     unsigned coefficient) {
    assert(coefficient > 1);
    if (find_ghost_var(x, coefficient) == boost::none ||
        find_ghost_var(y, coefficient) == boost::none ||
        find_ghost_var(z, coefficient) == boost::none) {
      return;
    }
    // rewrite("x := y op z") = "x*COEF := y*COEF op z*COEF"
    variable_t ghost_x = get_ghost_var(x, coefficient);
    variable_t ghost_y = get_ghost_var(y, coefficient);
    variable_t ghost_z = get_ghost_var(z, coefficient);
    m_ext_absval.apply(op, ghost_x, ghost_y, ghost_z);
  }

  // ============================================================
  // TVPI APIs
  // ============================================================

  bool add_tvpi_constraint(const variable_t &ax, const variable_t &by,
                           const number_t &c) {
    // add ax - by <= c
    auto oldcopt = m_ext_absval.difference_bound(by, ax);
    if (!oldcopt || c < *oldcopt) { // if new bound is tighter
      m_ext_absval += (ax - by <= c);
      return false;
    }
    return true;
  }

  bool add_utvpi_constraint(const variable_t &x, const variable_t &y,
                            const number_t &c) {
    // add x - y <= c
    auto oldcopt = m_base_absval.difference_bound(y, x);
    if (!oldcopt || c < *oldcopt) { // if new bound is tighter
      m_base_absval += (x - y <= c);
      return false;
    }
    return true;
  }

  bool add_ub_constraint(const variable_t &x, const number_t &ub) {
    // add x <= ub
    bound_t x_ub = m_base_absval.at(x).ub();
    bound_t new_ub = bound_t(ub);
    if (new_ub < x_ub) { // if new bound is tighter
      m_base_absval += (x <= ub);
      return false;
    }
    return true;
  }

  bool add_lb_constraint(const variable_t &x, const number_t &lb) {
    // add -x <= lb
    bound_t x_lb = m_base_absval.at(x).lb();
    bound_t new_lb = bound_t(-lb);
    if (new_lb > x_lb) { // if new bound is tighter
      m_base_absval += (-x <= lb);
      return false;
    }
    return true;
  }

  std::pair<unsigned, number_t> normalize_tvpi(const unsigned &a,
                                               const number_t &c) const {
    // Normalize TVPI constraints
    // For each constraint ax <= c, normalize it to x <= c' where c' = c / a
    return {1, c / number_t(a)};
  }

  std::pair<unsigned, number_t> normalize_tvpi(const unsigned &a,
                                               const unsigned &b,
                                               const number_t &c) const {
    // Normalize TVPI constraints
    // For each constraint ax - by <= c, normalize it to a'x - b'y <= c'
    // where a' = a / gcd(a, b), b' = b / gcd(a, b), c' = c / gcd(a, b)
    unsigned gcd = tvpi_utils::gcd(a, b);
    return {gcd, c / number_t(gcd)};
  }

  std::tuple<unsigned, unsigned, number_t>
  resultant(const unsigned &a, const variable_t &x, const unsigned &b,
            const variable_t &y, const number_t &c, const unsigned &d,
            const unsigned &e, const boost::optional<variable_t> &z,
            const number_t &f) const {

    // Given two inequalities ax - by <= c and dy - ez <= f, compute the
    // resultant inequalities by eliminating the middle variable y.
    unsigned gcd = tvpi_utils::gcd(b, d);
    unsigned lambda1 = d / gcd;
    unsigned lambda2 = b / gcd;
    // Now we have l1 * ax - l2 * ez <= l1 * c + l2 * f
    std::string z_name = z ? z->name().str() : "v0";
    CRAB_LOG("tvpi-dbm-resultant",
             crab::outs() << lambda1 << "*" << a << x << " - " << lambda2 << "*"
                          << e << z_name << "<=" << lambda1 << "*" << c << "+"
                          << lambda2 << "*" << f << "\n");

    auto c_p = c * number_t(lambda1) + f * number_t(lambda2);
    // Normalize the result
    if (z == boost::none || x == *z) {
      int a_p = (z == boost::none) ? lambda1 * a : lambda1 * a - lambda2 * e;
      if (a_p == 0) {
        return {0, 0, c_p};
      }
      unsigned abs_a = static_cast<unsigned>(std::abs(a_p));
      bool neg = (a_p) < 0;
      auto ret = normalize_tvpi(abs_a, c_p);
      unsigned new_a = neg ? 0 : 1;
      unsigned new_b = neg ? 1 : 0;
      number_t new_c = ret.second;
      return {new_a, new_b, new_c};
    } else {
      auto ret = normalize_tvpi(lambda1 * a, lambda2 * e, c_p);
      unsigned gcd2 = ret.first;
      unsigned new_a = lambda1 * a / gcd2;
      unsigned new_b = lambda2 * e / gcd2;
      number_t new_c = ret.second;
      return {new_a, new_b, new_c};
    }
  }

  std::tuple<unsigned, unsigned, number_t>
  resultant(const unsigned &b, const variable_t &y, const unsigned &a,
            const variable_t &x, const number_t &c, const unsigned &e,
            const boost::optional<variable_t> &z, const unsigned &d,
            const number_t &f) const {
    // Given two inequalities by - ax <= c and ez - dy <= f, compute the
    // resultant inequalities by eliminating the middle variable y.
    unsigned gcd = tvpi_utils::gcd(b, d);
    unsigned lambda1 = d / gcd;
    unsigned lambda2 = b / gcd;
    // Now we have l2 * ez - l1 * ax <= l1 * c + l2 * f
    std::string z_name = z ? z->name().str() : "v0";
    CRAB_LOG("tvpi-dbm-resultant",
             crab::outs() << lambda2 << "*" << e << z_name << " - " << lambda1
                          << "*" << a << x << "<=" << lambda1 << "*" << c << "+"
                          << lambda2 << "*" << f << "\n");

    auto c_p = c * number_t(lambda1) + f * number_t(lambda2);
    // Normalize the result
    if (z == boost::none || x == *z) {
      int a_p = (z == boost::none) ? -lambda1 * a : lambda2 * e - lambda1 * a;
      if (a_p == 0) {
        return {0, 0, c_p};
      }
      unsigned abs_a = static_cast<unsigned>(std::abs(a_p));
      bool neg = (a_p) < 0;
      auto ret = normalize_tvpi(abs_a, c_p);
      unsigned new_a = neg ? 0 : 1;
      unsigned new_b = neg ? 1 : 0;
      number_t new_c = ret.second;
      return {new_a, new_b, new_c};
    } else {
      auto ret = normalize_tvpi(lambda2 * e, lambda1 * a, c_p);
      unsigned gcd2 = ret.first;
      unsigned new_a = lambda2 * e / gcd2;
      unsigned new_b = lambda1 * a / gcd2;
      number_t new_c = ret.second;
      return {new_a, new_b, new_c};
    }
  }

#if TVPI_DBM_FIXED_COEFFICIENTS == 0
  tvpi_dbm_domain(base_domain_t &&base, base_domain_t &&extd,
                  coefficient_map_t &&coeff_map)
      : m_base_absval(std::move(base)), m_ext_absval(std::move(extd)),
        m_coeff_map(std::move(coeff_map)) {
    // tvpi_reduce();
  }
#else
  tvpi_dbm_domain(base_domain_t &&base, base_domain_t &&extd)
      : m_base_absval(std::move(base)), m_ext_absval(std::move(extd)) {
    // tvpi_reduce();
  }
#endif

public:
  DEFAULT_SELECT(tvpi_dbm_domain_t)
  BOOL_OPERATIONS_NOT_IMPLEMENTED(tvpi_dbm_domain_t)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(tvpi_dbm_domain_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(tvpi_dbm_domain_t)

  tvpi_dbm_domain() {}

  tvpi_dbm_domain(const tvpi_dbm_domain_t &o) = default;
  tvpi_dbm_domain(tvpi_dbm_domain_t &&o) = default;
  tvpi_dbm_domain_t &operator=(const tvpi_dbm_domain_t &o) = default;
  tvpi_dbm_domain_t &operator=(tvpi_dbm_domain_t &&o) = default;

  bool is_asc_phase() const override {
    return m_base_absval.is_asc_phase() && m_ext_absval.is_asc_phase();
  }

  void set_phase(bool is_ascending) override {
    m_base_absval.set_phase(is_ascending);
    m_ext_absval.set_phase(is_ascending);
  }

  void set_to_top() override {
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
    m_coeff_map.set_to_top();
#endif
    m_base_absval.set_to_top();
    m_ext_absval.set_to_top();
  }

  void set_to_bottom() override {
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
    m_coeff_map.set_to_bottom();
#endif
    m_base_absval.set_to_bottom();
    m_ext_absval.set_to_bottom();
  }

  tvpi_dbm_domain_t make_bottom() const override {
    tvpi_dbm_domain_t res;
    res.set_to_bottom();
    return res;
  }

  tvpi_dbm_domain_t make_top() const override {
    tvpi_dbm_domain_t res;
    return res;
  }

  bool is_bottom() const override {
    bool res = m_base_absval.is_bottom() || m_ext_absval.is_bottom();
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
    res = res || m_coeff_map.is_bottom();
#endif
    return res;
  }

  bool is_top() const override {
    bool res = m_base_absval.is_top() && m_ext_absval.is_top();
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
    res = res && m_coeff_map.is_top();
#endif
    return res;
  }

  // void normalize_dbms() {
  //   // This is a top level function to normalize tvpi constraints.
  //   // However, we wish not to run this process in the end since we always
  //   keep
  //   // inequalities after normalized.
  //   // Just in case we need in the future, this is the implementation.
  //   for (auto it1 = m_coeff_map.begin(); it1 != m_coeff_map.end(); ++it1) {
  //     const variable_t &x = it1->first;
  //     const coefficient_set_t &a_coeffs = it1->second;
  //     for (auto ita = a_coeffs.begin(); ita != a_coeffs.end(); ++ita) {
  //       auto gax = get_ghost_var(x, *ita);
  //       for (auto it2 = std::next(it1); it2 != m_coeff_map.end(); ++it2) {
  //         const variable_t &y = it2->first;
  //         const coefficient_set_t &b_coeffs = it2->second;
  //         for (auto itb = b_coeffs.begin(); itb != b_coeffs.end(); ++itb) {
  //           auto gby = get_ghost_var(y, *itb);
  //           auto copt = m_ext_absval.difference_bound(gby, gax);
  //           if (copt) {
  //             auto ret = normalize_tvpi(*ita, *itb, *copt);
  //             auto gcd = ret.first;
  //             unsigned new_a = *ita / gcd;
  //             unsigned new_b = *itb / gcd;
  //             number_t new_c = ret.second;
  //             if (new_a == 1 && new_b == 1) {
  //               m_base_absval += (x - y <= new_c);
  //             } else {
  //               auto gax_new = get_ghost_var(x, new_a);
  //               auto gby_new = get_ghost_var(y, new_b);
  //               m_ext_absval += (gax_new - gby_new <= new_c);
  //             }
  //             // TODO: Remove old constraint
  //           }
  //         }
  //       }
  //     }
  //   }
  // }

  void filter_dbms() {
    // This function is doing some inequalities removal based on the coefficient
    // map and redundant inequality (if any).

    // Propagate inferred inequalties from ext to base if base should track
    // them.
    auto ext_vars = m_ext_absval.vars();
    auto base_vars = m_base_absval.vars();
    variable_vector_t bases;
    bases.reserve(ext_vars.size());
    for (auto &ev : ext_vars) {
      if (tvpi_utils::find(base_vars, ev) != boost::none) {
        bases.push_back(ev);
      }
    }
    auto extd2 = m_ext_absval;
    extd2.project(bases);
    m_base_absval &= extd2;

// Remove inequalities if coefficients are not tracked.
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
    ext_vars = m_ext_absval.vars();
    base_vars = m_base_absval.vars();
    std::unordered_set<variable_t> to_remove(ext_vars.begin(), ext_vars.end());
    for (auto it1 = m_coeff_map.begin(); it1 != m_coeff_map.end(); ++it1) {
      const variable_t &x = it1->first;
      const coefficient_set_t &a_coeffs = it1->second;
      for (auto ita = a_coeffs.begin(); ita != a_coeffs.end(); ++ita) {
        auto gax = get_ghost_var(x, *ita);
        if (to_remove.find(gax) != to_remove.end()) {
          to_remove.erase(gax);
        }
      }
    }

    for (auto it = to_remove.begin(); it != to_remove.end();) {
      if (tvpi_utils::find(base_vars, *it) || tvpi_utils::find(ext_vars, *it)) {
        it = to_remove.erase(it);
      } else {
        ++it;
      }
    }
    CRAB_LOG("tvpi-dbm-reduce2", crab::outs() << "Remove dimensions: ";
             tvpi_utils::print_set(crab::outs(), to_remove);
             crab::outs() << "\n";);
    m_ext_absval.forget(variable_vector_t(to_remove.begin(), to_remove.end()));
#endif
  }

  void prune_coefficients() {
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
    // This function removes coefficients that current DBM lost those
    // dimensions.
    auto ext_vars = m_ext_absval.vars();
    for (auto it = m_coeff_map.begin(); it != m_coeff_map.end();) {
      const variable_t &x = it->first;
      coefficient_set_t &a_coeffs = it->second;
      for (auto ita = a_coeffs.begin(); ita != a_coeffs.end();) {
        auto gax = get_ghost_var(x, *ita);
        if (tvpi_utils::find(ext_vars, gax) == boost::none) {
          ita = a_coeffs.erase(ita);
        } else {
          ++ita;
        }
      }
      if (a_coeffs.empty()) {
        it = m_coeff_map.erase(it);
      } else {
        ++it;
      }
    }
#endif
  }

  void tvpi_reduce() {
    CRAB_LOG("tvpi-dbm-reduce", crab::outs()
                                    << "Before reduction: " << *this << "\n");

    auto for_each_coefficient =
        [](const coefficient_set_t &s,
           const std::function<void(unsigned)> &process) {
          process(1);
          for (const auto &i : s) {
            process(i);
          }
        };

    auto base_vars = m_base_absval.vars();
    CRAB_LOG("tvpi-dbm-reduce2", crab::outs() << "base dimensions: ";
             tvpi_utils::print_vector(crab::outs(), base_vars);
             crab::outs() << "\n";);
    auto ext_vars = m_ext_absval.vars();
    variable_vector_t allvars;
    allvars.reserve(base_vars.size() + ext_vars.size());
    allvars.insert(allvars.end(), base_vars.begin(), base_vars.end());
    allvars.insert(allvars.end(), ext_vars.begin(), ext_vars.end());
    CRAB_LOG("tvpi-dbm-reduce2", crab::outs() << "extend dimensions: ";
             tvpi_utils::print_vector(crab::outs(), ext_vars);
             crab::outs() << "\n";);
    auto e_set = std::unordered_set<variable_t>(allvars.begin(), allvars.end());
    for (auto &v : allvars) {
      auto ov = get_original_var(v);
      if (ov != v) {
        e_set.erase(v);
        e_set.insert(ov);
      }
    }
    auto traverse_vars = variable_vector_t(e_set.begin(), e_set.end());
    CRAB_LOG("tvpi-dbm-reduce2", crab::outs() << "Do reductions on: ";
             tvpi_utils::print_vector(crab::outs(), traverse_vars);
             crab::outs() << "\n";);
    coefficient_set_t no_coeff = {};
    for (auto it1 = traverse_vars.cbegin(); it1 != traverse_vars.cend();
         ++it1) {
      const variable_t &x = *it1;
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      auto itc1 = m_coeff_map.find(x);
      const coefficient_set_t &a_coeffs =
          itc1 != m_coeff_map.end() ? itc1->second : no_coeff;
#else
      const coefficient_set_t &a_coeffs =
          crab_domain_params_man::get().coefficients();
#endif
      for (auto it2 = traverse_vars.cbegin(); it2 != traverse_vars.cend();
           ++it2) {
        const variable_t &y = *it2;
        if (x == y) {
          continue;
        }
        if (tvpi_utils::find(base_vars, y) == boost::none)
          continue;
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
        auto itc2 = m_coeff_map.find(y);
        const coefficient_set_t &b_coeffs =
            itc2 != m_coeff_map.end() ? itc2->second : no_coeff;
#else
        const coefficient_set_t &b_coeffs =
            crab_domain_params_man::get().coefficients();
#endif
        for_each_coefficient(a_coeffs, [&](unsigned a) {
          for_each_coefficient(b_coeffs, [&](unsigned b) {
            if (a != 1 or
                b != 1) { // avoid infer inequalities that closure does
              auto gax = get_ghost_var(x, a);
              auto gby = get_ghost_var(y, b);
              auto copt = m_ext_absval.difference_bound(gby, gax);
              if (copt) {
                for (auto &z : base_vars) {
                  if (z == y) {
                    continue;
                  }
                  auto fopt = m_base_absval.difference_bound(z, y);
                  if (fopt == boost::none) {
                    continue;
                  }
                  // ax - by <= c && y - z <= f => a'x - b'z <= c'
                  CRAB_LOG("tvpi-dbm-reduce2",
                           crab::outs() << "resultant(" << gax << "-" << gby
                                        << " <= " << *copt << ", " << y << "-"
                                        << z << " <= " << *fopt
                                        << "), eliminating " << y << "\n");
                  auto ret = resultant(a, x, b, y, *copt, 1, 1, z, *fopt);
                  auto new_a = std::get<0>(ret);
                  auto new_b = std::get<1>(ret);
                  auto new_c = std::get<2>(ret);
                  bool skip = true;
                  if (new_a == 0 && new_b == 0) {
                    if (new_c < 0) { // 0 <= c where c < 0, UNSAT
                      set_to_bottom();
                      return;
                    }
                  } else if (new_a == 1 && new_b == 1) { // x - z <= c'
                    skip = add_utvpi_constraint(x, z, new_c);
                  } else if (new_a == 1 && new_b == 0) { // x <= c'
                    skip = add_ub_constraint(x, new_c);
                  } else if (new_a == 0 && new_b == 1) { // -z <= c'
                    skip = add_lb_constraint(z, new_c);
                  } else { // a'x - b'z <= c' with a' > 1, b' > 1
                    if (find_ghost_var(x, new_a) && find_ghost_var(z, new_b)) {
                      auto gax_new = get_ghost_var(x, new_a);
                      auto gbz_new = get_ghost_var(z, new_b);
                      if (tvpi_utils::find(ext_vars, gax_new) ||
                          tvpi_utils::find(ext_vars, gbz_new)) {
                        skip = add_tvpi_constraint(gax_new, gbz_new, new_c);
                      }
                    }
                  }
                  CRAB_LOG("tvpi-dbm-reduce2",
                           crab::outs() << "=>>>" << new_a << "*" << x << "-"
                                        << new_b << "*" << z << "<=" << new_c
                                        << (skip ? ", skip" : "") << "\n");
                }
                // ax - by <= c && y <= f => x <= c'
                auto bnds = m_base_absval.at(y);
                if (auto fopt = bnds.ub().number()) {
                  CRAB_LOG("tvpi-dbm-reduce2",
                           crab::outs() << "resultant(" << gax << "-" << gby
                                        << " <= " << *copt << ", " << y
                                        << " <= " << *fopt << "), eliminating "
                                        << y << "\n");
                  auto ret =
                      resultant(a, x, b, y, *copt, 1, 1, boost::none, *fopt);
                  auto new_a = std::get<0>(ret);
                  auto new_c = std::get<2>(ret);
                  bool skip = add_ub_constraint(x, new_c);
                  CRAB_LOG("tvpi-dbm-reduce2",
                           crab::outs()
                               << "=>>>" << new_a << "*" << x << " <= " << new_c
                               << (skip ? ", skip" : "") << "\n");
                }
                // ax - by <= c && -x <= f => -y <= c'
                bnds = m_base_absval.at(x);
                if (auto fopt = bnds.lb().number()) {
                  CRAB_LOG("tvpi-dbm-reduce2",
                           crab::outs() << "resultant(" << gax << "-" << gby
                                        << " <= " << *copt << ", -" << x
                                        << " <= " << (-*fopt)
                                        << "), eliminating " << x << "\n");
                  auto ret =
                      resultant(a, x, b, y, *copt, 1, boost::none, 1, (-*fopt));
                  auto new_b = std::get<1>(ret);
                  auto new_c = std::get<2>(ret);
                  bool skip = add_lb_constraint(y, new_c);
                  CRAB_LOG("tvpi-dbm-reduce2",
                           crab::outs()
                               << "=>>>"
                               << "-" << new_b << "*" << y << " <= " << new_c
                               << (skip ? ", skip" : "") << "\n");
                }
              } else if (x != y && !(a == 1 && b == 1) && counter &&
                         (x == *counter || y == *counter)) {
                // see whether we can recover this
                // x <= c && -y <= f => ax - by <= c'
                // Now, we skip lots of cases here since it will create many
                // inequalities
                auto bndsx = m_base_absval.at(x);
                auto bndsy = m_base_absval.at(y);
                if (bndsx.ub().is_finite() && bndsy.lb().is_finite()) {
                  auto xub = bndsx.ub().number();
                  auto ylb = bndsy.lb().number();
                  if (xub && ylb) {
                    CRAB_LOG("tvpi-dbm-reduce2",
                             crab::outs()
                                 << "create(" << a << "*" << x << " <= " << a
                                 << "*" << *xub << ", -" << b << "*" << y
                                 << " <= -" << b << "*" << *ylb << ")\n");
                    number_t c_p = number_t(a) * (*xub) - number_t(b) * (*ylb);
                    bool skip = add_tvpi_constraint(gax, gby, c_p);
                    CRAB_LOG("tvpi-dbm-reduce2",
                             crab::outs()
                                 << "=>>>" << gax << "-" << gby << " <= " << c_p
                                 << (skip ? ", skip" : "") << "\n");
                  }
                }
              }
              // explore more cases
              copt = m_ext_absval.difference_bound(gax, gby);
              if (copt) {
                for (auto &z : base_vars) {
                  if (z == y) {
                    continue;
                  }
                  auto fopt = m_base_absval.difference_bound(y, z);
                  if (fopt == boost::none) {
                    continue;
                  }
                  // by - ax <= c && z - y <= f => a'z - b'x <= c'
                  CRAB_LOG("tvpi-dbm-reduce2",
                           crab::outs() << "resultant(" << gby << "-" << gax
                                        << " <= " << *copt << ", " << z << "-"
                                        << y << " <= " << *fopt
                                        << "), eliminating " << y << "\n");
                  auto ret = resultant(b, y, a, x, *copt, 1, z, 1, *fopt);
                  auto new_a = std::get<0>(ret);
                  auto new_b = std::get<1>(ret);
                  auto new_c = std::get<2>(ret);
                  bool skip = true;
                  if (new_a == 0 && new_b == 0) {
                    if (new_c < 0) { // 0 <= c where c < 0, UNSAT
                      set_to_bottom();
                      return;
                    }
                  } else if (new_a == 1 && new_b == 1) { // z - x <= c'
                    skip = add_utvpi_constraint(z, x, new_c);
                  } else if (new_a == 1 && new_b == 0) { // z <= c'
                    skip = add_ub_constraint(z, new_c);
                  } else if (new_a == 0 && new_b == 1) { // -x <= c'
                    skip = add_lb_constraint(x, new_c);
                  } else { // a'z - b'x <= c' with a' > 1, b' > 1
                    if (find_ghost_var(z, new_a) && find_ghost_var(x, new_b)) {
                      auto gaz_new = get_ghost_var(z, new_a);
                      auto gbx_new = get_ghost_var(x, new_b);
                      if (tvpi_utils::find(ext_vars, gbx_new) ||
                          tvpi_utils::find(ext_vars, gaz_new)) {
                        skip = add_tvpi_constraint(gaz_new, gbx_new, new_c);
                      }
                    }
                  }
                  CRAB_LOG("tvpi-dbm-reduce2",
                           crab::outs() << "=>>>" << new_a << "*" << z << "-"
                                        << new_b << "*" << x << " <= " << new_c
                                        << (skip ? ", skip" : "") << "\n");
                }
                // by - ax <= c && x <= f => y <= c'
                auto bnds = m_base_absval.at(x);
                if (auto fopt = bnds.ub().number()) {
                  CRAB_LOG("tvpi-dbm-reduce2",
                           crab::outs() << "resultant(" << gby << "-" << gax
                                        << " <= " << *copt << ", " << x
                                        << " <= " << *fopt << "), eliminating "
                                        << x << "\n");
                  auto ret =
                      resultant(b, y, a, x, *copt, 1, 1, boost::none, *fopt);
                  auto new_a = std::get<0>(ret);
                  auto new_c = std::get<2>(ret);
                  bool skip = add_ub_constraint(y, new_c);
                  CRAB_LOG("tvpi-dbm-reduce2",
                           crab::outs()
                               << "=>>>" << new_a << "*" << y << "<=" << new_c
                               << (skip ? ", skip" : "") << "\n");
                }
                // by - ax <= c && -y <= f => -x <= c'
                bnds = m_base_absval.at(y);
                if (auto fopt = bnds.lb().number()) {
                  CRAB_LOG("tvpi-dbm-reduce2",
                           crab::outs() << "resultant(" << gby << "-" << gax
                                        << " <= " << *copt << ", -" << y
                                        << " <= " << (-*fopt)
                                        << "), eliminating " << y << "\n");
                  auto ret =
                      resultant(b, y, a, x, *copt, 1, boost::none, 1, (-*fopt));
                  auto new_b = std::get<1>(ret);
                  auto new_c = std::get<2>(ret);
                  bool skip = add_lb_constraint(x, new_c);
                  CRAB_LOG("tvpi-dbm-reduce2",
                           crab::outs()
                               << "=>>>"
                               << "-" << new_b << "*" << x << " <= " << new_c
                               << (skip ? ", skip" : "") << "\n");
                }
              } else if (x != y && !(a == 1 && b == 1) && counter &&
                         (x == *counter || y == *counter)) {
                // see whether we can recover this
                // -x <= c && y <= f => by - ax <= c'
                auto bndsx = m_base_absval.at(x);
                auto bndsy = m_base_absval.at(y);
                if (bndsx.lb().is_finite() && bndsy.ub().is_finite()) {
                  auto xlb = bndsx.lb().number();
                  auto yub = bndsy.ub().number();
                  if (xlb && yub) {
                    CRAB_LOG("tvpi-dbm-reduce2",
                             crab::outs()
                                 << "create(-" << a << "*" << x << " <= -" << a
                                 << "*" << *xlb << ", " << b << "*" << y
                                 << " <= " << b << "*" << *yub << ")\n");
                    number_t c_p = number_t(b) * (*yub) - number_t(a) * (*xlb);
                    bool skip = add_tvpi_constraint(gby, gax, c_p);
                    CRAB_LOG("tvpi-dbm-reduce2",
                             crab::outs()
                                 << "=>>>" << gby << "-" << gax << " <= " << c_p
                                 << (skip ? ", skip" : "") << "\n");
                  }
                }
              }
            }
          });
        });
      }
    }

    CRAB_LOG("tvpi-dbm-reduce", crab::outs()
                                    << "After resultant: " << *this << "\n");

    // ext dbm may includes some inequalities that base dbm can represented.
    // Adding a filter to propagate and remove these inequalities.
    filter_dbms();
    CRAB_LOG("tvpi-dbm-reduce", crab::outs()
                                    << "After filter: " << *this << "\n");
  }

  bool operator<=(const tvpi_dbm_domain_t &other) const override {
    if (is_bottom() || other.is_top()) {
      return true;
    } else if (is_top() || other.is_bottom()) {
      return false;
    } else if (m_base_absval <= other.m_base_absval &&
               m_ext_absval <= other.m_ext_absval) {
      return true;
    } else {
      CRAB_LOG("tvpi-dbm-leq", crab::outs() << "[leq]\n"
                                            << *this << "\n<=\n"
                                            << other << "\n");
      tvpi_dbm_domain_t this2 = *this;
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      this2.m_coeff_map.meet(other.m_coeff_map);
#endif
      this2.counter = this2.counter ? this2.counter : other.counter;
      this2.tvpi_reduce();
      tvpi_dbm_domain_t other2 = other;
      other2.counter = other2.counter ? other2.counter : counter;
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      other2.m_coeff_map.meet(m_coeff_map);
#endif
      other2.tvpi_reduce();
      CRAB_LOG("tvpi-dbm-leq", crab::outs() << "[leq reduced]\n"
                                            << this2 << "\n<=\n"
                                            << other2 << "\n");
      bool res = this2.m_base_absval <= other2.m_base_absval &&
                 this2.m_ext_absval <= other2.m_ext_absval;
      CRAB_LOG("tvpi-dbm-leq",
               crab::outs() << "[res]=" << (res ? "true" : "false") << "\n");
      return res;
    }
  }

  void operator|=(const tvpi_dbm_domain_t &other) override {
    if (is_bottom() || other.is_top()) {
      *this = other;
    } else if (other.is_bottom() || is_top()) {
      // do nothing
    } else {
      CRAB_LOG("tvpi-dbm-join", crab::outs() << "[join]\n"
                                             << *this << "\nwith\n"
                                             << other << "\n");
      counter = counter ? counter : other.counter;
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      m_coeff_map.join(other.m_coeff_map);
#endif
      tvpi_reduce();
      tvpi_dbm_domain_t other2 = other;
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      other2.m_coeff_map.join(m_coeff_map);
#endif
      other2.counter = other2.counter ? other2.counter : counter;
      other2.tvpi_reduce();
      CRAB_LOG("tvpi-dbm-join", crab::outs() << "[reduced join]\n"
                                             << *this << "\nwith\n"
                                             << other2 << "\n");
      m_base_absval |= other2.m_base_absval;
      m_ext_absval |= other2.m_ext_absval;
      filter_dbms();
      prune_coefficients();
      CRAB_LOG("tvpi-dbm-join", crab::outs() << "[res]\n" << *this << "\n");
    }
  }

  tvpi_dbm_domain_t operator|(const tvpi_dbm_domain_t &other) const override {
    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    } else {
      CRAB_LOG("tvpi-dbm-join", crab::outs() << "[join]\n"
                                             << *this << "\nwith\n"
                                             << other << "\n");
      tvpi_dbm_domain_t this2 = *this;
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      this2.m_coeff_map.join(other.m_coeff_map);
#endif
      this2.counter = this2.counter ? this2.counter : other.counter;
      this2.tvpi_reduce();
      tvpi_dbm_domain_t other2 = other;
      other2.counter = other2.counter ? other2.counter : this2.counter;
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      other2.m_coeff_map.join(m_coeff_map);
#endif
      other2.tvpi_reduce();
      CRAB_LOG("tvpi-dbm-join", crab::outs() << "[reduced join]\n"
                                             << this2 << "\nwith\n"
                                             << other2 << "\n");
      this2.m_base_absval |= other2.m_base_absval;
      this2.m_ext_absval |= other2.m_ext_absval;
      this2.filter_dbms();
      this2.prune_coefficients();
      CRAB_LOG("tvpi-dbm-join", crab::outs() << "[res]\n" << this2 << "\n");
      return this2;
    }
  }

  void operator&=(const tvpi_dbm_domain_t &other) override {
    if (is_bottom() || other.is_top()) {
      // do nothing
    } else if (other.is_bottom() || is_top()) {
      *this = other;
    } else {
      CRAB_LOG("tvpi-dbm-meet", crab::outs() << "[meet]\n"
                                             << *this << "\nwith\n"
                                             << other << "\n");
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      m_coeff_map.meet(other.m_coeff_map);
#endif
      counter = counter ? counter : other.counter;
      tvpi_reduce();
      tvpi_dbm_domain_t other2 = other;
      other2.counter = other2.counter ? other2.counter : counter;
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      other2.m_coeff_map.meet(m_coeff_map);
#endif
      other2.tvpi_reduce();
      CRAB_LOG("tvpi-dbm-meet", crab::outs() << "[reduced meet]\n"
                                             << *this << "\nwith\n"
                                             << other2 << "\n");
      m_base_absval &= other2.m_base_absval;
      m_ext_absval &= other2.m_ext_absval;
      if (m_base_absval.is_bottom() || m_ext_absval.is_bottom()) {
        set_to_bottom();
      }
      CRAB_LOG("tvpi-dbm-meet", crab::outs() << "[res]\n" << *this << "\n");
    }
  }

  tvpi_dbm_domain_t operator&(const tvpi_dbm_domain_t &other) const override {
    if (is_bottom() || other.is_top()) {
      return *this;
    } else if (other.is_bottom() || is_top()) {
      return other;
    } else {
      CRAB_LOG("tvpi-dbm-meet", crab::outs() << "[meet]\n"
                                             << *this << "\nwith\n"
                                             << other << "\n");
      tvpi_dbm_domain_t this2 = *this;
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      this2.m_coeff_map.meet(other.m_coeff_map);
#endif
      this2.counter = this2.counter ? this2.counter : other.counter;
      this2.tvpi_reduce();
      tvpi_dbm_domain_t other2 = other;
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      other2.m_coeff_map.meet(m_coeff_map);
#endif
      other2.counter = other2.counter ? other2.counter : this2.counter;
      other2.tvpi_reduce();
      CRAB_LOG("tvpi-dbm-meet", crab::outs() << "[reduced meet]\n"
                                             << this2 << "\nwith\n"
                                             << other2 << "\n");
      this2.m_base_absval &= other2.m_base_absval;
      this2.m_ext_absval &= other2.m_ext_absval;
      if (this2.m_base_absval.is_bottom() || this2.m_ext_absval.is_bottom()) {
        this2.set_to_bottom();
      }
      CRAB_LOG("tvpi-dbm-meet", crab::outs() << "[res]\n" << this2 << "\n");
      return this2;
    }
  }

  tvpi_dbm_domain_t operator||(const tvpi_dbm_domain_t &other) const override {
    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    } else {
      CRAB_LOG("tvpi-dbm-widen", crab::outs() << "[widen]\n"
                                              << *this << "\nwith\n"
                                              << other << "\n");
      tvpi_dbm_domain_t this2 = *this;
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      this2.m_coeff_map.join(other.m_coeff_map);
#endif
      this2.counter = this2.counter ? this2.counter : other.counter;
      this2.tvpi_reduce();
      tvpi_dbm_domain_t other2 = other;
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      other2.m_coeff_map.join(m_coeff_map);
#endif
      other2.counter = other2.counter ? other2.counter : this2.counter;
      other2.tvpi_reduce();
      CRAB_LOG("tvpi-dbm-widen", crab::outs() << "[reduced widen]\n"
                                              << this2 << "\nwith\n"
                                              << other2 << "\n");
      base_domain_t out_base_absval =
          this2.m_base_absval || other2.m_base_absval;
      base_domain_t out_ext_absval = this2.m_ext_absval || other2.m_ext_absval;
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      coefficient_map_t out_coeff_map = this2.m_coeff_map;
      tvpi_dbm_domain_t res(std::move(out_base_absval),
                            std::move(out_ext_absval),
                            std::move(out_coeff_map));
#else
      tvpi_dbm_domain_t res(std::move(out_base_absval),
                            std::move(out_ext_absval));
#endif
      res.counter = this2.counter ? this2.counter : other2.counter;
      res.filter_dbms();
      res.prune_coefficients();
      CRAB_LOG("tvpi-dbm-widen", crab::outs() << "[res]\n" << res << "\n");
      return res;
    }
  }

  tvpi_dbm_domain_t
  widening_thresholds(const tvpi_dbm_domain_t &other,
                      const thresholds<number_t> &ts) const override {
    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    } else {
      CRAB_LOG("tvpi-dbm-widen", crab::outs() << "[widen]\n"
                                              << *this << "\nwith\n"
                                              << other << "\n");
      tvpi_dbm_domain_t this2 = *this;
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      this2.m_coeff_map.join(other.m_coeff_map);
#endif
      this2.counter = this2.counter ? this2.counter : other.counter;
      this2.tvpi_reduce();
      tvpi_dbm_domain_t other2 = other;
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      other2.m_coeff_map.join(m_coeff_map);
#endif
      other2.counter = other2.counter ? other2.counter : this2.counter;
      other2.tvpi_reduce();
      CRAB_LOG("tvpi-dbm-widen", crab::outs() << "[reduced widen]\n"
                                              << this2 << "\nwith\n"
                                              << other2 << "\n");
      base_domain_t out_base_absval =
          this2.m_base_absval.widening_thresholds(other2.m_base_absval, ts);
      base_domain_t out_ext_absval =
          this2.m_ext_absval.widening_thresholds(other2.m_ext_absval, ts);
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      coefficient_map_t out_coeff_map = this2.m_coeff_map;
      tvpi_dbm_domain_t res(std::move(out_base_absval),
                            std::move(out_ext_absval),
                            std::move(out_coeff_map));
#else
      tvpi_dbm_domain_t res(std::move(out_base_absval),
                            std::move(out_ext_absval));
#endif
      res.counter = this2.counter ? this2.counter : other2.counter;
      res.filter_dbms();
      res.prune_coefficients();
      CRAB_LOG("tvpi-dbm-widen", crab::outs() << "[res]\n" << res << "\n");
      return res;
    }
  }

  tvpi_dbm_domain_t operator&&(const tvpi_dbm_domain_t &other) const override {
    if (is_bottom() || other.is_top()) {
      return *this;
    } else if (other.is_bottom() || is_top()) {
      return other;
    } else {
      CRAB_LOG("tvpi-dbm-narrow", crab::outs() << "[narrow]\n"
                                               << *this << "\nwith\n"
                                               << other << "\n");
      tvpi_dbm_domain_t this2 = *this;
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      this2.m_coeff_map.meet(other.m_coeff_map);
#endif
      this2.counter = this2.counter ? this2.counter : other.counter;
      this2.tvpi_reduce();
      tvpi_dbm_domain_t other2 = other;
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      other2.m_coeff_map.meet(m_coeff_map);
#endif
      other2.counter = other2.counter ? other2.counter : this2.counter;
      other2.tvpi_reduce();
      CRAB_LOG("tvpi-dbm-narrow", crab::outs() << "[reduced narrow]\n"
                                               << this2 << "\nwith\n"
                                               << other2 << "\n");
      base_domain_t out_base_absval =
          this2.m_base_absval && other2.m_base_absval;
      base_domain_t out_ext_absval = this2.m_ext_absval && other2.m_ext_absval;
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      coefficient_map_t out_coeff_map = this2.m_coeff_map;
      tvpi_dbm_domain_t res(std::move(out_base_absval),
                            std::move(out_ext_absval),
                            std::move(out_coeff_map));
#else
      tvpi_dbm_domain_t res(std::move(out_base_absval),
                            std::move(out_ext_absval));
#endif
      res.counter = this2.counter ? this2.counter : other2.counter;
      CRAB_LOG("tvpi-dbm-narrow", crab::outs() << "[res]\n" << res << "\n");
      return res;
    }
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    CRAB_LOG("tvpi-dbm-+=",
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

        CRAB_LOG("tvpi-dbm-+=", crab::outs()
                                    << "processing original: " << cst << "\n");

        m_base_absval += cst;
        if (m_base_absval.is_bottom()) {
          set_to_bottom();
          break;
        }
        auto ecst = rewrite_linear_constraint(cst);
        if (ecst.equal(cst)) {
          CRAB_LOG("tvpi-dbm-+=", crab::outs()
                                      << "cannot rewrite: " << cst << "\n");
          continue;
        }
        CRAB_LOG("tvpi-dbm-+=", crab::outs()
                                    << "processing rewritten " << ecst << "\n");
        m_ext_absval += ecst;
        if (m_ext_absval.is_bottom()) {
          set_to_bottom();
          break;
        }

      } // end for
      tvpi_reduce();
    }

    CRAB_LOG("tvpi-dbm-+=",
             crab::outs() << "After assume(" << csts << ")=" << *this << "\n");
  }

  bool entails(const linear_constraint_t &cst) const override {
    if (is_bottom()) {
      return true;
    } else if (cst.is_tautology()) {
      return true;
    } else if (cst.is_contradiction()) {
      return false;
    } else if (m_base_absval.entails(cst) || m_ext_absval.entails(cst)) {
      return true;
    } else {
      auto ecst = try_rewrite_linear_constraint(cst);
      if (ecst && !(*ecst).equal(cst)) {
        if (m_ext_absval.entails((*ecst))) {
          return true;
        }
      }

      bool need_reduce = true;

      if (need_reduce) {
        tvpi_dbm_domain_t tmp = *this;
        tmp.tvpi_reduce();
        CRAB_LOG("tvpi-dbm-entails", crab::outs()
                                         << "reduced: " << tmp << "\n");

        if (tmp.m_base_absval.entails(cst) ||
            (ecst && tmp.m_ext_absval.entails(*ecst))) {
          return true;
        }

        if (ecst && !(*ecst).equal(cst)) {
          if (tmp.m_ext_absval.entails(*ecst)) {
            return true;
          }
        }
      }

      return false;
    }
  }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    // x := c1*x1 + c2*x2 +... + k

    // For DBM, if we have complex linear expression, we can only approximate
    // values c1*x1 + c2*x2 +... + k Thus, to gain precision, we handle several
    // special cases: (1) e == cosntant (2) e == b*y +/- c (3) e == b*y +/- c*z
    // +/- d (4) e == general form Not necessary to rewrite: rewrite("x := e") =
    // "x*COEF := e*COEF" since this new constraint is equivalent to original if
    // necessary, the new constraint can be recover by using reduction.
    if (!is_bottom()) {
      CRAB_LOG("tvpi-dbm-assign", crab::outs()
                                      << "Before assign(" << x << " := " << e
                                      << ")=" << *this << "\n");
      CRAB_LOG("tvpi-dbm-assign", crab::outs() << "processing original " << x
                                               << " := " << e << "\n");
      m_base_absval.assign(x, e);
      rewrite_assign(x, e, false /*weak*/);

      CRAB_LOG("tvpi-dbm-assign", crab::outs() << "After assign(" << x << " := "
                                               << e << ")=" << *this << "\n");
    }
  }

  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    if (!is_bottom()) {
      m_base_absval.weak_assign(x, e);
      rewrite_assign(x, e, true /*weak*/);
    }
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t z) override {
    if (!is_bottom()) {
      m_base_absval.apply(op, x, y, z);
      bool neg = (z < 0);
      number_t sign_one = neg ? number_t(-1) : number_t(1);
      number_t z_abs = neg ? -z : z;

      // rewrite("x := y op z")
      switch (op) {
      case OP_ADDITION:
      case OP_SUBTRACTION:
        if (x == y && counter && *counter != x) { // x := x +/- z
// hint: this might be an index or a counter, add it to extended
// value
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
          m_coeff_map.insert({*counter, convert(z)});
#endif
        }
        break;
      case OP_MULTIPLICATION: // x := y * z
        if (z_abs > number_t(1)) {
// "x := zy for z > 1"
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
          m_ext_absval.apply(op, x, get_ghost_var(y, z_abs), sign_one);
#else
          if (find_ghost_var(y, z_abs)) {
            m_ext_absval.apply(op, x, get_ghost_var(y, z_abs), sign_one);
          }
#endif
        }
        break;
      case OP_SDIV:                // x := y /s z
      case OP_UDIV:                // x := y /u z
        if (z_abs > number_t(1)) { // "zx := y for z > 1"
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
          m_ext_absval.apply(op, get_ghost_var(x, z_abs), y, sign_one);
#else
          if (find_ghost_var(x, z_abs)) {
            m_ext_absval.apply(op, get_ghost_var(x, z_abs), y, sign_one);
          }
#endif
        }
        break;
      case OP_SREM: // x := y %/s z
      case OP_UREM: // x := y %/u z
        // This can be rewritten into b - c * floor(b / c) if
        // b / c in low level domain is equivalent to floor(b / c).
        break;
      default:
        break;
      }

      // add x := x op z if x in extended dbm
      if (x == y) {
        auto extvars = m_ext_absval.vars();
        if (tvpi_utils::find(extvars, x)) {
          m_ext_absval.apply(op, x, y, z);
        }
      }

#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      auto it = m_coeff_map.find(x);
      if (it == m_coeff_map.end()) {
        return;
      }
      auto &coeffs = it->second;
#else
      auto &coeffs = crab_domain_params_man::get().coefficients();
#endif
      for (auto coefficient : coeffs) {
        rewrite_apply(op, x, y, z, coefficient);
      }
    }
  }

  void eval_apply(arith_operation_t op, const variable_t &x,
                  const variable_t &y, const variable_t &z) {

    // check if y or z are constants
    auto itv_y = m_base_absval[y];
    if (auto c = itv_y.singleton()) { // x := y op z if eval(y) = c
      switch (op) {                   // x := c op z
      case OP_ADDITION:
      case OP_MULTIPLICATION:
        // x := z op c since these operations are commutative
        apply(op, x, z, *c);
        break;
      case OP_SUBTRACTION: // x := c - z => need an additional ghost variable
        break;
      case OP_SDIV: // x := c /s z
      case OP_UDIV: // x := c /u z
        break;
      default:
        break;
      }
    }
    auto itv_z = m_base_absval[z];
    if (auto c = itv_z.singleton()) {
      apply(op, x, y, *c);
    }
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    if (!is_bottom()) {
      m_base_absval.apply(op, x, y, z);
      eval_apply(op, x, y, z);
      eval_apply(op, x, z, y);
      if (x == y || x == z) {
        auto extvars = m_ext_absval.vars();
        if (tvpi_utils::find(extvars, x)) {
          m_ext_absval.apply(op, x, y, z);
        }
      }
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      auto it = m_coeff_map.find(x);
      if (it == m_coeff_map.end()) {
        return;
      }
      auto &coeffs = it->second;
#else
      auto &coeffs = crab_domain_params_man::get().coefficients();
#endif
      for (auto coefficient : coeffs) {
        rewrite_apply(op, x, y, z, coefficient);
      }
    }
  }

  // dst := trunc or extend src
  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    if (!is_bottom()) {
      m_base_absval.apply(op, dst, src);
      if (dst == src) {
        auto extvars = m_ext_absval.vars();
        if (tvpi_utils::find(extvars, dst)) {
          m_ext_absval.apply(op, dst, src);
        }
      }
    }
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    if (!is_bottom()) {
      m_base_absval.apply(op, x, y, z);
      if (x == y || x == z) {
        auto extvars = m_ext_absval.vars();
        if (tvpi_utils::find(extvars, x)) {
          m_ext_absval.apply(op, x, y, z);
        }
      }
    }
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t z) override {
    if (!is_bottom()) {
      m_base_absval.apply(op, x, y, z);
      if (x == y) {
        auto extvars = m_ext_absval.vars();
        if (tvpi_utils::find(extvars, x)) {
          m_ext_absval.apply(op, x, y, z);
        }
      }
    }
  }

  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const tvpi_dbm_domain_t &invariant) override {
    CRAB_WARN(domain_name(), "::backward_assign not implemented");
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const tvpi_dbm_domain_t &invariant) override {
    CRAB_WARN(domain_name(), "::backward_apply not implemented");
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const tvpi_dbm_domain_t &invariant) override {
    CRAB_WARN(domain_name(), "::backward_apply not implemented");
  }

  void callee_entry(const callsite_info<variable_t> &callsite,
                    const tvpi_dbm_domain_t &caller) override {
    inter_abstract_operations<
        tvpi_dbm_domain_t,
        Params::implement_inter_transformers>::callee_entry(callsite, caller,
                                                            *this);
  }

  void caller_continuation(const callsite_info<variable_t> &callsite,
                           const tvpi_dbm_domain_t &callee) override {
    inter_abstract_operations<
        tvpi_dbm_domain_t,
        Params::implement_inter_transformers>::caller_continuation(callsite,
                                                                   callee,
                                                                   *this);
  }

  linear_constraint_system_t to_linear_constraint_system() const override {

    if (is_bottom()) {
      return linear_constraint_system_t(linear_constraint_t::get_false());
    }

    if (is_top()) {
      return linear_constraint_system_t(linear_constraint_t::get_true());
    }

    // TODO: convet ghost variables to real variables
    return m_base_absval.to_linear_constraint_system();
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    CRAB_WARN(domain_name(),
              "::to_disjunctive_linear_constraint_system not implemented");
    disjunctive_linear_constraint_system_t res;
    return res;
  }

  void operator-=(const variable_t &var) override {
    if (!(is_bottom() || is_top())) {
      CRAB_LOG("tvpi-dbm-forget", crab::outs() << "Before forget[" << var
                                               << "]=" << *this << "\n");
      tvpi_reduce();
      m_base_absval -= var;
      m_ext_absval -= var;

#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      auto it = m_coeff_map.find(var);
      if (it == m_coeff_map.end()) {
        return;
      }
      auto &coeffs = it->second;
#else
      auto &coeffs = crab_domain_params_man::get().coefficients();
#endif
      for (auto coefficient : coeffs) {
        variable_t ghost_var = get_ghost_var(var, coefficient);
        m_ext_absval -= ghost_var;
      }
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      m_coeff_map.remove(var);
#endif
      CRAB_LOG("tvpi-dbm-forget", crab::outs()
                                      << "After forget=" << *this << "\n");
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

  void forget(const variable_vector_t &variables) override {
    if (!(is_bottom() || is_top())) {
      CRAB_LOG("tvpi-dbm-forget", crab::outs() << "Before forget";
               tvpi_utils::print_vector(crab::outs(), variables);
               crab::outs() << "=" << *this << "\n");
      tvpi_reduce();
      m_base_absval.forget(variables);
      variable_vector_t allvars(variables);
      for (auto const &v : variables) {
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
        auto it = m_coeff_map.find(v);
        if (it != m_coeff_map.end()) {
          for (auto coefficient : it->second) {
            variable_t gv = get_ghost_var(v, coefficient);
            allvars.push_back(gv);
          }
          m_coeff_map.remove(v);
        }
#else
        auto &coeffs = crab_domain_params_man::get().coefficients();
        for (auto coefficient : coeffs) {
          variable_t gv = get_ghost_var(v, coefficient);
          allvars.push_back(gv);
        }
#endif
      }
      m_ext_absval.forget(allvars);
      CRAB_LOG("tvpi-dbm-forget", crab::outs()
                                      << "After forget=" << *this << "\n");
    }
  }

  void project(const variable_vector_t &variables) override {
    if (!is_bottom()) {
      m_base_absval.project(variables);
      variable_vector_t allvars(variables);
      for (auto const &v : variables) {
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
        auto it = m_coeff_map.find(v);
        if (it != m_coeff_map.end()) {
          for (auto coefficient : it->second) {
            variable_t gv = get_ghost_var(v, coefficient);
            allvars.push_back(gv);
          }
        }
#else
        auto &coeffs = crab_domain_params_man::get().coefficients();
        for (auto coefficient : coeffs) {
          variable_t gv = get_ghost_var(v, coefficient);
          allvars.push_back(gv);
        }
#endif
      }
      m_ext_absval.project(allvars);
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
      m_coeff_map.keep(variables);
#endif
    }
  }

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    if (!is_bottom()) {
      m_base_absval.rename(from, to);
      variable_vector_t extd_from(from);
      variable_vector_t extd_to(to);
      for (unsigned i = 0, sz = from.size(); i < sz; ++i) {
        const variable_t &f = from[i];
        const variable_t &t = to[i];
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
        auto it = m_coeff_map.find(f);
        if (it != m_coeff_map.end()) {
          for (auto coefficient : it->second) {
            extd_from.push_back(get_ghost_var(f, coefficient));
            extd_to.push_back(get_ghost_var(t, coefficient));
          }
          auto set2 = it->second;
          m_coeff_map.insert({t, std::move(set2)});
          m_coeff_map.remove(f);
        }
#else
        auto &coeffs = crab_domain_params_man::get().coefficients();
        for (auto coefficient : coeffs) {
          extd_from.push_back(get_ghost_var(f, coefficient));
          extd_to.push_back(get_ghost_var(t, coefficient));
        }
#endif
      }
      m_ext_absval.rename(extd_from, extd_to);
    }
  }

  void expand(const variable_t &var, const variable_t &new_var) override {
    if (is_bottom() || is_top()) {
      return;
    }

    m_base_absval.expand(var, new_var);
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
    auto it = m_coeff_map.find(var);
    if (it != m_coeff_map.end()) {
      for (auto coefficient : it->second) {
        variable_t gv = get_ghost_var(var, coefficient);
        variable_t gnv = get_ghost_var(new_var, coefficient);
        m_ext_absval.expand(gv, gnv);
      }
      auto set2 = it->second;
      m_coeff_map.insert({new_var, std::move(set2)});
    }
#else
    auto &coeffs = crab_domain_params_man::get().coefficients();
    for (auto coefficient : coeffs) {
      variable_t gv = get_ghost_var(var, coefficient);
      variable_t gnv = get_ghost_var(new_var, coefficient);
      m_ext_absval.expand(gv, gnv);
    }
#endif
  }

  void normalize() override {}
  void minimize() override {}

  void intrinsic(std::string name, const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    auto error_if_not_variable =
        [&, func = __func__](const variable_or_constant_t &vc) {
          if (!vc.is_variable()) {
            CRAB_ERROR(domain_name(), "::", func, " ", name,
                       " expected a variable input");
          }
        };
    if (name == "loop_counter") {
      assert(inputs.size() == 1);
      error_if_not_variable(inputs[0]);
      variable_t i = inputs[0].get_variable();
      counter = i;
    }
  }

  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const tvpi_dbm_domain_t &invariant) override {
    CRAB_WARN(domain_name(), "::backward_intrinsic for ", name,
              " not implemented");
  }

  void write(crab_os &o) const override {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "top";
    } else {
      o << "{";
      {
        o << "coeffs:";
#if TVPI_DBM_FIXED_COEFFICIENTS == 0
        m_coeff_map.write(o);
#else
        tvpi_utils::print_vector(o,
                                 crab_domain_params_man::get().coefficients());
#endif
        o << ", ";
      }
      o << "base:";
      m_base_absval.write(o);
      o << ", extend:";
      m_ext_absval.write(o);
      o << "}";
    }
  }

  friend crab_os &operator<<(crab_os &o, const tvpi_dbm_domain_t &val) {
    val.write(o);
    return o;
  }

  std::string domain_name() const override {
    base_domain_t absval;
    std::string base_name = absval.domain_name();
    const char *prefix = "TVPI";
    std::string name;
    name.reserve(base_name.size() + 7);
    name.append(prefix);
    name.append("(");
    name.append(base_name);
    name.append(")");
    return name;
  }
};

template <typename OctLikeDomain, typename Params>
struct abstract_domain_traits<tvpi_dbm_domain<OctLikeDomain, Params>> {
  using number_t = typename OctLikeDomain::number_t;
  using varname_t = typename OctLikeDomain::varname_t;
};

} // end namespace domains
} // end namespace crab
