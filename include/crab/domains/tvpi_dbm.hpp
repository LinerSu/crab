#pragma once

#include <algorithm>
#include <boost/optional.hpp>
#include <chrono>
#include <string>

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/inter_abstract_operations.hpp>
#include <crab/numbers/bignums.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>

namespace crab {
namespace domains {

#define tvpi_dbm_domain_SCOPED_STATS(NAME)                                     \
  CRAB_DOMAIN_SCOPED_STATS(this, NAME, 1)
#define tvpi_dbm_domain_SCOPED_STATS_ASSIGN_CTOR(NAME)                         \
  CRAB_DOMAIN_SCOPED_STATS(&o, NAME, 0)

namespace tvpi_imple {
/// @brief a special log method to print vector
/// @tparam TType
/// @param o crab ostream
/// @param vec the vector for printing
template <typename TType>
void print_vector(crab::crab_os &o, const std::vector<TType> &vec) {
  typename std::vector<TType>::const_iterator it;
  o << "[";
  for (it = vec.begin(); it != vec.end(); it++) {
    if (it != vec.begin())
      o << ",";
    o << (*it);
  }
  o << "]";
}

// based on Farkas Lemma, we can remove some constraints:
// given two constraints a1*x1 + b1*x2 <= c1 and a2*x1 + b2*x2 <= c2
// check is the conjunction of the two constraints implies the third one
// a3*x1 + b3*x2 <= c3
// if so, then we can safely remove the third constraint.
// otherwise, we encounter unsat result.
// Based on Farkas Lemma, let A be:
//  | a1 a2 |
//  | b1 b2 |
/// if we found lambda1 and lambda2 such that:
//  a1*lambda1 + a2*lambda2 = a3
//  b1*lambda1 + b2*lambda2 = b3
//   => Determinants for
//   | a1 a2 | | a3 a2 | | a1 a3 |
//   | b1 b2 | | b3 b2 | | b1 b3 |
//   det = a1*b2 - a2*b1, det' = a3*b2 - a2*b3, det'' = a1*b3 - b1*a3
//   Cramer's rule:
//    lambda1 = det'  / det
//    lambda2 = det'' / det
//   check (1) lambda1 < 0 || lambda2 < 0
//   if not (1), check (2) c3 >= lambda1*c1 + lambda2*c2
//   if (2) is true, then the third constraint is entailed by the first two
template <typename number>
bool check_entailement(const std::tuple<int, int, number> &cst1,
                       const std::tuple<int, int, number> &cst2,
                       const std::tuple<int, int, number> &cst3) {
  const int &a1 = std::get<0>(cst1), &b1 = std::get<1>(cst1);
  const int &a2 = std::get<0>(cst2), &b2 = std::get<1>(cst2);
  const int &a3 = std::get<0>(cst3), &b3 = std::get<1>(cst3);

  const int64_t &c1 = static_cast<int64_t>(std::get<2>(cst1));
  const int64_t &c2 = static_cast<int64_t>(std::get<2>(cst2));
  const int64_t &c3 = static_cast<int64_t>(std::get<2>(cst3));

  int det = a1 * b2 - a2 * b1;
  if (det == 0) {
    return false;
  }
  int det1 = a3 * b2 - a2 * b3;
  int det2 = a1 * b3 - b1 * a3;

  // Cramer's rule
  double lambda1 = static_cast<double>(det1) / det;
  double lambda2 = static_cast<double>(det2) / det;

  // Check if lambda1 and lambda2 are non-negative
  if (lambda1 < 0 || lambda2 < 0) {
    return false;
  }

  // Check if c3 >= lambda1 * c1 + lambda2 * c2
  return static_cast<double>(c3) >= lambda1 * c1 + lambda2 * c2;
}

template <typename key, typename value> class coefficient_map {
  using key_t = key;
  using value_t = value;
  using key_value_t = std::pair<const key_t, value_t>;
  using value_set_t = std::set<value_t>;
  using key_value_set_t = std::pair<const key_t, value_set_t>;
  using map_t = std::unordered_map<key_t, value_set_t>;
  using const_iterator_t = typename map_t::const_iterator;
  using iterator_t = typename map_t::iterator;
  using shared_map_t = std::shared_ptr<map_t>;

public:
  using coefficient_map_t = coefficient_map<key, value>;
  using coefficient_set_t = value_set_t;
  // initializes but no map created
  coefficient_map() : m_map(nullptr) {}
  coefficient_map(const coefficient_map_t &o) = default;
  coefficient_map(coefficient_map_t &&o) = default;
  coefficient_map_t &operator=(const coefficient_map_t &o) = default;
  coefficient_map_t &operator=(coefficient_map_t &&o) = default;

  void insert(const key_value_t &kv) {
    ensure_exists();
    auto &m = *m_map;
    auto it = m.find(kv.first);
    if (it != m.end()) {
      it->second.insert(kv.second);
    } else {
      m.insert({kv.first, {kv.second}});
    }
  }

  void insert(key_value_t &&kv) {
    ensure_exists();
    auto &m = *m_map;
    auto it = m.find(kv.first);
    if (it != m.end()) {
      it->second.insert(std::forward<key_value_t>(kv).second);
    } else {
      m.insert({std::forward<key_value_t>(kv).first,
                {std::forward<key_value_t>(kv).second}});
    }
  }

  void insert(const key_value_set_t &kvs) {
    ensure_exists();
    auto &m = *m_map;
    auto it = m.find(kvs.first);
    if (it != m.end()) {
      it->second.insert(kvs.second.begin(), kvs.second.end());
    } else {
      m.insert(kvs);
    }
  }

  void insert(key_value_set_t &&kvs) {
    ensure_exists();
    auto &m = *m_map;
    auto it = m.find(kvs.first);
    if (it != m.end()) {
      it->second.insert(std::make_move_iterator(kvs.second.begin()),
                        std::make_move_iterator(kvs.second.end()));
    } else {
      m.insert(std::move(kvs));
    }
  }

  const_iterator_t find(const key_t &k) const {
    ensure_exists();
    return m_map->find(k);
  }

  iterator_t find(const key_t &k) {
    ensure_exists();
    return m_map->find(k);
  }

  const_iterator_t begin() const {
    ensure_exists();
    return m_map->begin();
  }

  iterator_t begin() {
    ensure_exists();
    return m_map->begin();
  }

  const_iterator_t end() const {
    ensure_exists();
    return m_map->end();
  }

  iterator_t end() {
    ensure_exists();
    return m_map->end();
  }

  iterator_t erase(iterator_t it) {
    ensure_exists();
    return m_map->erase(it);
  }

  void merge(const coefficient_map_t &o) {
    if (!exists()) {
      m_map = o.m_map;
    }
    if (!o.exists()) {
      return;
    }

    for (auto it = o.begin(), et = o.end(); it != et; ++it) {
      const auto &k = it->first;
      auto it2 = find(k);
      if (it2 == end()) {
        m_map->insert({k, it->second});
      } else {
        it2->second.insert(it->second.begin(), it->second.end());
      }
    }
  }

  void write(crab_os &o) const {
    if (!exists()) {
      o << "not exists";
      return;
    } else {
      auto &m = *m_map;
      o << "{";
      for (auto it = m.begin(), et = m.end(); it != et;) {
        o << it->first << " => (";
        for (auto sit = it->second.begin(), set = it->second.end();
             sit != set;) {
          o << *sit;
          ++sit;
          if (sit != set) {
            o << ", ";
          }
        }
        o << ")";
        ++it;
        if (it != et) {
          o << ", ";
        }
      }
      o << "}";
    }
  }

  friend crab_os &operator<<(crab_os &o, const coefficient_map_t &m) {
    m.write(o);
    return o;
  }

  template <typename Key, typename Value>
  void print_unordered_map(crab::crab_os &o,
                           std::unordered_map<Key, Value> const &m) {
    o << "{";
    for (auto it = m.begin(), et = m.end(); it != et;) {
      o << it->first << " => " << it->second;
      ++it;
      if (it != et) {
        o << ", ";
      }
    }
    o << "}";
  }

private:
  shared_map_t m_map;
  void ensure_exists() {
    if (!m_map) {
      m_map = std::make_shared<map_t>();
    }
  }
  void ensure_exists() const {
    if (!m_map) {
      throw std::runtime_error("shared map does not exist!");
    }
  }
  bool exists() const { return m_map != nullptr; }
};

} // namespace tvpi_imple

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
  using coefficient_map_t =
      typename tvpi_imple::coefficient_map<variable_t,
                                           unsigned>::coefficient_map_t;
  using coefficient_set_t = typename coefficient_map_t::coefficient_set_t;

  static_assert(std::is_same<typename abstract_domain_api_t::number_t,
                             ikos::z_number>::value,
                "abstract_domain_api_t::number_t must be the ikos::z_number");
  // design, a dbm domain with original variables
  // an extended dbm domain with coefficient ghost variables
  // original variables may be shared with both domains
  // the question is how to handle reduction between originals shared in two
  // domains the trival solution is meet the two domains and then project to the
  // scope of domains.

private:
  using base_domain_t = OctLikeDomain;
  using coeff_map_t = std::unordered_map<variable_t, std::set<unsigned>>;

  base_domain_t m_base_absval;
  base_domain_t m_ext_absval;
  coefficient_map_t m_coeff_map;

  variable_t get_ghost_var(const variable_t &v, unsigned coefficient) {
    if (coefficient == 0) {
      CRAB_ERROR("Coefficient must be > 0");
    } else if (coefficient == 1) {
      return v;
    }

    auto it = m_coeff_map.find(v);
    if (it != m_coeff_map.end()) {
      it->second.insert(coefficient);
    } else {
      m_coeff_map.insert({v, coefficient});
    }
    auto &vfac = const_cast<varname_t *>(&(v.name()))->get_var_factory();

    variable_t coeff_v(vfac.get_or_insert_varname(v.name(), coefficient),
                       v.get_type());
    return coeff_v;
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
  }

  boost::optional<variable_t>
  find_ghost_var(const variable_t &v, const number_t &coefficient) const {
    return find_ghost_var(v, convert(coefficient));
  }

  linear_expression_t rewrite_linear_expression(const linear_expression_t &e) {
    /**
     *
     * Given c1*x1 + c2*x2 +... + k, rewrite each ci*xi into
     *
     * cixi         if ci is one of tracked coefficients then
     *    c1x1 + c2x2 +... + k
     *      ci*xi  otherwise
     **/
    linear_expression_t res;
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const variable_t &v = (*it).second;
      const number_t &coeff = (*it).first;
      if (coeff == 0) {
        continue;
      } else if (coeff > 0) {
        res = res + get_ghost_var(v, coeff);
      } else if (coeff < 0) {
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
    return linear_constraint_t(rewrite_linear_expression(cst.expression()),
                               cst.kind());
  }

  linear_expression_t rewrite_linear_expression(const linear_expression_t &e,
                                                unsigned coefficient,
                                                bool divd) const {
    /**
     *
     * Given c1*x1 + c2*x2 +... + k, and a coefficient c. rewrite into
     *       (c op c1)x1 + (c op c2)x2 +... + c op k
     * where op is either division or multiplication
     **/

    number_t tracked_coeff(coefficient);
    linear_expression_t res;
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const variable_t &v = (*it).second;
      const number_t &coeff = (*it).first;
      bool neg = coeff < 0;
      if (coeff == 0) {
        continue;
      }
      if (divd) {
        if (coeff % coefficient == 0) {
          number_t new_coeff =
              neg ? -coeff / tracked_coeff : coeff / tracked_coeff;
          if (auto vgvar = find_ghost_var(v, new_coeff)) {
            res = neg ? res - *vgvar : res + *vgvar;
          } else {
            res = neg ? res - new_coeff * v : res + new_coeff * v;
          }
        } else { // give up, cannot rewrite it
          return e;
        }
      } else { // multiple
        number_t new_coeff =
            neg ? -coeff * tracked_coeff : coeff * tracked_coeff;
        if (auto vgvar = find_ghost_var(v, new_coeff)) {
          res = neg ? res - *vgvar : res + *vgvar;
        } else {
          res = neg ? res - new_coeff * v : res + new_coeff * v;
        }
      }
    }
    if (divd) {
      if (e.constant() % tracked_coeff == 0) {
        res = res + e.constant() / tracked_coeff;
      } else { // give up, cannot rewrite it
        return e;
      }
    } else {
      res = res + e.constant() * tracked_coeff;
    }
    return res;
  }

  linear_constraint_t rewrite_linear_constraint(const linear_constraint_t &cst,
                                                unsigned coefficient,
                                                bool divd) const {
    return linear_constraint_t(
        rewrite_linear_expression(cst.expression(), coefficient, divd),
        cst.kind());
  }

  boost::optional<linear_expression_t>
  try_rewrite_linear_expression(const linear_expression_t &e) const {
    linear_expression_t res;
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const variable_t &v = (*it).second;
      const number_t &coeff = (*it).first;
      if (coeff == 0) {
        continue;
      } else if (coeff > 0) {
        auto gv = find_ghost_var(v, coeff);
        if (gv == boost::none) {
          return boost::none;
        }
        res = res + *gv;
      } else if (coeff < 0) {
        auto gv = find_ghost_var(v, -coeff);
        if (gv == boost::none) {
          return boost::none;
        }
        res = res - *gv;
      } else { // give up
        return boost::none;
      }
    }
    res = res + e.constant();
    return res;
  }

  boost::optional<linear_constraint_t>
  try_rewrite_linear_constraint(const linear_constraint_t &cst) const {
    auto e_opt = try_rewrite_linear_expression(cst.expression());
    if (e_opt == boost::none) {
      return boost::none;
    }
    return linear_constraint_t(*e_opt, cst.kind());
  }

  void rewrite_assign(const variable_t &x, const linear_expression_t &e,
                      unsigned coefficient, bool weak) {
    assert(coefficient > 1);

    variable_t ghost_x = get_ghost_var(x, coefficient);
    number_t tracked_coefficient(coefficient);
    if (e.is_constant()) {
      // rewrite("x := n") = "x * COEF := n * COEF"
      if (!weak) {
        m_ext_absval.assign(ghost_x, e * tracked_coefficient);
      } else {
        m_ext_absval.weak_assign(ghost_x, e * tracked_coefficient);
      }
    } else if (boost::optional<variable_t> y = e.get_variable()) {
      // rewrite("x := y") = "x * COEF := y * COEF"
      variable_t ghost_y = get_ghost_var(*y, coefficient);
      if (!weak) {
        m_ext_absval.assign(ghost_x, ghost_y);
      } else {
        m_ext_absval.weak_assign(ghost_x, ghost_y);
      }
    } else {
      if (!weak) {
        linear_expression_t e1 =
            rewrite_linear_expression(e, coefficient, true);
        if (!e1.equal(e)) {
          m_ext_absval.assign(ghost_x, e1);
        }
        e1 = rewrite_linear_expression(e, coefficient, false);
        if (!e1.equal(e)) {
          m_ext_absval.assign(ghost_x, e1);
        }
      } else {
        linear_expression_t e1 =
            rewrite_linear_expression(e, coefficient, true);
        if (!e1.equal(e)) {
          m_ext_absval.weak_assign(ghost_x, e1);
        }
        e1 = rewrite_linear_expression(e, coefficient, false);
        if (!e1.equal(e)) {
          m_ext_absval.weak_assign(ghost_x, e1);
        }
      }
    }
  }

  void rewrite_apply(arith_operation_t op, const variable_t &x,
                     const variable_t &y, number_t z, unsigned coefficient) {
    assert(coefficient > 1);

    variable_t ghost_x = get_ghost_var(x, coefficient);
    variable_t ghost_y = get_ghost_var(y, coefficient);
    number_t tracked_coefficient(coefficient);
    switch (op) {
    case OP_MULTIPLICATION:
    case OP_SDIV:
    case OP_UDIV: // x := y * z or x := y / z
      // rewrite to x*COEF := y*COEF op z
      m_ext_absval.apply(op, ghost_x, ghost_y, z);
      break;
    case OP_ADDITION:
    case OP_SUBTRACTION: // x := y + z or x := y - z
      // rewrite to x*COEF := y*COEF op z*COEF
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
    // rewrite("x := y op z") = "x*COEF := y*COEF op z*COEF"
    variable_t ghost_x = get_ghost_var(x, coefficient);
    variable_t ghost_y = get_ghost_var(y, coefficient);
    variable_t ghost_z = get_ghost_var(z, coefficient);
    m_ext_absval.apply(op, ghost_x, ghost_y, ghost_z);
  }

  tvpi_dbm_domain(base_domain_t &&base, base_domain_t &&extd,
                  coefficient_map_t &&coeff_map, bool apply_reduction = true)
      : m_base_absval(std::move(base)), m_ext_absval(std::move(extd)),
        m_coeff_map(std::move(coeff_map)) {
    // reduce();
  }

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
    m_base_absval.set_to_top();
    m_ext_absval.set_to_top();
  }

  void set_to_bottom() override {
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
    return m_base_absval.is_bottom() || m_ext_absval.is_bottom();
  }

  bool is_top() const override {
    return m_base_absval.is_top() && m_ext_absval.is_top();
  }

  void scaling() {
    auto add_range = [](unsigned ratio, base_domain_t &absval,
                        const variable_t &var, interval_t range, bool div) {
      auto r = interval_t(number_t(ratio));
      interval_t new_range = div ? range / r : range * r;
      if (auto mx = new_range.ub().number()) {
        absval += (var <= *mx);
      }
      if (auto mn = new_range.lb().number()) {
        absval += (var >= *mn);
      }
    };
    auto handle_scaling_down = [](unsigned ratio, base_domain_t &absval,
                                  const variable_t &var1,
                                  const variable_t &var2, number_t c1) {
      boost::optional<number_t> copt2 = absval.difference_bound(var2, var1);
      copt2 = copt2 ? std::min(*copt2, (c1 % ratio == 0 ? c1 / ratio : *copt2))
                    : (c1 % ratio == 0 ? c1 / number_t(ratio)
                                       : boost::optional<number_t>());

      if (copt2) {
        auto cons2 = (var1 - var2 <= *copt2);
        CRAB_LOG("fixed-tvpi-scaling",
                 crab::outs() << "Scaling found: " << cons2 << "\n");
        absval += cons2;
      }
    };
    auto handle_scaling_up = [](unsigned ratio, base_domain_t &absval,
                                const variable_t &var1, const variable_t &var2,
                                number_t c1) {
      auto copt2 = absval.difference_bound(var2, var1);
      auto c2 = copt2 ? std::min(*copt2, c1 * ratio) : c1 * ratio;
      auto cons2 = (var1 - var2 <= c2);
      CRAB_LOG("fixed-tvpi-scaling", crab::outs()
                                         << "Scaling found: " << cons2 << "\n");
      absval += cons2;
    };
    // scaling between a1 * x - b1 * y <= c1 and a2 * x - b2 * y <= c2
    // if a1 / a2 == b1 / b2 == k then we can scaling constraints down / up
    // based on k to minimize c2 and c1.
    for (auto it1 = m_coeff_map.begin(); it1 != m_coeff_map.end(); ++it1) {
      const variable_t &x = it1->first;
      const coefficient_set_t &a_coeffs = it1->second;

      // scaling down
      // for lb <= ax <= ub => lb / a <= x <= ub / a
      for (auto ita = a_coeffs.rbegin(); ita != a_coeffs.rend(); ++ita) {
        auto gax = get_ghost_var(x, *ita);
        interval_t rax = m_ext_absval.at(gax);
        add_range(*ita, m_base_absval, x, rax, true);
      }
      for (auto it2 = std::next(it1); it2 != m_coeff_map.end(); ++it2) {
        const variable_t &y = it2->first;
        const coefficient_set_t &b_coeffs = it2->second;
        // for ax - by <= c => a'x - b'y <= min(c', c / k) where a' = a / k and
        // b' = b / k c' may exist if a'x - b'y <= c' exists.
        for (auto ita = a_coeffs.rbegin(); ita != a_coeffs.rend(); ++ita) {
          for (auto itb = b_coeffs.rbegin(); itb != b_coeffs.rend(); ++itb) {
            auto gax = get_ghost_var(x, *ita);
            auto gby = get_ghost_var(y, *itb);
            auto copt = m_ext_absval.difference_bound(gby, gax);
            if (!copt)
              continue;
            for (auto ita2 = std::next(ita); ita2 != a_coeffs.rend(); ++ita2) {
              if (*ita % *ita2 != 0)
                continue;
              unsigned ratio = *ita / *ita2;

              for (auto itb2 = std::next(itb); itb2 != b_coeffs.rend();
                   ++itb2) {
                if (*itb % *itb2 != 0 || (*itb / *itb2) != ratio)
                  continue;
                auto gax2 = get_ghost_var(x, *ita2);
                auto gby2 = get_ghost_var(y, *itb2);
                auto cons1 = (gax - gby <= *copt);
                handle_scaling_down(ratio, m_ext_absval, gax2, gby2, *copt);
              }
            }
            // case for base_absval
            if (*itb == *ita) {
              auto ratio = *ita;
              handle_scaling_down(ratio, m_base_absval, x, y, *copt);
              // add range:
              interval_t rax = m_ext_absval.at(gax);
              interval_t rby = m_ext_absval.at(gby);
              add_range(ratio, m_base_absval, x, rax, true);
              add_range(ratio, m_base_absval, y, rby, true);
            }
          }
        }
        // scaling up
        for (auto ita = a_coeffs.begin(); ita != a_coeffs.end(); ++ita) {
          auto gax = get_ghost_var(x, *ita);
          for (auto itb = b_coeffs.begin(); itb != b_coeffs.end(); ++itb) {
            auto gby = get_ghost_var(y, *itb);
            if (*itb == *ita) { // case for base_absval
              auto copt = m_base_absval.difference_bound(y, x);
              if (copt) {
                handle_scaling_up(*ita, m_ext_absval, gax, gby, *copt);
                interval_t rx = m_base_absval.at(x);
                interval_t ry = m_base_absval.at(y);
                add_range(*ita, m_ext_absval, gax, rx, false);
                add_range(*ita, m_ext_absval, gby, ry, false);
              }
            }
            auto copt = m_ext_absval.difference_bound(gby, gax);
            if (!copt)
              continue;

            for (auto ita2 = std::next(ita); ita2 != a_coeffs.end(); ++ita2) {
              if (*ita2 % *ita != 0)
                continue;
              unsigned ratio = *ita2 / *ita;

              for (auto itb2 = std::next(itb); itb2 != b_coeffs.end(); ++itb2) {
                if (*itb2 % *itb != 0 || *itb2 / *itb != ratio)
                  continue;
                auto gax2 = get_ghost_var(x, *ita2);
                auto gby2 = get_ghost_var(y, *itb2);
                handle_scaling_up(ratio, m_ext_absval, gax2, gby2, *copt);
              }
            }
          }
        }
      }
      for (auto ita = a_coeffs.rbegin(); ita != a_coeffs.rend(); ++ita) {
        // for lb <= x <= ub => lb * a <= ax <= ub * a
        // add range:
        auto gax = get_ghost_var(x, *ita);
        interval_t rx = m_base_absval.at(x);
        add_range(*ita, m_ext_absval, gax, rx, false);
      }
    }
  }

  void adding() {
    // adding between a1 * x - b1 * y <= c1 and a2 * x - b2 * y <= c2
    // if a1 + a2 (or b1 + b2) exists in coefficient set.
    for (auto it1 = m_coeff_map.begin(); it1 != m_coeff_map.end(); ++it1) {
      for (auto it2 = std::next(it1); it2 != m_coeff_map.end(); ++it2) {
        const variable_t &x = it1->first;
        const coefficient_set_t &a_coeffs = it1->second;
        const variable_t &y = it2->first;
        const coefficient_set_t &b_coeffs = it2->second;
        for (auto ita = a_coeffs.rbegin(); ita != a_coeffs.rend(); ++ita) {
          for (auto itb = b_coeffs.rbegin(); itb != b_coeffs.rend(); ++itb) {
            auto gax = get_ghost_var(x, *ita);
            auto gby = get_ghost_var(y, *itb);
            if (auto copt = m_ext_absval.difference_bound(gby, gax)) {
              for (auto ita2 = std::next(ita); ita2 != a_coeffs.rend();
                   ++ita2) {
                auto new_a = *ita + *ita2;
                if (a_coeffs.find(new_a) == a_coeffs.end()) {
                  continue;
                }
                for (auto itb2 = std::next(itb); itb2 != b_coeffs.rend();
                     ++itb2) {
                  auto new_b = *itb + *itb2;
                  if (b_coeffs.find(new_b) == b_coeffs.end()) {
                    continue;
                  }
                  auto gax2 = get_ghost_var(x, *ita2);
                  auto gby2 = get_ghost_var(y, *itb2);
                  if (auto copt2 = m_ext_absval.difference_bound(gby2, gax2)) {
                    auto gax3 = get_ghost_var(x, new_a);
                    auto gby3 = get_ghost_var(y, new_b);
                    auto cons = (gax3 - gby3 <= *copt + *copt2);
                    CRAB_LOG("fixed-tvpi-adding",
                             crab::outs() << "Adding found: " << cons << "\n");
                    m_ext_absval += cons;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  void simplifying() {
    // Check if any constraint in the form a1 * x - a2 * x <= c1 exists.
    // It can be simplified to (a1 - a2) * x <= c1 as an interval.
    // However, we do not remove the original constraint, because it may be
    // an implicit result of the closure algorithm.
    for (auto it1 = m_coeff_map.begin(); it1 != m_coeff_map.end(); ++it1) {
      const variable_t &x = it1->first;
      const coefficient_set_t &a_coeffs = it1->second;
      coefficient_set_t a_coeffs_with_one = a_coeffs;
      a_coeffs_with_one.insert(1);
      for (auto ita = a_coeffs.rbegin(); ita != a_coeffs.rend(); ++ita) {
        auto gax = get_ghost_var(x, *ita);
        for (auto ita2 = std::next(ita); ita2 != a_coeffs.rend(); ++ita2) {
          auto gax2 = get_ghost_var(x, *ita2);
          unsigned new_a = *ita2 > *ita ? *ita2 - *ita : *ita - *ita2;
          bool a2_big = *ita2 > *ita;
          if (a_coeffs_with_one.find(new_a) != a_coeffs.end()) {
            variable_t gax_new = get_ghost_var(x, new_a);
            // a1 * x - a2 * x
            if (auto copt = m_ext_absval.difference_bound(gax2, gax)) {
              auto new_cons = a2_big ? (-gax_new <= *copt) : (gax_new <= *copt);
              CRAB_LOG("fixed-tvpi-simplifying",
                       crab::outs()
                           << "Simplifying found: " << new_cons << "\n");
              if (new_a == 1) {
                m_base_absval += new_cons;
              } else {
                m_ext_absval += new_cons;
              }
            }
            // a2 * x - a1 * x
            if (auto copt = m_ext_absval.difference_bound(gax, gax2)) {
              auto new_cons = a2_big ? (gax_new <= *copt) : (-gax_new <= *copt);
              CRAB_LOG("fixed-tvpi-simplifying",
                       crab::outs()
                           << "Simplifying found: " << new_cons << "\n");
              if (new_a == 1) {
                m_base_absval += new_cons;
              } else {
                m_ext_absval += new_cons;
              }
            }
          }
        }
      }
    }
  }

  void elimination() {
    // Based on tvpi domain, filter out constraints (C) represented by DBM
    // requires: filter(C) <= C, this does not mean filter(DBM) <= DBM.
    // filter(C) \equiv C, this indicates constraints after filtering can
    // still prove the same properties. However, for DBM, when remove redundant
    // constraints, the entails operation requires to check if the constraints
    // can be proved by the remaining constraints through the check_entailement
    // operation. It cannot simpliy by checking the DBM.
    // Thus, it needs to check whether eliminating constraints can derive
    // efficiency; otherwise, ignore this process.
    //
    // check if a1 * x - b1 * y <= c1 and a2 * x - b2 * y <= c2 implies
    // a3 * x - b3 * y <= c3
    for (auto it1 = m_coeff_map.begin(); it1 != m_coeff_map.end(); ++it1) {
      for (auto it2 = std::next(it1); it2 != m_coeff_map.end(); ++it2) {
        const variable_t &x = it1->first;
        const coefficient_set_t &a_coeffs = it1->second;
        const variable_t &y = it2->first;
        const coefficient_set_t &b_coeffs = it2->second;
      }
    }
  }

  /*
      for (auto it1 = m_coeff_map.begin(); it1 != m_coeff_map.end(); ++it1) {
        const variable_t &a = it1->first;
        const coefficient_set_t &a_coeffs = it1->second;
        for (auto it2 = std::next(it1); it2 != m_coeff_map.end(); ++it2) {
          const variable_t &b = it2->first;
          const coefficient_set_t &b_coeffs = it2->second;
          if (auto copt = m_base_absval.difference_bound(b, a)) {
            coefficient_set_t s;
            std::set_intersection(a_coeffs.begin(), a_coeffs.end(),
                                  b_coeffs.begin(), b_coeffs.end(),
                                  std::inserter(s, s.begin()));
            for (auto &c : s) {
              number_t mc = number_t(c) * *copt;
              m_ext_absval += (get_ghost_var(a, c) - get_ghost_var(b, c) <= mc);
            }
          }
        }
        for (auto c : a_coeffs) {
          m_ext_absval += (get_ghost_var(a, c) == c * a);
        }
      }
  */

  void reduce() {
    CRAB_LOG("fixed-tvpi-reduce", crab::outs()
                                      << "Before reduction: " << *this << "\n");

    // for (auto it1 = m_coeff_map.begin(); it1 != m_coeff_map.end(); ++it1) {
    //   const variable_t &a = it1->first;
    //   const coefficient_set_t &a_coeffs = it1->second;
    //   for (auto c : a_coeffs) {
    //     m_ext_absval += (get_ghost_var(a, c) == c * a);
    //   }
    // }
    scaling();
    adding();
    simplifying();
    // propagate constraints for shared variables
    auto vars = m_base_absval.vars();
    auto varst = m_ext_absval.vars();
    std::vector<variable_t> shared_vars;
    std::set_intersection(vars.begin(), vars.end(), varst.begin(), varst.end(),
                          std::back_inserter(shared_vars));
    tvpi_imple::print_vector(crab::outs(), shared_vars);
    linear_constraint_system_t cst1, cst2;
    for (const auto &v : shared_vars) {
      m_base_absval.extract(v, cst1, false);
    }
    m_ext_absval += cst1;
    for (const auto &v : shared_vars) {
      m_ext_absval.extract(v, cst2, false);
    }
    m_base_absval += cst2;
    m_ext_absval.project(varst);
    m_base_absval.project(vars);
    if (m_base_absval.is_bottom() || m_ext_absval.is_bottom()) {
      set_to_bottom();
    }
    CRAB_LOG("fixed-tvpi-reduce", crab::outs()
                                      << "After reduction: " << *this << "\n");
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
      tvpi_dbm_domain_t this2 = *this;
      this2.reduce();
      tvpi_dbm_domain_t other2 = other;
      other2.reduce();
      CRAB_LOG("fixed-tvpi-leq", crab::outs() << "is\n"
                                              << this2 << "\n<=\n"
                                              << other2 << "\n");
      return this2.m_base_absval <= other2.m_base_absval &&
             this2.m_ext_absval <= other2.m_ext_absval;
    }
  }

  void operator|=(const tvpi_dbm_domain_t &other) override {
    if (is_bottom() || other.is_top()) {
      *this = other;
    } else if (other.is_bottom() || is_top()) {
      // do nothing
    } else {
      CRAB_LOG("fixed-tvpi-join", crab::outs() << "join\n"
                                               << *this << "\nwith\n"
                                               << other << "\n");
      reduce();
      tvpi_dbm_domain_t other2 = other;
      other2.reduce();
      CRAB_LOG("fixed-tvpi-join", crab::outs() << "reduced join\n"
                                               << *this << "\nwith\n"
                                               << other2 << "\n");
      m_base_absval |= other2.m_base_absval;
      m_ext_absval |= other2.m_ext_absval;
      m_coeff_map.merge(other2.m_coeff_map);
      CRAB_LOG("fixed-tvpi-join", crab::outs() << *this << "\n");
    }
  }

  tvpi_dbm_domain_t operator|(const tvpi_dbm_domain_t &other) const override {
    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    } else {
      CRAB_LOG("fixed-tvpi-join", crab::outs() << "join\n"
                                               << *this << "\nwith\n"
                                               << other << "\n");
      tvpi_dbm_domain_t this2 = *this;
      this2.reduce();
      tvpi_dbm_domain_t other2 = other;
      other2.reduce();
      CRAB_LOG("fixed-tvpi-join", crab::outs() << "reduced join\n"
                                               << this2 << "\nwith\n"
                                               << other2 << "\n");
      this2.m_base_absval |= other2.m_base_absval;
      this2.m_ext_absval |= other2.m_ext_absval;
      this2.m_coeff_map.merge(other2.m_coeff_map);
      CRAB_LOG("fixed-tvpi-join", crab::outs() << this2 << "\n");
      return this2;
    }
  }

  void operator&=(const tvpi_dbm_domain_t &other) override {
    if (is_bottom() || other.is_top()) {
      // do nothing
    } else if (other.is_bottom() || is_top()) {
      *this = other;
    } else {
      reduce();
      tvpi_dbm_domain_t other2 = other;
      other2.reduce();
      m_base_absval &= other2.m_base_absval;
      m_ext_absval &= other2.m_ext_absval;
      m_coeff_map.merge(other2.m_coeff_map);
      if (m_base_absval.is_bottom() || m_ext_absval.is_bottom()) {
        set_to_bottom();
      }
    }
  }

  tvpi_dbm_domain_t operator&(const tvpi_dbm_domain_t &other) const override {
    if (is_bottom() || other.is_top()) {
      return *this;
    } else if (other.is_bottom() || is_top()) {
      return other;
    } else {
      tvpi_dbm_domain_t this2 = *this;
      this2.reduce();
      tvpi_dbm_domain_t other2 = other;
      other2.reduce();
      this2.m_base_absval &= other2.m_base_absval;
      this2.m_ext_absval &= other2.m_ext_absval;
      this2.m_coeff_map.merge(other2.m_coeff_map);
      if (this2.m_base_absval.is_bottom() || this2.m_ext_absval.is_bottom()) {
        this2.set_to_bottom();
      }
      return this2;
    }
  }

  tvpi_dbm_domain_t operator||(const tvpi_dbm_domain_t &other) const override {
    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    } else {
      tvpi_dbm_domain_t this2 = *this;
      this2.reduce();
      tvpi_dbm_domain_t other2 = other;
      other2.reduce();
      base_domain_t out_base_absval =
          this2.m_base_absval || other2.m_base_absval;
      base_domain_t out_ext_absval = this2.m_ext_absval || other2.m_ext_absval;
      coefficient_map_t out_coeff_map = m_coeff_map;
      out_coeff_map.merge(other2.m_coeff_map);
      return tvpi_dbm_domain_t(std::move(out_base_absval),
                               std::move(out_ext_absval),
                               std::move(out_coeff_map));
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
      base_domain_t out_base_absval =
          m_base_absval.widening_thresholds(other.m_base_absval, ts);
      base_domain_t out_ext_absval =
          m_ext_absval.widening_thresholds(other.m_ext_absval, ts);
      coefficient_map_t out_coeff_map = m_coeff_map;
      out_coeff_map.merge(other.m_coeff_map);
      return tvpi_dbm_domain_t(std::move(out_base_absval),
                               std::move(out_ext_absval),
                               std::move(out_coeff_map));
    }
  }

  tvpi_dbm_domain_t operator&&(const tvpi_dbm_domain_t &other) const override {
    if (is_bottom() || other.is_top()) {
      return *this;
    } else if (other.is_bottom() || is_top()) {
      return other;
    } else {
      base_domain_t out_base_absval = m_base_absval && other.m_base_absval;
      base_domain_t out_ext_absval = m_ext_absval && other.m_ext_absval;
      coefficient_map_t out_coeff_map = m_coeff_map;
      out_coeff_map.merge(other.m_coeff_map);
      return tvpi_dbm_domain_t(std::move(out_base_absval),
                               std::move(out_ext_absval),
                               std::move(out_coeff_map));
    }
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    CRAB_LOG("fixed-tvpi", crab::outs() << "assume(" << csts
                                        << ")\nBefore: " << *this << "\n");
    if (!is_bottom()) {
      for (auto const &cst : csts) {
        if (cst.is_contradiction()) {
          set_to_bottom();
          break;
        }

        if (cst.is_tautology()) {
          continue;
        }

        CRAB_LOG("fixed-tvpi", crab::outs()
                                   << "processing original: " << cst << "\n");

        m_base_absval += cst;
        if (m_base_absval.is_bottom()) {
          set_to_bottom();
          break;
        }
        auto ecst = rewrite_linear_constraint(cst);
        if (ecst.equal(cst)) {
          continue;
        }
        CRAB_LOG("fixed-tvpi", crab::outs()
                                   << "processing rewritten: " << ecst << "\n");
        m_ext_absval += ecst;
        if (m_ext_absval.is_bottom()) {
          set_to_bottom();
          break;
        }

      } // end for
    }

    CRAB_LOG("fixed-tvpi", crab::outs() << "assume(" << csts
                                        << ")\nAfter: " << *this << "\n");
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
        tmp.reduce();
        CRAB_LOG("fixed-tvpi", crab::outs() << "reduced: " << tmp << "\n");

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
    if (!is_bottom()) {
      CRAB_LOG("fixed-tvpi", crab::outs() << "Before assign(" << x << " := "
                                          << e << ")=" << *this << "\n");
      CRAB_LOG("fixed-tvpi", crab::outs() << "processing original: " << x
                                          << " := " << e << "\n");
      m_base_absval.assign(x, e);
      auto ex = e.get_variable();
      if (ex && *(ex) == x) {
        // heuristics, based on syntax of expression, x := x +/- c
        // probably some index or counter, add it to extended dbm
        m_ext_absval.assign(x, e);
      }

      linear_expression_t e1 = rewrite_linear_expression(e);
      if (!e1.equal(e)) {
        CRAB_LOG("fixed-tvpi", crab::outs() << "processing rewritten" << x
                                            << " := " << e1 << "\n");
        m_ext_absval.assign(x, e1);
      }

      auto it = m_coeff_map.find(x);
      if (it == m_coeff_map.end()) {
        return;
      }
      auto &coeffs = it->second;
      for (auto coefficient : coeffs) {
        rewrite_assign(x, e, coefficient, false /*!weak*/);
      }

      CRAB_LOG("fixed-tvpi", crab::outs() << "After assign(" << x << " := " << e
                                          << ")=" << *this << "\n");
    }
  }

  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    if (!is_bottom()) {
      m_base_absval.weak_assign(x, e);

      linear_expression_t e1 = rewrite_linear_expression(e);
      if (!e1.equal(e)) {
        m_ext_absval.weak_assign(x, e1);
      }

      auto it = m_coeff_map.find(x);
      if (it == m_coeff_map.end()) {
        return;
      }
      auto &coeffs = it->second;
      for (auto coefficient : coeffs) {
        rewrite_assign(x, e, coefficient, true /*weak*/);
      }
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
        if (x == y) { // x := x +/- z
          // hint: this might be an index or a counter, add it to extended
          // value
          m_ext_absval.apply(op, x, y, z);
        }
        break;
      case OP_MULTIPLICATION: // x := y * z
        if (z_abs > number_t(1)) {
          // "x := zy for z > 1"
          m_ext_absval.apply(op, x, get_ghost_var(y, z_abs), sign_one);
        }
        break;
      case OP_SDIV: // x := y /s z
      case OP_UDIV: // x := y /u z
        if (z_abs > number_t(1)) {
          // "zx := y for z > 1"
          m_ext_absval.apply(op, get_ghost_var(x, z_abs), y, sign_one);
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

      auto it = m_coeff_map.find(x);
      if (it == m_coeff_map.end()) {
        return;
      }
      auto &coeffs = it->second;
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
      auto it = m_coeff_map.find(x);
      if (it == m_coeff_map.end()) {
        return;
      }
      auto &coeffs = it->second;
      for (auto coefficient : coeffs) {
        rewrite_apply(op, x, y, z, coefficient);
      }
    }
  }

  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    if (!is_bottom()) {
      m_base_absval.apply(op, dst, src);
    }
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    if (!is_bottom()) {
      m_base_absval.apply(op, x, y, z);
    }
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    if (!is_bottom()) {
      m_base_absval.apply(op, x, y, k);
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
      m_base_absval -= var;
      m_ext_absval -= var;

      auto it = m_coeff_map.find(var);
      if (it == m_coeff_map.end()) {
        return;
      }
      auto &coeffs = it->second;
      for (auto coefficient : coeffs) {
        variable_t ghost_var = get_ghost_var(var, coefficient);
        m_ext_absval -= ghost_var;
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

  void forget(const variable_vector_t &variables) override {
    if (!(is_bottom() || is_top())) {
      m_base_absval.forget(variables);
      variable_vector_t allvars(variables);
      for (auto const &v : variables) {
        auto it = m_coeff_map.find(v);
        if (it != m_coeff_map.end()) {
          for (auto coefficient : it->second) {
            variable_t gv = get_ghost_var(v, coefficient);
            allvars.push_back(gv);
          }
        }
      }
      m_ext_absval.forget(allvars);
    }
  }

  void project(const variable_vector_t &variables) override {
    if (!is_bottom()) {
      m_base_absval.project(variables);
      variable_vector_t allvars(variables);
      for (auto const &v : variables) {
        auto it = m_coeff_map.find(v);
        if (it != m_coeff_map.end()) {
          for (auto coefficient : it->second) {
            variable_t gv = get_ghost_var(v, coefficient);
            allvars.push_back(gv);
          }
        }
      }
      m_ext_absval.project(allvars);
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
        auto it = m_coeff_map.find(f);
        if (it != m_coeff_map.end()) {
          for (auto coefficient : it->second) {
            extd_from.push_back(get_ghost_var(f, coefficient));
            extd_to.push_back(get_ghost_var(t, coefficient));
          }
          auto set2 = it->second;
          // m_coeff_map.erase(it);
          m_coeff_map.insert({t, std::move(set2)});
        }
      }
      m_ext_absval.rename(extd_from, extd_to);
    }
  }

  void expand(const variable_t &var, const variable_t &new_var) override {
    if (is_bottom() || is_top()) {
      return;
    }

    m_base_absval.expand(var, new_var);

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
        o << "map:";
        m_coeff_map.write(o);
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
