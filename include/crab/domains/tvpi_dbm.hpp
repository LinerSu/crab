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

namespace tvpi_utils {
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

/// @brief a special method for find a value in a vector
/// @tparam T a type of the vector item
/// @param vec the vector for searching
/// @param value the value for searching
/// @return item if found; otherwise, boost::none
template <typename T>
boost::optional<T> find(const std::vector<T> &vec, const T &value) {
  auto it = std::find(vec.begin(), vec.end(), value);
  if (it != vec.end()) {
    return *it;
  } else {
    return boost::none;
  }
}

/// @brief compute gcd value by two unsigned ints
/// @tparam T a type of the number, we used unsigned int only
/// @param a unsigned int
/// @param b unsigned int
/// @return the gcd value; if a or b is 0, return another value
template <typename T> T gcd(T a, T b) {
  // Continue until b becomes zero
  while (b != 0) {
    T temp = b;
    b = a % b;
    a = temp;
  }
  return a;
}

template <typename T> void intersect(std::set<T> &a, const std::set<T> &b) {
  for (auto it = a.begin(); it != a.end();) {
    if (b.find(*it) == b.end()) { // not found, remove
      it = a.erase(it);
    } else {
      ++it;
    }
  }
}

/**
 * @class coefficient_map
 * @brief A class template that manages a map of keys to sets of values.
 *
 * This class provides a map where each key is associated with a set of values.
 * It supports insertion, merging, and various other operations on the map.
 *
 * @tparam key The type of the keys in the map.
 * @tparam value The type of the values in the sets.
 */
template <typename key, typename value> class coefficient_map {
  // This is an abstract domain for tracking coefficient set
  // for each variable. The map is similar to separate domain where key is
  // variable and value is a set of coefficients.
  /* clang-format off */
  // Domain hierarchy:
  // For each coefficient set U, the lattice diagram is:
  //               top : a special case for top
  //              /   \
  //           ...    ...
  //           {a,b,c,d,e,...}
  //            /   |   \
  //       {a,b,c,d} ... {b,c,d,e}
  //          /  \         /   \
  //       {a,b} ...     {d,e} ...
  //       /  \           /  \
  //      {a} {b}  {c}  {d} {e}
  //       ...  \   |   /    ...
  //               bot

  // For coefficient map, the lattice diagram is:
  //               top: { x: top, y: top }
  //              /   \
  //           ...    ...
  //       {x: U_x, y: U_y}
  //        /      \
  //   {x: U_x} {y: U_y} ...
  //       \       /     /
  //           bot

  /* clang-format on */
  using key_t = key;
  using value_t = value;
  using key_value_t = std::pair<const key_t, value_t>;
  using value_set_t = std::set<value_t>;
  using key_value_set_t = std::pair<const key_t, value_set_t>;
  using map_t = std::unordered_map<key_t, value_set_t>;
  using const_iterator_t = typename map_t::const_iterator;
  using iterator_t = typename map_t::iterator;
  using shared_map_t = std::shared_ptr<map_t>;
  enum class lattice_val { bottom, top, neither_top_nor_bot };

public:
  using coefficient_map_t = coefficient_map<key, value>;
  using coefficient_set_t = value_set_t;

  coefficient_map() : m_map(nullptr), m_val(lattice_val::neither_top_nor_bot) {}
  coefficient_map(const coefficient_map_t &o) = default;
  coefficient_map(coefficient_map_t &&o) = default;
  coefficient_map_t &operator=(const coefficient_map_t &o) = default;
  coefficient_map_t &operator=(coefficient_map_t &&o) = default;

  /**
   * @brief Returns map is empty of not.
   * @return true if map is not exist or empty.
   */
  bool empty() const { return !exists() || m_map->empty(); }

  void clear() {
    if (exists() && m_map.unique()) {
      m_map->clear();
    }
    m_map = nullptr;
  }

  bool is_bottom() const { return m_val == lattice_val::bottom; }

  void set_to_bottom() {
    clear();
    m_val = lattice_val::bottom;
  }

  bool is_top() const { return m_val == lattice_val::top; }

  void set_to_top() {
    clear();
    m_val = lattice_val::top;
  }

  /**
   * @brief Inserts a key-value pair into the map.
   * @param kv The key-value pair to insert.
   */
  void insert(const key_value_t &kv) {
    ensure_exists();
    ensure_unique();
    auto &m = *m_map;
    auto it = m.find(kv.first);
    if (it != m.end()) {
      it->second.insert(kv.second);
    } else {
      m.insert({kv.first, {kv.second}});
    }
  }

  /**
   * @brief Inserts a key-value pair into the map using move semantics.
   * @param kv The key-value pair to insert.
   */
  void insert(key_value_t &&kv) {
    ensure_exists();
    ensure_unique();
    auto &m = *m_map;
    auto it = m.find(kv.first);
    if (it != m.end()) {
      it->second.insert(std::forward<key_value_t>(kv).second);
    } else {
      m.insert({std::forward<key_value_t>(kv).first,
                {std::forward<key_value_t>(kv).second}});
    }
  }

  /**
   * @brief Inserts a key and its associated set of values into the map.
   * @param kvs The key and its associated set of values to insert.
   */
  void insert(const key_value_set_t &kvs) {
    ensure_exists();
    ensure_unique();
    auto &m = *m_map;
    auto it = m.find(kvs.first);
    if (it != m.end()) {
      it->second.insert(kvs.second.begin(), kvs.second.end());
    } else {
      m.insert(kvs);
    }
  }

  /**
   * @brief Inserts a key and its associated set of values into the map using
   * move semantics.
   * @param kvs The key and its associated set of values to insert.
   */
  void insert(key_value_set_t &&kvs) {
    ensure_exists();
    ensure_unique();
    auto &m = *m_map;
    auto it = m.find(kvs.first);
    if (it != m.end()) {
      it->second.insert(std::make_move_iterator(kvs.second.begin()),
                        std::make_move_iterator(kvs.second.end()));
    } else {
      m.insert(std::move(kvs));
    }
  }

  /**
   * @brief Finds a key in the map.
   * @param k The key to find.
   * @return A constant iterator to the key-value pair if found, otherwise
   * end().
   */
  const_iterator_t find(const key_t &k) const {
    ensure_exists();
    return m_map->find(k);
  }

  /**
   * @brief Finds a key in the map.
   * @param k The key to find.
   * @return An iterator to the key-value pair if found, otherwise end().
   */
  iterator_t find(const key_t &k) {
    ensure_exists();
    ensure_unique();
    return m_map->find(k);
  }

  /**
   * @brief Returns a constant iterator to the beginning of the map.
   * @return A constant iterator to the beginning of the map.
   */
  const_iterator_t begin() const {
    ensure_exists();
    return m_map->begin();
  }

  /**
   * @brief Returns an iterator to the beginning of the map.
   * @return An iterator to the beginning of the map.
   */
  iterator_t begin() {
    ensure_exists();
    ensure_unique();
    return m_map->begin();
  }

  /**
   * @brief Returns a constant iterator to the end of the map.
   * @return A constant iterator to the end of the map.
   */
  const_iterator_t end() const {
    ensure_exists();
    return m_map->end();
  }

  /**
   * @brief Returns an iterator to the end of the map.
   * @return An iterator to the end of the map.
   */
  iterator_t end() {
    ensure_exists();
    ensure_unique();
    return m_map->end();
  }

  void remove(const key_t &k) {
    if (!exists() || m_map->find(k) == m_map->end()) {
      return;
    }
    ensure_unique();
    m_map->erase(k);
  }

  void keep(const std::vector<key_t> &keys) {
    if (!exists()) {
      return;
    }
    ensure_unique();
    for (auto it = m_map->begin(); it != m_map->end();) {
      if (tvpi_utils::find(keys, it->first)) { // keep key if found
        ++it;
      } else { // remove
        it = m_map->erase(it);
      }
    }
  }

  /**
   * @brief Merges another coefficient_map by inserting missing key value pair
   * from another map.
   * @param o The coefficient_map to merge from.
   */
  void unions(const coefficient_map_t &o) {
    if (!exists()) { // current is uninitialized, make it initialized
      if (o.exists()) {
        m_map = std::make_shared<map_t>(*o.m_map);
        ;
      }
      return;
    }
    if (!o.exists()) {
      return;
    }

    if (m_map == o.m_map) { // if they are the same, ignore
      return;
    }

    ensure_unique();

    for (auto it = o.begin(), et = o.end(); it != et; ++it) {
      const auto &k = it->first;
      insert({k, it->second});
    }
  }

  void meet(const coefficient_map_t &o) {
    if (is_bottom() || o.is_top()) {
      // do nothing
    } else if (o.is_bottom() || is_top()) {
      *this = o;
    } else {
      unions(o);
    }
  }

  /**
   * @brief Merges another coefficient_map by keep only the common key value
   * pair from another map.
   * @param o The coefficient_map to merge from.
   */
  void intersects(const coefficient_map_t &o) {
    if (!exists()) { // current is uninitialized, make it initialized
      if (o.exists()) {
        m_map = std::make_shared<map_t>(*o.m_map);
        ;
      }
      return;
    }
    if (!o.exists()) {
      return;
    }

    if (m_map == o.m_map) { // if they are the same, ignore
      return;
    }

    ensure_unique();

    for (auto it = m_map->begin(), et = m_map->end(); it != et;) {
      auto oit = o.find(it->first);
      if (oit == o.end()) {
        it = m_map->erase(it);
      } else { // common key
        intersect(it->second, oit->second);
        if (it->second.empty()) {
          it = m_map->erase(it);
        } else {
          ++it;
        }
      }
    }
  }

  void join(const coefficient_map_t &o) {
    if (is_bottom() || o.is_top()) {
      *this = o;
    } else if (o.is_bottom() || is_top()) {
      // do nothing
    } else {
      intersects(o);
    }
  }

  /**
   * @brief Writes the map to an output stream.
   * @param o The output stream to write to.
   */
  void write(crab_os &o) const {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "top";
    } else if (!exists()) {
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

  /**
   * @brief Output stream operator for coefficient_map.
   * @param o The output stream.
   * @param m The coefficient_map to write.
   * @return The output stream.
   */
  friend crab_os &operator<<(crab_os &o, const coefficient_map_t &m) {
    m.write(o);
    return o;
  }

private:
  shared_map_t
      m_map; // Shared pointer to the underlying map, enabling state sharing.
  lattice_val m_val; // The lattice value to indicate the state of the map.

  /**
   * @brief Ensures that the map exists, creating it if necessary.
   */
  void ensure_exists() {
    if (!m_map) {
      m_map = std::make_shared<map_t>();
    }
  }

  /**
   * @brief Ensures that the map exists, throwing an error if it does not.
   */
  void ensure_exists() const {
    if (!m_map) {
      CRAB_ERROR("shared map does not exist!");
    }
  }

  /**
   * @brief Checks if the map exists.
   * @return True if the map exists, false otherwise.
   */
  bool exists() const { return m_map != nullptr; }

  /**
   * @brief Ensures that the map is written when it is copied.
   */
  void ensure_unique() {
    if (!m_map.unique()) {
      m_map = std::make_shared<map_t>(*m_map);
    }
  }
};

} // namespace tvpi_utils

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
      typename tvpi_utils::coefficient_map<variable_t,
                                           unsigned>::coefficient_map_t;
  using coefficient_set_t = typename coefficient_map_t::coefficient_set_t;

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

  base_domain_t m_base_absval;
  base_domain_t m_ext_absval;
  coefficient_map_t m_coeff_map;
  boost::optional<variable_t> counter;

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
                  coefficient_map_t &&coeff_map)
      : m_base_absval(std::move(base)), m_ext_absval(std::move(extd)),
        m_coeff_map(std::move(coeff_map)) {
    // tvpi_reduce();
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
    m_coeff_map.set_to_top();
    m_base_absval.set_to_top();
    m_ext_absval.set_to_top();
  }

  void set_to_bottom() override {
    m_coeff_map.set_to_bottom();
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
    return m_coeff_map.is_bottom() || m_base_absval.is_bottom() ||
           m_ext_absval.is_bottom();
  }

  bool is_top() const override {
    return m_coeff_map.is_top() && m_base_absval.is_top() &&
           m_ext_absval.is_top();
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

  void normalize_dbms() {
    // This is a top level function to normalize tvpi constraints.
    // However, we wish not to run this process in the end since we always keep
    // inequalities after normalized.
    // Just in case we need in the future, this is the implementation.
    for (auto it1 = m_coeff_map.begin(); it1 != m_coeff_map.end(); ++it1) {
      const variable_t &x = it1->first;
      const coefficient_set_t &a_coeffs = it1->second;
      for (auto ita = a_coeffs.begin(); ita != a_coeffs.end(); ++ita) {
        auto gax = get_ghost_var(x, *ita);
        for (auto it2 = std::next(it1); it2 != m_coeff_map.end(); ++it2) {
          const variable_t &y = it2->first;
          const coefficient_set_t &b_coeffs = it2->second;
          for (auto itb = b_coeffs.begin(); itb != b_coeffs.end(); ++itb) {
            auto gby = get_ghost_var(y, *itb);
            auto copt = m_ext_absval.difference_bound(gby, gax);
            if (copt) {
              auto ret = normalize_tvpi(*ita, *itb, *copt);
              auto gcd = ret.first;
              unsigned new_a = *ita / gcd;
              unsigned new_b = *itb / gcd;
              number_t new_c = ret.second;
              if (new_a == 1 && new_b == 1) {
                m_base_absval += (x - y <= new_c);
              } else {
                auto gax_new = get_ghost_var(x, new_a);
                auto gby_new = get_ghost_var(y, new_b);
                m_ext_absval += (gax_new - gby_new <= new_c);
              }
              // TODO: Remove old constraint
            }
          }
        }
      }
    }
  }

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
      if (tvpi_utils::find(base_vars, *it)) {
        it = to_remove.erase(it);
      } else {
        ++it;
      }
    }
    m_ext_absval.forget(variable_vector_t(to_remove.begin(), to_remove.end()));
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
    CRAB_LOG("tvpi-dbm-reduce2", crab::outs() << "extend dimensions: ";
             tvpi_utils::print_vector(crab::outs(), ext_vars);
             crab::outs() << "\n";);
    auto e_set = std::set<variable_t>(base_vars.begin(), base_vars.end());
    for (auto &v : base_vars) {
      auto itc = m_coeff_map.find(v);
      if (itc != m_coeff_map.end()) {
        for (auto &c : itc->second) {
          auto gv = get_ghost_var(v, c);
          if (tvpi_utils::find(ext_vars, v)) {
            e_set.insert(v);
          }
        }
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
      auto itc1 = m_coeff_map.find(x);
      const coefficient_set_t &a_coeffs =
          itc1 != m_coeff_map.end() ? itc1->second : no_coeff;
      for (auto it2 = traverse_vars.cbegin(); it2 != traverse_vars.cend();
           ++it2) {
        const variable_t &y = *it2;
        if (x == y) {
          continue;
        }
        if (tvpi_utils::find(base_vars, y) == boost::none)
          continue;
        auto itc2 = m_coeff_map.find(y);
        const coefficient_set_t &b_coeffs =
            itc2 != m_coeff_map.end() ? itc2->second : no_coeff;
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
                  CRAB_LOG("tvpi-dbm-reduce2",
                           crab::outs()
                               << "=>>>" << new_a << "*" << x << "-" << new_b
                               << "*" << z << "<=" << new_c << "\n");
                  if (new_a == 0 && new_b == 0) {
                    if (new_c < 0) { // 0 <= c where c < 0, UNSAT
                      set_to_bottom();
                      return;
                    }
                  } else if (new_a == 1 && new_b == 1) {
                    m_base_absval += (x - z <= new_c);
                  } else if (new_a == 1 && new_b == 0) {
                    m_base_absval += (x <= new_c);
                  } else if (new_a == 0 && new_b == 1) {
                    m_base_absval += (-z <= new_c);
                  } else {
                    if (find_ghost_var(x, new_a) && find_ghost_var(z, new_b)) {
                      auto gax_new = get_ghost_var(x, new_a);
                      auto gbz_new = get_ghost_var(z, new_b);
                      if (tvpi_utils::find(ext_vars, gax_new) &&
                          tvpi_utils::find(ext_vars, gbz_new)) {
                        m_ext_absval += (gax_new - gbz_new <= new_c);
                      }
                    }
                  }
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
                  CRAB_LOG("tvpi-dbm-reduce2",
                           crab::outs() << "=>>>" << new_a << "*" << x
                                        << " <= " << new_c << "\n");
                  m_base_absval += (x <= new_c);
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
                  CRAB_LOG("tvpi-dbm-reduce2", crab::outs()
                                                   << "=>>>"
                                                   << "-" << new_b << "*" << y
                                                   << " <= " << new_c << "\n");
                  m_base_absval += (-y <= new_c);
                }
              } else if (x != y && !(a == 1 && b == 1) && counter &&
                         (x == *counter || y == *counter)) {
                // see whether we can recover this
                // x <= c && -y <= f => ax - by <= c'
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
                    CRAB_LOG("tvpi-dbm-reduce2",
                             crab::outs() << "=>>>" << gax << "-" << gby
                                          << " <= " << c_p << "\n");
                    m_ext_absval += (gax - gby <= c_p);
                  }
                }
              }
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
                  CRAB_LOG("tvpi-dbm-reduce2",
                           crab::outs()
                               << "=>>>" << new_a << "*" << z << "-" << new_b
                               << "*" << x << " <= " << new_c << "\n");
                  if (new_a == 0 && new_b == 0) {
                    if (new_c < 0) { // 0 <= c where c < 0, UNSAT
                      set_to_bottom();
                      return;
                    }
                  } else if (new_a == 1 && new_b == 1) {
                    m_base_absval += (z - x <= new_c);
                  } else if (new_a == 1 && new_b == 0) {
                    m_base_absval += (z <= new_c);
                  } else if (new_a == 0 && new_b == 1) {
                    m_base_absval += (-x <= new_c);
                  } else {
                    if (find_ghost_var(z, new_a) && find_ghost_var(x, new_b)) {
                      auto gaz_new = get_ghost_var(z, new_a);
                      auto gbx_new = get_ghost_var(x, new_b);
                      if (tvpi_utils::find(ext_vars, gbx_new) &&
                          tvpi_utils::find(ext_vars, gaz_new)) {
                        m_ext_absval += (gaz_new - gbx_new <= new_c);
                      }
                    }
                  }
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
                  CRAB_LOG("tvpi-dbm-reduce2",
                           crab::outs() << "=>>>" << new_a << "*" << y
                                        << "<=" << new_c << "\n");
                  m_base_absval += (y <= new_c);
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
                  CRAB_LOG("tvpi-dbm-reduce2", crab::outs()
                                                   << "=>>>"
                                                   << "-" << new_b << "*" << x
                                                   << " <= " << new_c << "\n");
                  m_base_absval += (-x <= new_c);
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
                    CRAB_LOG("tvpi-dbm-reduce2",
                             crab::outs() << "=>>>" << gby << "-" << gax
                                          << " <= " << c_p << "\n");
                    m_ext_absval += (gby - gax <= c_p);
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
      this2.m_coeff_map.meet(other.m_coeff_map);
      this2.tvpi_reduce();
      tvpi_dbm_domain_t other2 = other;
      other2.m_coeff_map.meet(m_coeff_map);
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
      m_coeff_map.join(other.m_coeff_map);
      tvpi_reduce();
      tvpi_dbm_domain_t other2 = other;
      other2.m_coeff_map.join(m_coeff_map);
      other2.tvpi_reduce();
      CRAB_LOG("tvpi-dbm-join", crab::outs() << "[reduced join]\n"
                                             << *this << "\nwith\n"
                                             << other2 << "\n");
      m_base_absval |= other2.m_base_absval;
      m_ext_absval |= other2.m_ext_absval;
      filter_dbms();
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
      this2.m_coeff_map.join(other.m_coeff_map);
      this2.tvpi_reduce();
      tvpi_dbm_domain_t other2 = other;
      other2.m_coeff_map.join(m_coeff_map);
      other2.tvpi_reduce();
      CRAB_LOG("tvpi-dbm-join", crab::outs() << "[reduced join]\n"
                                             << this2 << "\nwith\n"
                                             << other2 << "\n");
      this2.m_base_absval |= other2.m_base_absval;
      this2.m_ext_absval |= other2.m_ext_absval;
      this2.filter_dbms();
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
      m_coeff_map.meet(other.m_coeff_map);
      tvpi_reduce();
      tvpi_dbm_domain_t other2 = other;
      other2.m_coeff_map.meet(m_coeff_map);
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
      this2.m_coeff_map.meet(other.m_coeff_map);
      this2.tvpi_reduce();
      tvpi_dbm_domain_t other2 = other;
      other2.m_coeff_map.meet(m_coeff_map);
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
      this2.m_coeff_map.join(other.m_coeff_map);
      this2.tvpi_reduce();
      tvpi_dbm_domain_t other2 = other;
      other2.m_coeff_map.join(m_coeff_map);
      other2.tvpi_reduce();
      CRAB_LOG("tvpi-dbm-widen", crab::outs() << "[reduced widen]\n"
                                              << this2 << "\nwith\n"
                                              << other2 << "\n");
      base_domain_t out_base_absval =
          this2.m_base_absval || other2.m_base_absval;
      base_domain_t out_ext_absval = this2.m_ext_absval || other2.m_ext_absval;
      coefficient_map_t out_coeff_map = this2.m_coeff_map;
      tvpi_dbm_domain_t res(std::move(out_base_absval),
                            std::move(out_ext_absval),
                            std::move(out_coeff_map));
      res.filter_dbms();
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
      this2.m_coeff_map.join(other.m_coeff_map);
      this2.tvpi_reduce();
      tvpi_dbm_domain_t other2 = other;
      other2.m_coeff_map.join(m_coeff_map);
      other2.tvpi_reduce();
      CRAB_LOG("tvpi-dbm-widen", crab::outs() << "[reduced widen]\n"
                                              << this2 << "\nwith\n"
                                              << other2 << "\n");
      base_domain_t out_base_absval =
          this2.m_base_absval.widening_thresholds(other2.m_base_absval, ts);
      base_domain_t out_ext_absval =
          this2.m_ext_absval.widening_thresholds(other2.m_ext_absval, ts);
      coefficient_map_t out_coeff_map = this2.m_coeff_map;
      tvpi_dbm_domain_t res(std::move(out_base_absval),
                            std::move(out_ext_absval),
                            std::move(out_coeff_map));
      res.filter_dbms();
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
      this2.m_coeff_map.meet(other.m_coeff_map);
      this2.tvpi_reduce();
      tvpi_dbm_domain_t other2 = other;
      other2.m_coeff_map.meet(m_coeff_map);
      other2.tvpi_reduce();
      CRAB_LOG("tvpi-dbm-narrow", crab::outs() << "[reduced narrow]\n"
                                               << this2 << "\nwith\n"
                                               << other2 << "\n");
      base_domain_t out_base_absval =
          this2.m_base_absval && other2.m_base_absval;
      base_domain_t out_ext_absval = this2.m_ext_absval && other2.m_ext_absval;
      coefficient_map_t out_coeff_map = this2.m_coeff_map;
      tvpi_dbm_domain_t res(std::move(out_base_absval),
                            std::move(out_ext_absval),
                            std::move(out_coeff_map));
      CRAB_LOG("tvpi-dbm-narrow", crab::outs() << "[res]\n" << res << "\n");
      return res;
    }
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    CRAB_LOG("tvpi-dbm", crab::outs() << "assume(" << csts
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

        CRAB_LOG("tvpi-dbm", crab::outs()
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
        CRAB_LOG("tvpi-dbm", crab::outs()
                                 << "processing rewritten: " << ecst << "\n");
        m_ext_absval += ecst;
        if (m_ext_absval.is_bottom()) {
          set_to_bottom();
          break;
        }

      } // end for
    }

    CRAB_LOG("tvpi-dbm", crab::outs() << "assume(" << csts
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
    if (!is_bottom()) {
      CRAB_LOG("tvpi-dbm", crab::outs() << "Before assign(" << x << " := " << e
                                        << ")=" << *this << "\n");
      CRAB_LOG("tvpi-dbm", crab::outs() << "processing original: " << x
                                        << " := " << e << "\n");
      m_base_absval.assign(x, e);
      auto ex = e.get_variable();
      // if (ex && *(ex) == x) {
      //   // heuristics, based on syntax of expression, x := x +/- c
      //   // probably some index or counter, add it to extended dbm
      //   m_ext_absval.assign(x, e);
      // }

      linear_expression_t e1 = rewrite_linear_expression(e);
      if (!e1.equal(e)) {
        CRAB_LOG("tvpi-dbm", crab::outs() << "processing rewritten" << x
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

      CRAB_LOG("tvpi-dbm", crab::outs() << "After assign(" << x << " := " << e
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
        if (x == y && counter && *counter != x) { // x := x +/- z
          // hint: this might be an index or a counter, add it to extended
          // value
          m_coeff_map.insert({*counter, convert(z)});
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
      m_coeff_map.remove(var);
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
          m_coeff_map.remove(v);
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
      m_coeff_map.keep(variables);
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
          m_coeff_map.insert({t, std::move(set2)});
          m_coeff_map.remove(f);
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
