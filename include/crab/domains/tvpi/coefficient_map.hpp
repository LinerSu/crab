#pragma once
#include <algorithm>
#include <boost/optional.hpp>
#include <chrono>
#include <string>

#include <crab/numbers/bignums.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>

namespace crab {
namespace domains {

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

/// @brief a special log method to print unordered set
/// @tparam TType
/// @param o crab ostream
/// @param vec the vector for printing
template <typename TType>
void print_set(crab::crab_os &o, const std::unordered_set<TType> &s) {
  typename std::unordered_set<TType>::const_iterator it;
  o << "(";
  for (it = s.begin(); it != s.end(); it++) {
    if (it != s.begin())
      o << ",";
    o << (*it);
  }
  o << ")";
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

  iterator_t erase(iterator_t it) {
    ensure_exists();
    ensure_unique();
    return m_map->erase(it);
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
      unions(o);
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
} // end namespace domains
} // end namespace crab