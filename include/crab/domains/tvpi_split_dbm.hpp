#pragma once

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_params.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/backward_assign_operations.hpp>
#include <crab/domains/graphs/graph_config.hpp>
#include <crab/domains/graphs/graph_ops.hpp>
#include <crab/domains/graphs/graph_views.hpp>
#include <crab/domains/inter_abstract_operations.hpp>
#include <crab/domains/interval.hpp>
#include <crab/domains/tvpi/coefficient_map.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

#include <boost/optional.hpp>
#include <unordered_set>

#define JOIN_CLOSE_AFTER_MEET
// #define CHECK_POTENTIAL
// #define TVPIDBM_NO_NORMALIZE
#define USE_FLAT_MAP

#ifdef USE_FLAT_MAP
#include <boost/container/flat_map.hpp>
#else
// Operations like rename are much faster using unordered_map
#include <unordered_map>
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)
#define LOCATION_STRING (__FILE__ ":" STR(__LINE__))

namespace crab {
namespace domains {

namespace tvpi_utils {
#pragma region CoeffKey
template <class Variable> class var_coeff_key {
  Variable v; // program variable, no change, just keep it
  unsigned c; // coefficient for the variable
  using class_t = var_coeff_key<Variable>;

public:
  var_coeff_key(const Variable &_v, unsigned _c) : v(_v), c(_c) {}

  var_coeff_key(const Variable &_v) : v(_v), c(1) {}

  var_coeff_key(const class_t &o) = default;

  var_coeff_key(class_t &&o) = default;

  class_t &operator=(const class_t &o) = default;

  class_t &operator=(class_t &&o) = default;

  const Variable &var() const { return v; }

  Variable var() { return v; }

  unsigned coeff() const { return c; }

  bool operator==(const class_t &other) const {
    return v == other.v && c == other.c;
  }

  bool operator<(const class_t &other) const {
    if (v < other.v)
      return true;
    else if (v == other.v)
      return c < other.c;
    else
      return false;
  }

  friend crab_os &operator<<(crab_os &o, const class_t &k) {
    if (k.coeff() == 1) {
      o << k.var();
    } else {
      o << k.coeff() << "*" << k.var();
    }
    return o;
  }
};

// Custom hash function for var_coeff_key
template <class Variable> struct var_coeff_hash {
  std::size_t operator()(const var_coeff_key<Variable> &key) const {
    size_t hash_var = key.var().hash();
    size_t hash_coeff = std::hash<unsigned>()(key.coeff());
    return hash_var ^ (hash_coeff << 1);
  }
};
#pragma endregion CoeffKey

#pragma region ConvertNumber
template <typename Number> struct NtoC {
  static unsigned convert(const Number &coeff) {
    if (coeff < std::numeric_limits<int64_t>::min() ||
        coeff > std::numeric_limits<int64_t>::max()) {
      CRAB_ERROR("Coefficient out of bounds");
    }
    int64_t tmp = static_cast<int64_t>(coeff);
    return static_cast<unsigned>((tmp > 0 ? tmp : -tmp));
  }
};

template <typename Weight> struct WtoC {
  static unsigned convert(const Weight &coeff) {
    int64_t tmp = static_cast<int64_t>(coeff);
    return static_cast<unsigned>((tmp > 0 ? tmp : -tmp));
  }
};
#pragma endregion ConvertNumber

#pragma region TVPIUtils
template <class Number, class Variable> struct tvpi_op {
  static std::pair<unsigned, Number> normalize_tvpi(const unsigned &a,
                                                    const Number &c) {
    // Normalize TVPI constraints
    // For each constraint ax <= c, normalize it to x <= c' where c' = c / a
    return {1, c / Number(a)};
  }

  static std::pair<unsigned, Number>
  normalize_tvpi(const unsigned &a, const unsigned &b, const Number &c) {
    // Normalize TVPI constraints
    // For each constraint ax - by <= c, normalize it to a'x - b'y <= c'
    // where a' = a / gcd(a, b), b' = b / gcd(a, b), c' = c / gcd(a, b)
    unsigned gcd = tvpi_utils::gcd(a, b);
    return {gcd, c / Number(gcd)};
  }

  static std::tuple<unsigned, unsigned, Number>
  resultant(const unsigned &a, const Variable &x, const unsigned &b,
            const Variable &y, const Number &c, const unsigned &d,
            const unsigned &e, const boost::optional<Variable> &z,
            const Number &f) {
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

    auto c_p = c * Number(lambda1) + f * Number(lambda2);
    // Normalize the result
    if (z == boost::none || x == *z) {
      int a_p = (z == boost::none) ? lambda1 * a : lambda1 * a - lambda2 * e;
      if (a_p == 0) {
        return {0, 0, c_p};
      }
      unsigned abs_a = static_cast<unsigned>(std::abs(a_p));
      bool neg = (a_p) < 0;
      auto ret = tvpi_op::normalize_tvpi(abs_a, c_p);
      unsigned new_a = neg ? 0 : 1;
      unsigned new_b = neg ? 1 : 0;
      Number new_c = ret.second;
      return {new_a, new_b, new_c};
    } else {
      auto ret = tvpi_op::normalize_tvpi(lambda1 * a, lambda2 * e, c_p);
      unsigned gcd2 = ret.first;
      unsigned new_a = lambda1 * a / gcd2;
      unsigned new_b = lambda2 * e / gcd2;
      Number new_c = ret.second;
      return {new_a, new_b, new_c};
    }
  }

  static std::tuple<unsigned, unsigned, Number>
  resultant(const unsigned &b, const Variable &y, const unsigned &a,
            const Variable &x, const Number &c, const unsigned &e,
            const boost::optional<Variable> &z, const unsigned &d,
            const Number &f) {
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

    auto c_p = c * Number(lambda1) + f * Number(lambda2);
    // Normalize the result
    if (z == boost::none || x == *z) {
      int a_p = (z == boost::none) ? -lambda1 * a : lambda2 * e - lambda1 * a;
      if (a_p == 0) {
        return {0, 0, c_p};
      }
      unsigned abs_a = static_cast<unsigned>(std::abs(a_p));
      bool neg = (a_p) < 0;
      auto ret = tvpi_op::normalize_tvpi(abs_a, c_p);
      unsigned new_a = neg ? 0 : 1;
      unsigned new_b = neg ? 1 : 0;
      Number new_c = ret.second;
      return {new_a, new_b, new_c};
    } else {
      auto ret = tvpi_op::normalize_tvpi(lambda2 * e, lambda1 * a, c_p);
      unsigned gcd2 = ret.first;
      unsigned new_a = lambda2 * e / gcd2;
      unsigned new_b = lambda1 * a / gcd2;
      Number new_c = ret.second;
      return {new_a, new_b, new_c};
    }
  }
};
#pragma endregion TVPIUtils
} // namespace tvpi_utils

class TVPISplitDBMDefaultParams {
public:
  enum { implement_inter_transformers = 0 };
};

#define TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(NAME)                               \
  CRAB_DOMAIN_SCOPED_STATS(this, NAME, 0)
#define TVPI_TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS_ASSIGN_CTOR(NAME)              \
  CRAB_DOMAIN_SCOPED_STATS(&o, NAME, 0)

template <class Number, class VariableName,
          class DBMParams = DBM_impl::DefaultParams<Number>,
          class DomainParams = TVPISplitDBMDefaultParams>
class tvpi_split_dbm_domain final
    : public abstract_domain_api<tvpi_split_dbm_domain<
          Number, VariableName, DBMParams, DomainParams>> {
  using DBM_t =
      tvpi_split_dbm_domain<Number, VariableName, DBMParams, DomainParams>;
  using abstract_domain_t = abstract_domain_api<DBM_t>;

public:
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
  using number_t = Number;
  using varname_t = VariableName;
  using constraint_kind_t = typename linear_constraint_t::kind_t;

private:
  using bound_t = ikos::bound<number_t>;
  using Wt = typename DBMParams::Wt;
  using graph_t = typename DBMParams::graph_t;
  using ntow = DBM_impl::NtoW<number_t, Wt>;
  using ntoc = tvpi_utils::NtoC<number_t>;
  using wtoc = tvpi_utils::WtoC<Wt>;
  using tvpi_op = tvpi_utils::tvpi_op<number_t, variable_t>;
  using vert_id = typename graph_t::vert_id;
  using wt_ref_t = typename graph_t::wt_ref_t;
  using key_t = tvpi_utils::var_coeff_key<variable_t>;
#ifdef USE_FLAT_MAP
  using vert_map_t = boost::container::flat_map<key_t, vert_id>;
#else
  using vert_map_t = std::unordered_map<key_t, vert_id>;
#endif
  using vmap_elt_t = typename vert_map_t::value_type;
  using rev_map_t = std::vector<boost::optional<key_t>>;
  using GrOps = GraphOps<graph_t>;
  using GrPerm = GraphPerm<graph_t>;
  using edge_vector = typename GrOps::edge_vector;
  // < <ax, by>, k> == ax - by <= k.
  using diffcst_t = std::pair<std::pair<key_t, key_t>, Wt>;
  using vert_set_t = std::unordered_set<vert_id>;

protected:
  //================
  // Domain data
  //================
  // GKG: ranges are now maintained in the graph

  // ---- mapping ----
  vert_map_t vert_map; // Mapping from variables to vertices
  rev_map_t rev_map;   // Reverse Mapping from vertices to variables

  // ---- direct graph ----
  graph_t g; // Underlying relation graph representing UTVPI constraints.
             // This is the original, untransformed graph.
             // When a reformulated is needed (e.g., Dijkstra),
             // we convert it by applying the potential function.
  //  Q: Is there any better way to keep direct graph?
  //  The answer is yes, our current implementation is not optimal with respect
  //  to resultant for edges required aligning coefficients.
  //  In short, such resultant computation beyond to what original direct graph
  //  provides.
  //  Now the implementation is looking back to the original form, compute new
  //  form and then convert to new edge to the graph. All lookup and backward
  //  lookup requires map searching.
  //  The average cost for one edge computation is O(1).

  // ---- potential function ----
  // In the paper it computes the bound between two variables
  // slack(y, x) = \pi(x) + k - \pi(y) for y - x <= k
  // for all x and y, slack(y, x) >= 0 iff graph is sat;
  // Otherwise, the value is bot.
  // the pi maps each vertex to a weight
  std::vector<Wt> potential; // Stored potential for the vertex

  // ---- others ----
  // A set of vertices that has this property:
  //  for all x in the set, if there exists an edge x -> y with weight w1 ==> w2
  //  that is weakened during iterations (w1 <= w2).
  // Then the edge has been removed during widening.
  // However, since we remove the edge, the graph becomes not closed.
  // So we use this vertex set to perform close-after-widen operation to restore
  // closure.
  vert_set_t unstable;
  bool _is_bottom;

  // ============================================================
  // Vertex Helpers
  // ============================================================
#pragma region Vertex
  class Wt_max {
  public:
    Wt_max() {}
    Wt apply(const Wt &x, const Wt &y) { return std::max(x, y); }
    bool default_is_absorbing() { return true; }
  };

  class Wt_min {
  public:
    Wt_min() {}
    Wt apply(const Wt &x, const Wt &y) { return std::min(x, y); }
    bool default_is_absorbing() { return false; }
  };

  vert_id get_vert(key_t k) {
    auto it = vert_map.find(k);
    if (it != vert_map.end())
      return (*it).second;

    Wt w(0);
    // if (coefficient != 1) {
    //   boost::optional<vert_id> v1 = get_vert(v);
    //   if (v1) {
    //     crab::outs() << *v1 << "\n";
    //     w = potential[*v1] * Wt(coefficient);
    //   }
    // }

    vert_id vert(g.new_vertex());
    // Initialize potential and reverse mapping
    assert(vert <= rev_map.size());
    if (vert < rev_map.size()) {
      potential[vert] = w;
      rev_map[vert] = k;
    } else {
      potential.push_back(w);
      rev_map.push_back(k);
    }
    vert_map.insert(vmap_elt_t(k, vert));

    assert(vert != 0);

    return vert;
  }

  vert_id get_vert(const variable_t &v) { return get_vert(v, 1); }

  vert_id get_vert(const variable_t &v, unsigned coefficient) {
    return get_vert(key_t(v, coefficient));
  }

  boost::optional<vert_id> get_vert_opt(const key_t &k) const {
    auto it = vert_map.find(k);
    if (it != vert_map.end()) {
      return (*it).second;
    } else {
      return boost::none;
    }
  }

  boost::optional<vert_id> get_vert_opt(const variable_t &v,
                                        const unsigned &coefficient) const {
    return get_vert_opt(key_t(v, coefficient));
  }

  boost::optional<vert_id> get_vert_opt(const variable_t &v) const {
    return get_vert_opt(v, 1);
  }

#pragma endregion Vertex

  // ============================================================
  // Potential Helpers
  // ============================================================
#pragma region Potential
  template <class G, class P>
  static void check_potential(const G &g, const P &p, unsigned line) {
#ifdef CHECK_POTENTIAL
    for (vert_id v : g.verts()) {
      for (vert_id d : g.succs(v)) {
        if (p[v] + g.edge_val(v, d) - p[d] < Wt(0)) {
          CRAB_ERROR("Invalid potential at line ", line, ":", "pot[", v,
                     "]=", p[v], " ", "pot[", d, "]=", p[d], " ", "edge(", v,
                     ",", d, ")=", g.edge_val(v, d));
        }
      }
    }
#endif
  }

  // TODO: need a check tvpi function to make sure inequalities are satisfied?

  class vert_set_wrap_t {
  public:
    vert_set_wrap_t(const vert_set_t &_vs) : vs(_vs) {}

    bool operator[](vert_id v) const { return vs.find(v) != vs.end(); }
    const vert_set_t &vs;
  };

  // Evaluate the potential value of a variable.
  Wt pot_value(const variable_t &v) {
    auto it = vert_map.find(key_t(v));
    if (it != vert_map.end())
      return potential[(*it).second];
    return ((Wt)0); // default, 0
  }

  // Evaluate an expression under the chosen potentials
  Wt eval_expression(const linear_expression_t &e, bool &overflow) {
    overflow = false;
    Wt v(ntow::convert(e.constant(), overflow));
    if (overflow) {
      return Wt(0);
    }

    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      Wt coef = ntow::convert(it->first, overflow);
      if (overflow) {
        return Wt(0);
      }
      v += (pot_value(it->second) - potential[0]) * coef;
    }
    return v;
  }

  interval_t eval_interval(const linear_expression_t &e) {
    interval_t r = e.constant();
    for (auto it = e.begin(), et = e.end(); it != et; ++it)
      r += it->first * operator[](it->second);
    return r;
  }

  interval_t compute_residual(const linear_expression_t &e,
                              const variable_t &pivot) {
    interval_t residual(-e.constant());
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const variable_t &v = (*it).second;
      if (v.index() != pivot.index()) {
        residual = residual - (interval_t((*it).first) * this->operator[](v));
      }
    }
    return residual;
  }
#pragma endregion Potential

  // ============================================================
  // Assignment Helpers
  // ============================================================
#pragma region Assignment
  /**
   *  Turn an assignment into a set of difference constraints.
   *
   *  Given x := a*y + b*z + k, where a,b >= 0, we generate the
   *  difference constraints:
   *
   *  if extract_upper_bounds
   *     x - y <= ub((a-1)*y + b*z + k)
   *     x - z <= ub(a*y + (b-1)*z + k)
   *  else
   *     y - x <= lb((a-1)*y + b*z + k)
   *     z - x <= lb(a*y + (b-1)*z + k)
   *
   * additionally, for tvpi:
   *  if extract_upper_bounds
   *     x - ay <= ub(b*z + k)
   *     x - bz <= ub(a*y + k)
   *  else
   *     ay - x <= lb(b*z + k)
   *     bz - x <= lb(a*y + k)
   **/
  void diffcsts_of_assign(const variable_t &x, const linear_expression_t &exp,
                          /* if true then process the upper
                             bounds, else the lower bounds */
                          bool extract_upper_bounds,
                          /* output: foreach {cy, k} \in diff_csts we have
                             the difference constraint UTVPI(cy, k) */
                          std::vector<std::pair<key_t, Wt>> &diff_csts) {

    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".diffcsts_of_assign");
    boost::optional<key_t> unbounded_var;
    std::vector<std::pair<key_t, Wt>> terms;
    bool overflow;

    Wt residual(ntow::convert(exp.constant(), overflow));
    if (overflow) {
      return;
    }

    for (auto it = exp.begin(), et = exp.end(); it != et; ++it) {
      const variable_t &y = (*it).second;
      const number_t &nc = (*it).first;
      Wt coeff(ntow::convert(nc, overflow));
      unsigned c = ntoc::convert(nc);
      if (overflow) {
        continue;
      }
      if (coeff < Wt(0)) { // e.g. x := -3y + ...
        bound_t y_val =
            (extract_upper_bounds ? operator[](y).lb() : operator[](y).ub());

        if (y_val.is_infinite()) {
          return;
        }
        residual += ntow::convert(*(y_val.number()), overflow) * coeff;
        if (overflow) {
          continue;
        }
      } else { // e.g. x := 2y + ...
        bound_t y_val =
            (extract_upper_bounds ? operator[](y).ub() : operator[](y).lb());

        if (y_val.is_infinite()) {
          if (unbounded_var) {
            return;
          }
          unbounded_var = key_t(y, c);
        } else {
          Wt ymax(ntow::convert(*(y_val.number()), overflow));
          if (overflow) {
            continue;
          }
          residual += ymax * coeff;
          terms.push_back({key_t(y, 1), ymax});
          if (tvpi_utils::find(crab_domain_params_man::get().coefficients(),
                               c)) {
            terms.push_back({key_t(y, c), ymax});
            // TODO: need to pay attention on ax - cy * a <= a * k
          }
        }
      }
    }

    if (unbounded_var) {
      // There is exactly one unbounded variable
      diff_csts.push_back({*unbounded_var, residual});
    } else {
      for (auto &p : terms) {
        if (p.first.coeff() == 1) {
          diff_csts.push_back({p.first, residual - p.second});
        } else {
          diff_csts.push_back({p.first, residual - p.second * p.first.coeff()});
        }
      }
    }
    CRAB_LOG(
        "tvpi-dbm-assign2", crab::outs()
                                << "new TVPI constraints for " << x << ":\n";
        for (auto &kv
             : diff_csts) {
          auto &cy = kv.first;
          auto k = kv.second;
          if (extract_upper_bounds) {
            crab::outs() << x << " - " << cy << " <= " << k << "\n";
          } else {
            crab::outs() << cy << " - " << x << " <= " << -k << "\n";
          }
        } crab::outs()
        << "\n";);
  }

  // Turn an assignment into a set of difference constraints.
  void diffcsts_of_assign(const variable_t &x, const linear_expression_t &exp,
                          std::vector<std::pair<key_t, Wt>> &lb,
                          std::vector<std::pair<key_t, Wt>> &ub) {
    diffcsts_of_assign(x, exp, true, ub);
    diffcsts_of_assign(x, exp, false, lb);
  }

  /**
   * Turn a linear inequality into a set of difference
   * constraints.
   **/
  void diffcsts_of_lin_leq(const linear_expression_t &exp,
                           /* difference constraints */
                           std::vector<diffcst_t> &csts,
                           /* x >= lb for each {x,lb} in lbs */
                           std::vector<std::pair<variable_t, Wt>> &lbs,
                           /* x <= ub for each {x,ub} in ubs */
                           std::vector<std::pair<variable_t, Wt>> &ubs) const {

    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".diffcsts_of_lin_leq");
    Wt unbounded_lbcoeff;
    Wt unbounded_ubcoeff;
    boost::optional<variable_t> unbounded_lbvar; // var with unknow lower bound
    boost::optional<variable_t> unbounded_ubvar; // var with unknow upper bound
    bool underflow, overflow;

    Wt exp_ub = -(ntow::convert(exp.constant(), overflow));
    if (overflow) {
      return;
    }

    // temporary hack
    ntow::convert(exp.constant() - 1, underflow);
    if (underflow) {
      // We don't like MIN either because the code will compute
      // minus MIN and it will silently overflow.
      return;
    }

    std::vector<std::pair<std::pair<Wt, variable_t>, Wt>> pos_terms, neg_terms;
    for (auto it = exp.begin(), et = exp.end(); it != et; ++it) {
      const variable_t &y = (*it).second;
      const number_t &nc = (*it).first;
      Wt coeff(ntow::convert(nc, overflow));
      if (overflow) {
        continue;
      }
      if (coeff > Wt(0)) {
        bound_t y_lb = at(y).lb();
        if (y_lb.is_infinite()) {
          if (unbounded_lbvar) {
            return;
          }
          unbounded_lbvar = y;
          unbounded_lbcoeff = coeff;
        } else {
          Wt ymin(ntow::convert(*(y_lb.number()), overflow));
          if (overflow) {
            continue;
          }
          exp_ub -= ymin * coeff;
          pos_terms.push_back({{coeff, y}, ymin});
        }
      } else {
        bound_t y_ub = at(y).ub();
        if (y_ub.is_infinite()) {
          if (unbounded_ubvar) {
            return;
          }
          unbounded_ubvar = y;
          unbounded_ubcoeff = -coeff;
        } else {
          Wt ymax(ntow::convert(*(y_ub.number()), overflow));
          if (overflow) {
            continue;
          }
          exp_ub -= ymax * coeff;
          neg_terms.push_back({{-coeff, y}, ymax});
        }
      }
    }

    if (unbounded_lbvar) {
      const variable_t &x = *unbounded_lbvar;
      unsigned a = wtoc::convert(unbounded_lbcoeff);
      if (unbounded_ubvar) {
        unsigned b = wtoc::convert(unbounded_ubcoeff);
        if (!(tvpi_utils::find(crab_domain_params_man::get().coefficients(),
                               a) ||
              a == 1) ||
            !(tvpi_utils::find(crab_domain_params_man::get().coefficients(),
                               b) ||
              b == 1)) {
          return;
        }
        const variable_t &y = *unbounded_ubvar;
        csts.push_back({{key_t(x, a), key_t(y, b)}, exp_ub});
      } else {
        if (tvpi_utils::find(crab_domain_params_man::get().coefficients(), a) ||
            a == 1) {
          for (auto &p : neg_terms) {
            csts.push_back(
                {{key_t(x, a), key_t(p.first.second)}, exp_ub - p.second});
            unsigned c = ntoc::convert(p.first.first);
            if (tvpi_utils::find(crab_domain_params_man::get().coefficients(),
                                 c)) {
              csts.push_back({{key_t(x, a), key_t(p.first.second, c)},
                              exp_ub - p.second * c});
            }
          }
        }
        // Add bounds for x
        ubs.push_back({x, exp_ub / unbounded_lbcoeff});
      }
    } else {
      if (unbounded_ubvar) {
        const variable_t &y = *unbounded_ubvar;
        unsigned b = wtoc::convert(unbounded_ubcoeff);
        if (tvpi_utils::find(crab_domain_params_man::get().coefficients(), b) ||
            b == 1) {
          for (auto &p : pos_terms) {
            csts.push_back(
                {{key_t(p.first.second), key_t(y, b)}, exp_ub + p.second});
            unsigned c = ntoc::convert(p.first.first);
            if (tvpi_utils::find(crab_domain_params_man::get().coefficients(),
                                 c)) {
              csts.push_back({{key_t(p.first.second, c), key_t(y, b)},
                              exp_ub + p.second * c});
            }
          }
        }
        // Add bounds for y
        lbs.push_back({y, -exp_ub / unbounded_ubcoeff});
      } else {
        for (auto &pl : neg_terms) {
          for (auto &pu : pos_terms) {
            csts.push_back({{key_t(pu.first.second), key_t(pl.first.second)},
                            exp_ub - pl.second + pu.second});
            unsigned c1 = ntoc::convert(pu.first.first);
            unsigned c2 = ntoc::convert(pl.first.first);
            if (tvpi_utils::find(crab_domain_params_man::get().coefficients(),
                                 c1)) {
              csts.push_back(
                  {{key_t(pu.first.second, c1), key_t(pl.first.second)},
                   exp_ub - pl.second + pu.second * c1});
              if (tvpi_utils::find(crab_domain_params_man::get().coefficients(),
                                   c2)) {
                csts.push_back(
                    {{key_t(pu.first.second, c1), key_t(pl.first.second, c2)},
                     exp_ub - pl.second * c2 + pu.second * c1});
              }
            } else if (tvpi_utils::find(
                           crab_domain_params_man::get().coefficients(), c2)) {
              csts.push_back(
                  {{key_t(pu.first.second), key_t(pl.first.second, c2)},
                   exp_ub - pl.second * c2 + pu.second});
            }
          }
        }
        for (auto &pl : neg_terms) {
          lbs.push_back(
              {pl.first.second, -exp_ub / pl.first.first + pl.second});
        }
        for (auto &pu : pos_terms) {
          ubs.push_back({pu.first.second, exp_ub / pu.first.first + pu.second});
        }
      }
    }
    CRAB_LOG(
        "tvpi-dbm-+=2", if (!csts.empty()) {
          crab::outs() << "new TVPI constraints:\n";
          for (auto &cst : csts) {
            auto &cy = cst.first.first;
            auto &cx = cst.first.second;
            auto k = cst.second;
            crab::outs() << cy << " - " << cx << " <= " << k << "\n";
          }
          crab::outs() << "\n";
        });
  }

  boost::optional<linear_expression_t>
  try_rewrite_linear_expression(const linear_expression_t &e,
                                const unsigned &coefficient) const {
    /**
     *
     * Given c1*x1 + c2*x2 + ... + k and a coefficient c, rewrite into:
     *  c1' = c1 * c, c2' = c2 * c, ..., k' = k * c
     **/
    if (e.is_constant()) {
      return boost::none;
    }
    linear_expression_t res;
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const variable_t &v = (*it).second;
      const number_t &coeff = (*it).first;
      bool neg = coeff < 0;
      number_t abs_coeff = coeff < 0 ? -coeff : coeff;
      if (abs_coeff == number_t(0)) {
        continue;
      }
      number_t new_coeff = abs_coeff * number_t(coefficient);
      boost::optional<vert_id> gv = get_vert_opt(v, ntoc::convert(new_coeff));
      if (gv == boost::none) {
        return boost::none;
      }
      if (neg) {
        res = res - new_coeff * v;
      } else {
        res = res + new_coeff * v;
      }
    }
    const number_t k = e.constant();
    if (k != number_t(0)) {
      res = res + k * number_t(coefficient);
    }
    return res;
  }

  boost::optional<linear_constraint_t>
  try_rewrite_linear_constraint(const linear_constraint_t &cst,
                                const unsigned &coefficient) const {
    auto ce = try_rewrite_linear_expression(cst.expression(), coefficient);
    if (ce == boost::none) {
      return boost::none;
    } else {
      return linear_constraint_t(*ce, cst.kind());
    }
  }

  bool add_linear_leq(const linear_expression_t &exp) {
    // given a linear expression of the form: a*x + b*y + c*z <= k
    // convert into a*x + b*y + c*z - k <= 0
    CRAB_LOG("tvpi-dbm-+=", linear_expression_t exp_tmp(exp);
             crab::outs() << "Adding: " << exp_tmp << "<= 0"
                          << "\n");
    std::vector<std::pair<variable_t, Wt>> lbs, ubs;
    std::vector<diffcst_t> csts;
    // same procedure as assign, try computes TVPI constraints
    diffcsts_of_lin_leq(exp, csts, lbs, ubs);

    check_potential(g, potential, __LINE__);

    Wt_min min_op;
    wt_ref_t w;
    for (auto &p : lbs) {
      CRAB_LOG("tvpi-dbm-+=", crab::outs()
                                  << p.first << ">=" << p.second << "\n");
      auto c_list = crab_domain_params_man::get().coefficients();
      c_list.push_back(1); // add 1 to coefficients
      for (auto &c : c_list) {
        boost::optional<vert_id> v_opt =
            c == 1 ? get_vert(p.first, c) : get_vert_opt(p.first, c);
        if (v_opt == boost::none) {
          continue;
        }
        vert_id v = *v_opt;
        Wt p_val = p.second * Wt(c);
        CRAB_LOG("tvpi-dbm-+=2",
                 crab::outs() << key_t(p.first, c) << ">=" << p_val << "\n");
        if (g.lookup(v, 0, w) && w.get() <= -p_val)
          continue;
        g.set_edge(v, -p_val, 0);

        if (!repair_potential(v, 0)) {
          set_to_bottom();
          return false;
        }
        check_potential(g, potential, __LINE__);
      }
    }
    for (auto &p : ubs) {
      CRAB_LOG("tvpi-dbm-+=", crab::outs()
                                  << p.first << "<=" << p.second << "\n");
      vert_id v = get_vert(p.first);
      if (g.lookup(0, v, w) && w.get() <= p.second)
        continue;
      g.set_edge(0, p.second, v);
      if (!repair_potential(0, v)) {
        set_to_bottom();
        return false;
      }
      check_potential(g, potential, __LINE__);
    }

    for (auto &diff : csts) {
      CRAB_LOG("tvpi-dbm-+=", crab::outs() << diff.first.first << "-"
                                           << diff.first.second
                                           << "<=" << diff.second << "\n");
      vert_id src = get_vert(diff.first.second);
      vert_id dest = get_vert(diff.first.first);

      // Check if the edge (src,dest) via bounds already exists
      wt_ref_t w1, w2;
      if (g.lookup(src, 0, w1) && g.lookup(0, dest, w2) &&
          (w1.get() + w2.get()) <= diff.second) {
        continue;
      }

      g.update_edge(src, diff.second, dest, min_op);
      if (!repair_potential(src, dest)) {
        set_to_bottom();
        return false;
      }
      check_potential(g, potential, __LINE__);
      close_over_edge(src, dest);
      reduce_tvpi_edge(src, dest, boost::none, vert_map);
      check_potential(g, potential, __LINE__);
    }
    // Collect bounds
    // GKG: Now done in close_over_edge

    edge_vector delta;
    GrOps::close_after_assign(g, potential, 0, delta);
    GrOps::apply_delta(g, delta);

    check_potential(g, potential, __LINE__);
    return true;
  }

  // x != n
  bool add_univar_disequation(const variable_t &x, number_t n) {
    CRAB_LOG("tvpi-dbm-+=", crab::outs() << x << "!=" << n << "\n");
    bool overflow;
    interval_t i = get_interval(x);
    interval_t ni(n);
    interval_t new_i =
        ikos::linear_interval_solver_impl::trim_interval<interval_t>(i, ni);
    if (new_i.is_bottom()) {
      set_to_bottom();
      return false;
    } else if (!new_i.is_top() && (new_i <= i)) {
      vert_id v = get_vert(x);
      wt_ref_t w;
      Wt_min min_op;
      if (new_i.lb().is_finite()) {
        // strenghten lb
        Wt lb_val = ntow::convert(-(*(new_i.lb().number())), overflow);
        if (overflow) {
          return true;
        }

        if (g.lookup(v, 0, w) && lb_val < w.get()) {
          g.set_edge(v, lb_val, 0);
          if (!repair_potential(v, 0)) {
            set_to_bottom();
            return false;
          }
          check_potential(g, potential, __LINE__);
          // Update other bounds
          for (auto e : g.e_preds(v)) {
            if (e.vert == 0)
              continue;
            g.update_edge(e.vert, e.val + lb_val, 0, min_op);
            if (!repair_potential(e.vert, 0)) {
              set_to_bottom();
              return false;
            }
            check_potential(g, potential, __LINE__);
          }
        }
      }
      if (new_i.ub().is_finite()) {
        // strengthen ub
        Wt ub_val = ntow::convert(*(new_i.ub().number()), overflow);
        if (overflow) {
          return true;
        }

        if (g.lookup(0, v, w) && (ub_val < w.get())) {
          g.set_edge(0, ub_val, v);
          if (!repair_potential(0, v)) {
            set_to_bottom();
            return false;
          }
          check_potential(g, potential, __LINE__);
          // Update other bounds
          for (auto e : g.e_succs(v)) {
            if (e.vert == 0)
              continue;
            g.update_edge(0, e.val + ub_val, e.vert, min_op);
            if (!repair_potential(0, e.vert)) {
              set_to_bottom();
              return false;
            }
            check_potential(g, potential, __LINE__);
          }
        }
      }
    }
    return true;
  }

  void add_disequation(const linear_expression_t &e) {
    // XXX: similar precision as the interval domain

    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const variable_t &pivot = (*it).second;
      interval_t i = compute_residual(e, pivot) / interval_t((*it).first);
      if (auto k = i.singleton()) {
        if (!add_univar_disequation(pivot, *k)) {
          // set_to_bottom() was already called
          return;
        }
      }
    }
  }
#pragma endregion Assignment

  // ============================================================
  // Interval Helpers
  // ============================================================
#pragma region Interval
  interval_t get_interval(const variable_t &x) const {
    return get_interval(vert_map, g, x);
  }

  interval_t get_interval(const vert_map_t &m, const graph_t &g,
                          const variable_t &x) const {
    key_t k(x);
    auto it = m.find(k);
    if (it == m.end()) {
      return interval_t::top();
    }
    vert_id v = (*it).second;
    interval_t x_out = interval_t(
        g.elem(v, 0) ? -number_t(g.edge_val(v, 0)) : bound_t::minus_infinity(),
        g.elem(0, v) ? number_t(g.edge_val(0, v)) : bound_t::plus_infinity());
    return x_out;
  }

  // Restore potential after an edge addition
  // Source: Fast and Flexible Difference Constraint Propagation for DPLL(T)
  // Fig. 1. algorithm
  bool repair_potential(vert_id src, vert_id dest) {
    return GrOps::repair_potential(g, potential, src, dest);
  }
#pragma endregion Interval

  // ============================================================
  // Graph Helpers
  // ============================================================
#pragma region Graph
  // Restore closure after a single edge addition
  void close_over_edge(vert_id ii, vert_id jj) {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".close_edge");
    assert(ii != 0 && jj != 0);

    Wt_min min_op;
    wt_ref_t w;
    SubGraph<graph_t> g_excl(
        g, 0); // another way to interpret graph but avoiding v0
    Wt c = g_excl.edge_val(ii, jj);

    std::vector<std::pair<vert_id, Wt>> src_dec;
    // we add in delta so that we don't invalidate graph iterators
    edge_vector delta;
    for (auto edge : g_excl.e_preds(ii)) {
      vert_id se = edge.vert;
      Wt wt_sij = edge.val + c;
      assert(g_excl.succs(se).begin() != g_excl.succs(se).end());
      if (se != jj) {
        if (g_excl.lookup(se, jj, w)) {
          if (w.get() <= wt_sij) {
            continue;
          }
          // REVISIT(PERFORMANCE): extra call to lookup
          g.set_edge(se, wt_sij, jj);
        } else {
          delta.push_back({{se, jj}, wt_sij});
        }
        src_dec.push_back(std::make_pair(se, edge.val));
      }
    }

    GrOps::apply_delta(g, delta);
    delta.clear();
    std::vector<std::pair<vert_id, Wt>> dest_dec;
    for (auto edge : g_excl.e_succs(jj)) {
      vert_id de = edge.vert;
      Wt wt_ijd = edge.val + c;
      if (de != ii) {
        if (g_excl.lookup(ii, de, w)) {
          if (w.get() <= wt_ijd) {
            continue;
          }
          // REVISIT(PERFORMANCE): extra call to lookup
          g.set_edge(ii, wt_ijd, de);
        } else {
          delta.push_back({{ii, de}, {wt_ijd}});
        }
        dest_dec.push_back(std::make_pair(de, edge.val));
      }
    }
    GrOps::apply_delta(g, delta);

    for (auto s_p : src_dec) {
      vert_id se = s_p.first;
      Wt wt_sij = c + s_p.second;
      for (auto d_p : dest_dec) {
        vert_id de = d_p.first;
        Wt wt_sijd = wt_sij + d_p.second;
        if (g.lookup(se, de, w)) {
          if (w.get() <= wt_sijd) {
            continue;
          }
          // REVISIT(PERFORMANCE): extra call to lookup
          g.set_edge(se, wt_sijd, de);
        } else {
          g.add_edge(se, wt_sijd, de);
        }
      }
    }

    // Closure is now updated.
  }

  // Restore TVPI closure based on the current edge
  void reduce_tvpi_edge(vert_id ii, vert_id jj, boost::optional<vert_id> new_v,
                        vert_map_t &tmp_vert_map) {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".reduce_tvpi_edge");
    assert(ii != 0 && jj != 0);
    // get all vertices relate to ii
    std::pair<std::vector<key_t>, bool>
        xs; // all vertices related to variable x, so dx
    std::pair<std::vector<key_t>, bool>
        ys; // all vertices related to variable y, so dy
    SubGraph<graph_t> g_excl(
        g, 0); // another way to interpret graph but avoiding v0
    Wt c = g_excl.edge_val(ii, jj);
    CRAB_LOG("tvpi-dbm-tvpi", crab::outs() << "\n===\nBefore:\n";
             tvpi_utils::print_map(crab::outs(), vert_map);
             tvpi_utils::print_map(crab::outs(), tmp_vert_map);
             tvpi_utils::print_vector(crab::outs(), rev_map);
             crab::outs() << "\n"; crab::outs() << *this << "\n";
             /*print_details(crab::outs());*/);

    std::set<vert_id> new_verts;

    auto is_vert_new = [this, &tmp_vert_map](vert_id v) -> bool {
      if (v >= rev_map.size()) {
        CRAB_ERROR(LOCATION_STRING, " Invalid vertex id: ", v,
                   " out of range.");
      }
      if (v == 0) {
        return false;
      }
      if (rev_map[v] == boost::none) {
        CRAB_ERROR(LOCATION_STRING, "vertex ", v, " has no reverse mapping.");
      }
      return tmp_vert_map.find(rev_map[v].get()) != tmp_vert_map.end();
    };

    auto find_vert = [this](key_t k, vert_map_t &vert_map) -> vert_id {
      auto it = vert_map.find(k);
      if (it != vert_map.end()) {
        return (*it).second;
      }

      vert_id v = g.new_vertex(); // graph tracks empty vertex id
      assert(v <= rev_map.size());
      if (v == rev_map.size()) {
        rev_map.push_back(k);
        potential.push_back(Wt(0));
      } else {
        potential[v] = Wt(0);
        rev_map[v] = k;
      }
      vert_map.insert(vmap_elt_t(k, v)); // insert map for quick access
      return v;
    };

    auto find_keys_for_var = [&new_verts](const variable_t &v,
                                          std::vector<key_t> &keys, bool is_new,
                                          const vert_map_t &vert_map) {
      for (const auto &kv : vert_map) {
        const key_t &k = kv.first;
        if (k.var() == v) {
          keys.push_back(k);
          if (is_new) {
            new_verts.insert(kv.second);
          }
        }
      }
    };

    auto find_keys_with_same_var =
        [&tmp_vert_map, this,
         &find_keys_for_var](bool is_new, const variable_t &var,
                             std::pair<std::vector<key_t>, bool> &res) {
          res.second = is_new;
          find_keys_for_var(var, res.first, res.second,
                            res.second ? tmp_vert_map : vert_map);
        };

    bool overflow = false;
    Wt_min min_op;
    auto coeffs = crab_domain_params_man::get().coefficients();

    bool is_new_jj = new_v && *new_v == jj;
    bool is_new_ii = new_v && *new_v == ii;
    key_t ax =
        (rev_map[jj] ? key_t(*rev_map[jj]) : key_t((*rev_map[ii]).var(), 0));
    key_t by =
        (rev_map[ii] ? key_t(*rev_map[ii]) : key_t((*rev_map[jj]).var(), 0));
    // crab::outs() << "\nNew: " << ax << " - " << by << " <= " << c << "\n";
    // crab::outs() << "is_new_ax: " << is_new_jj << ", is_new_by: " <<
    // is_new_ii << "\n";

    variable_t x = ax.var();
    if (ax.coeff() > 0) {
      find_keys_with_same_var(is_new_jj, ax.var(),
                              xs); // find all keys relate to x
    }
    variable_t y = by.var();
    if (by.coeff() > 0) {
      find_keys_with_same_var(is_new_ii, by.var(),
                              ys); // find all keys relate to y
    }

    // A generic helper to process incremental saturation
    auto process_edge_loop = [&](bool forward) {
      // we add in delta so that we don't invalidate graph iterators
      edge_vector delta;
      std::vector<std::pair<vert_id, vert_id>> need_to_close;
      std::vector<std::pair<vert_id, vert_id>> need_to_close2;
      const std::vector<key_t> &keys = forward ? xs.first : ys.first;
      // layer 1: process edges by eliminating x or y
      for (auto &k : keys) {
        vert_id v = find_vert(
            k, ((forward ? xs.second : ys.second) ? tmp_vert_map : vert_map));

        for (vert_id vez : (forward ? g.succs(v) : g.preds(v))) {
          if (rev_map[vez] == boost::none) {
            continue;
          }
          key_t ez = *rev_map[vez];
          if (ez.var() == k.var()) {
            continue;
          }
          if (forward ? ez.var() == by.var() : ez.var() == ax.var()) {
            continue;
          }
          if (forward ? k.coeff() == ax.coeff() : k.coeff() == by.coeff()) {
            continue;
          }
          Wt f = forward ? g_excl.edge_val(v, vez) : g_excl.edge_val(vez, v);
          CRAB_LOG("tvpi-dbm-tvpi",
                   forward
                       ? crab::outs()
                             << "resultant(" << ez << "-" << k << " <= " << f
                             << ", " << ax << "-" << by << " <= " << c
                             << "), eliminating " << ax.var() << "\n"
                       : crab::outs()
                             << "resultant(" << ax << "-" << by << " <= " << c
                             << ", " << k << "-" << ez << " <= " << f
                             << "), eliminating " << by.var() << "\n");

          auto [new_a, new_b, new_c] =
              forward
                  ? tvpi_op::resultant(ez.coeff(), ez.var(), k.coeff(), k.var(),
                                       f, ax.coeff(), by.coeff(), by.var(), c)
                  : tvpi_op::resultant(ax.coeff(), ax.var(), by.coeff(),
                                       by.var(), c, k.coeff(), ez.coeff(),
                                       ez.var(), f);
          bool skip = false;
          key_t new_key_src =
              forward ? key_t(by.var(), new_b) : key_t(ez.var(), new_b);
          key_t new_key_dest =
              forward ? key_t(ez.var(), new_a) : key_t(ax.var(), new_a);
          Wt t(ntow::convert(new_c, overflow));
          if (overflow) {
            skip = true;
          } else if (new_a == 0 && new_b == 0) {
            skip = true;
          } else if (!((new_a == 1 && new_b != 1 &&
                        tvpi_utils::find(coeffs, new_b)) ||
                       (new_b == 1 && new_a != 1 &&
                        tvpi_utils::find(coeffs, new_a)))) {
            skip = true;
          }
          if (!skip) {
            wt_ref_t w;
            bool src_is_new = forward ? is_new_ii : is_vert_new(vez);
            vert_id src = new_b == 0
                              ? 0
                              : find_vert(new_key_src,
                                          src_is_new ? tmp_vert_map : vert_map);
            bool dest_is_new = forward ? is_vert_new(vez) : is_new_jj;
            vert_id dest =
                new_a == 0 ? 0
                           : find_vert(new_key_dest,
                                       dest_is_new ? tmp_vert_map : vert_map);
            if (g.lookup(src, dest, w) && w.get() <= t) {
              skip = true;
            }
            if (!skip) {
              delta.push_back({{src, dest}, {t}});
              // need_to_close.push_back({src, dest});
            }
          }
          CRAB_LOG("tvpi-dbm-tvpi",
                   crab::outs()
                       << "=>>>" << new_key_dest << "-" << new_key_src
                       << "<=" << t << (skip ? ", skip" : ", added") << "\n");
        }
      }
      check_potential(g, potential, __LINE__);
      GrOps::apply_delta(g, delta); // add new edges
      delta.clear();
      check_potential(g, potential, __LINE__);
      // need_to_close.clear();
      if (true) {
        CRAB_LOG("tvpi-dbm-tvpi", crab::outs() << "layer 2\n";);
        // layer 2: process new edges by eliminating remaining y or x
        const std::vector<key_t> &keys2 = forward ? ys.first : xs.first;
        for (auto &k : keys2) {
          vert_id v = find_vert(
              k, ((forward ? ys.second : xs.second) ? tmp_vert_map : vert_map));
          for (vert_id vgw : (forward ? g.preds(v) : g.succs(v))) {
            if (rev_map[vgw] == boost::none) {
              continue;
            }
            key_t gw = *rev_map[vgw];
            if (gw.var() == k.var()) {
              continue;
            }
            if (k == ax && gw == by) {
              continue;
            }
            for (auto sd : need_to_close) {
              vert_id src = sd.first, dst = sd.second;
              key_t key_s = *rev_map[src], key_d = *rev_map[dst];
              if (forward ? gw.var() == key_d.var() : k.var() == key_s.var()) {
                continue;
              }
              if (forward ? k.coeff() == key_s.coeff()
                          : gw.coeff() == key_d.coeff()) {
                continue;
              }
              Wt h =
                  forward ? g_excl.edge_val(vgw, v) : g_excl.edge_val(v, vgw);
              Wt c = g_excl.edge_val(src, dst);
              CRAB_LOG("tvpi-dbm-tvpi",
                       forward ? crab::outs()
                                     << "resultant2(" << key_d << "-" << key_s
                                     << " <= " << c << ", " << k << "-" << gw
                                     << " <= " << h << "), eliminating "
                                     << key_s.var() << "\n"
                               : crab::outs() << "resultant2(" << gw << "-" << k
                                              << " <= " << h << ", " << key_d
                                              << "-" << key_s << " <= " << c
                                              << "), eliminating "
                                              << key_d.var() << "\n");

              auto [new_a, new_b, new_c] =
                  forward
                      ? tvpi_op::resultant(key_d.coeff(), key_d.var(),
                                           key_s.coeff(), key_s.var(), c,
                                           k.coeff(), gw.coeff(), gw.var(), h)
                      : tvpi_op::resultant(gw.coeff(), gw.var(), k.coeff(),
                                           k.var(), h, key_d.coeff(),
                                           key_s.coeff(), key_s.var(), c);
              bool skip = false;
              key_t new_key_src =
                  forward ? key_t(gw.var(), new_b) : key_t(key_s.var(), new_b);
              key_t new_key_dest =
                  forward ? key_t(key_d.var(), new_a) : key_t(gw.var(), new_a);
              Wt t(ntow::convert(new_c, overflow));
              if (overflow) {
                skip = true;
              } else if (new_a == 0 && new_b == 0) {
                skip = true;
              } else if (!((new_a == 1 && new_b != 1 &&
                            tvpi_utils::find(coeffs, new_b)) ||
                           (new_b == 1 && new_a != 1 &&
                            tvpi_utils::find(coeffs, new_a)))) {
                skip = true;
              }
              if (!skip) {
                wt_ref_t w;
                bool src_is_new = forward ? is_vert_new(vgw) : is_vert_new(src);
                vert_id src = new_b == 0 ? 0
                                         : find_vert(new_key_src,
                                                     src_is_new ? tmp_vert_map
                                                                : vert_map);
                bool dest_is_new = forward ? is_vert_new(dst) : is_vert_new(v);
                vert_id dest = new_a == 0 ? 0
                                          : find_vert(new_key_dest,
                                                      dest_is_new ? tmp_vert_map
                                                                  : vert_map);
                if (g.lookup(src, dest, w) && w.get() <= t) {
                  skip = true;
                }
                if (!skip) {
                  delta.push_back({{src, dest}, {t}});
                  // need_to_close2.push_back({src, dest});
                }
              }
              CRAB_LOG("tvpi-dbm-tvpi",
                       crab::outs() << "=>>>" << new_key_dest << "-"
                                    << new_key_src << "<=" << t
                                    << (skip ? ", skip" : ", added") << "\n");
            }
          }
        }
        check_potential(g, potential, __LINE__);
        GrOps::apply_delta(g, delta); // add new edges
        check_potential(g, potential, __LINE__);
        need_to_close2.clear();
      }
      need_to_close.clear();
    };

    // process all vertices relate to y
    if (!(is_new_ii && rev_map[ii] == rev_map[jj])) {
      process_edge_loop(false); // y's predecessors
    }
    // process all vertices relate to x
    if (!(is_new_jj && rev_map[ii] == rev_map[jj])) {
      process_edge_loop(true); // x's successors
    }

    CRAB_LOG("tvpi-dbm-tvpi", crab::outs() << "\n===\nAfter:\n";
             tvpi_utils::print_map(crab::outs(), vert_map);
             tvpi_utils::print_map(crab::outs(), tmp_vert_map);
             crab::outs() << "\n"; crab::outs() << *this << "\n";
             /*print_details(crab::outs());*/);
  }

  // return true if edge from x to y with weight k is unsatisfiable
  bool is_unsat_edge(vert_id x, vert_id y, Wt k) const {
    wt_ref_t w;
    if (g.lookup(y, x, w)) {
      return ((w.get() + k) < Wt(0));
    } else {
      interval_t intv_x = interval_t::top();
      interval_t intv_y = interval_t::top();
      if (g.elem(0, x) || g.elem(x, 0)) {
        intv_x = interval_t(g.elem(x, 0) ? -number_t(g.edge_val(x, 0))
                                         : bound_t::minus_infinity(),
                            g.elem(0, x) ? number_t(g.edge_val(0, x))
                                         : bound_t::plus_infinity());
      }
      if (g.elem(0, y) || g.elem(y, 0)) {
        intv_y = interval_t(g.elem(y, 0) ? -number_t(g.edge_val(y, 0))
                                         : bound_t::minus_infinity(),
                            g.elem(0, y) ? number_t(g.edge_val(0, y))
                                         : bound_t::plus_infinity());
      }
      if (intv_x.is_top() || intv_y.is_top()) {
        return false;
      } else {
        return (!((intv_y - intv_x).lb() <= (number_t)k));
      }
    }
  }

  // return true iff cst is unsatisfiable without modifying the DBM
  bool is_unsat(const linear_constraint_t &cst) const {
    if (is_bottom() || cst.is_contradiction()) {
      return true;
    }

    if (is_top() || cst.is_tautology()) {
      return false;
    }

    std::vector<std::pair<variable_t, Wt>> lbs, ubs;
    std::vector<diffcst_t> diffcsts;

    if (cst.is_inequality()) {
      diffcsts_of_lin_leq(cst.expression(), diffcsts, lbs, ubs);
    } else if (cst.is_strict_inequality()) {
      auto nc =
          ikos::linear_constraint_impl::strict_to_non_strict_inequality(cst);
      if (nc.is_inequality()) {
        diffcsts_of_lin_leq(cst.expression(), diffcsts, lbs, ubs);
      } else {
        // we couldn't convert the strict into a non-strict
        return false;
      }
    } else if (cst.is_equality()) {
      diffcsts_of_lin_leq(cst.expression(), diffcsts, lbs, ubs);
      diffcsts_of_lin_leq(-cst.expression(), diffcsts, lbs, ubs);
    } else if (cst.is_disequation()) {
      CRAB_WARN("disequalities ", cst, " not implemented by ", domain_name(),
                "::is_unsat");
      return false;
    } else {
      return false;
    }

    // check difference constraints
    for (auto &diffcst : diffcsts) {
      key_t x = diffcst.first.first;
      key_t y = diffcst.first.second;
      Wt k = diffcst.second;

      auto vy = get_vert_opt(y);
      auto vx = get_vert_opt(x);
      if (vx && vy && is_unsat_edge(*vy, *vx, k)) {
        return true;
      }
    }

    // check interval constraints
    for (auto &ub : ubs) {
      auto vx = get_vert_opt(ub.first);
      if (vx && is_unsat_edge(0, *vx, ub.second)) {
        return true;
      }
    }
    for (auto &lb : lbs) {
      auto vx = get_vert_opt(lb.first);
      if (vx && is_unsat_edge(*vx, 0, -lb.second)) {
        return true;
      }
    }

    return false;
  }
#pragma endregion Graph

  // ============================================================
  // Domain Operation Helpers
  // ============================================================
#pragma region Helpers
  // Join of gx and gy.
  static graph_t join(GrPerm &gx, GrPerm &gy, unsigned sz,
                      std::vector<Wt> &pot_rx, std::vector<Wt> &pot_ry) {

    // Compute the deferred relations
    graph_t g_ix_ry;
    wt_ref_t ws, wd;
    g_ix_ry.growTo(sz);
    SubGraph<GrPerm> gy_excl(gy, 0);
    for (vert_id s : gy_excl.verts()) {
      for (vert_id d : gy_excl.succs(s)) {
        if (gx.lookup(s, 0, ws) && gx.lookup(0, d, wd)) {
          g_ix_ry.add_edge(s, ws.get() + wd.get(), d);
        }
      }
    }
    // Apply the deferred relations, and re-close.
    bool is_closed;
    graph_t g_rx(GrOps::meet(gx, g_ix_ry, is_closed));
    check_potential(g_rx, pot_rx, __LINE__);
#ifdef JOIN_CLOSE_AFTER_MEET
    // Conjecture: g_rx is closed
    if (!is_closed) {
      edge_vector delta;
      SubGraph<graph_t> g_rx_excl(g_rx, 0);
      GrOps::close_after_meet(g_rx_excl, pot_rx, gx, g_ix_ry, delta);
      GrOps::apply_delta(g_rx, delta);
    }
#endif

    graph_t g_rx_iy;
    g_rx_iy.growTo(sz);
    SubGraph<GrPerm> gx_excl(gx, 0);
    for (vert_id s : gx_excl.verts()) {
      for (vert_id d : gx_excl.succs(s)) {
        // Assumption: gx.mem(s, d) -> gx.edge_val(s, d) <=
        //             ranges[var(s)].ub() - ranges[var(d)].lb()
        // That is, if the relation exists, it's at least as strong as the
        // bounds.
        if (gy.lookup(s, 0, ws) && gy.lookup(0, d, wd))
          g_rx_iy.add_edge(s, ws.get() + wd.get(), d);
      }
    }
    // Similarly, should use a SubGraph view.
    graph_t g_ry(GrOps::meet(gy, g_rx_iy, is_closed));
    check_potential(g_ry, pot_ry, __LINE__);
#ifdef JOIN_CLOSE_AFTER_MEET
    // Conjecture: g_ry is closed
    if (!is_closed) {
      edge_vector delta;
      SubGraph<graph_t> g_ry_excl(g_ry, 0);
      GrOps::close_after_meet(g_ry_excl, pot_ry, gy, g_rx_iy, delta);
      GrOps::apply_delta(g_ry, delta);
    }
#endif

    // We now have the relevant set of relations. Because g_rx
    // and g_ry are closed, the result is also closed.
    Wt_min min_op;
    graph_t join_g(GrOps::join(g_rx, g_ry));

    // Now reapply the missing independent relations.
    // Need to derive vert_ids from lb_up/lb_down, and make sure the vertices
    // exist
    std::vector<vert_id> lb_up;
    std::vector<vert_id> lb_down;
    std::vector<vert_id> ub_up;
    std::vector<vert_id> ub_down;

    wt_ref_t wx, wy;
    for (vert_id v : gx_excl.verts()) {
      if (gx.lookup(0, v, wx) && gy.lookup(0, v, wy)) {
        if (wx.get() < wy.get())
          ub_up.push_back(v);
        if (wy.get() < wx.get())
          ub_down.push_back(v);
      }
      if (gx.lookup(v, 0, wx) && gy.lookup(v, 0, wy)) {
        if (wx.get() < wy.get())
          lb_down.push_back(v);
        if (wy.get() < wx.get())
          lb_up.push_back(v);
      }
    }

    for (vert_id s : lb_up) {
      Wt dx_s = gx.edge_val(s, 0);
      Wt dy_s = gy.edge_val(s, 0);
      for (vert_id d : ub_up) {
        if (s == d)
          continue;
        join_g.update_edge(
            s, std::max(dx_s + gx.edge_val(0, d), dy_s + gy.edge_val(0, d)), d,
            min_op);
      }
    }

    for (vert_id s : lb_down) {
      Wt dx_s = gx.edge_val(s, 0);
      Wt dy_s = gy.edge_val(s, 0);
      for (vert_id d : ub_down) {
        if (s == d)
          continue;
        join_g.update_edge(
            s, std::max(dx_s + gx.edge_val(0, d), dy_s + gy.edge_val(0, d)), d,
            min_op);
      }
    }
    return join_g;
  }

  template <class G1, class G2>
  static graph_t split_widen(G1 &l, G2 &r, std::vector<vert_id> &unstable,
                             const rev_map_t &revmap) {
    assert(l.size() == r.size());
    size_t sz = l.size();
    graph_t g;
    g.growTo(sz);

    Wt_min min_op;
    wt_ref_t wx, wy;

    auto update_edge_widen_g = [&g](vert_id src, vert_id dst, Wt val) {
      wt_ref_t wz;
      if (g.lookup(src, dst, wz)) {
        if (wz.get() > val) {
          g.set_edge(src, val, dst);
          return true;
        }
      } else {
        g.add_edge(src, val, dst);
        return true;
      }
      return false;
    };

    /**
     * Check for stable implicit relationships in r
     **/
    for (auto edge_pred : r.e_preds(0)) {
      vert_id s = edge_pred.vert;
      for (auto edge_succ : r.e_succs(0)) {
        vert_id d = edge_succ.vert;
        if (s == d)
          continue;
        /* for each edge(s,0,d) in r check if exists edge(s,d) in l */
        if (l.lookup(s, d, wx) &&
            ((edge_pred.val + edge_succ.val) <= wx.get())) {
          bool res = update_edge_widen_g(s, d, wx.get());
          if (res) {
            CRAB_LOG("tvpi-dbm-widening", auto vs = revmap[s];
                     auto vd = revmap[d];
                     crab::outs() << "Widening 1: added " << *vd << "-" << *vs
                                  << "<=" << wx.get() << "\n";);
          }
        }
      }
    }

    /**
     * Check for stable explicit relationships in r
     **/
    for (vert_id s : r.verts()) {
      for (auto e : r.e_succs(s)) {
        vert_id d = e.vert;
        /* for each edge(s,d) in r check if exists edge(s,d) in l */
        if (l.lookup(s, d, wx) && e.val <= wx.get()) {
          bool res = update_edge_widen_g(s, d, wx.get());
          if (res) {
            CRAB_LOG(
                "tvpi-dbm-widening", auto vs = revmap[s]; auto vd = revmap[d];
                if (s == 0 && d != 0) {
                  crab::outs() << "Widening 2: added " << *vd
                               << "<=" << wx.get() << "\n";
                } else if (s != 0 && d == 0) {
                  crab::outs() << "Widening 2: added "
                               << "-" << *vs << "<=" << wx.get() << "\n";
                } else {
                  crab::outs() << "Widening 2: added " << *vd << "-" << *vs
                               << "<=" << wx.get() << "\n";
                });
          }
        }
      }

      // Check if this vertex is stable
      for (vert_id d : l.succs(s)) {
        if (!g.elem(s, d)) {
          unstable.push_back(s);
          CRAB_LOG(
              "tvpi-dbm-widening",
              if (s == 0) {
                crab::outs() << "Widening 5: added v0"
                             << " in the normalization queue\n";
              } else {
                auto vs = revmap[s];
                crab::outs() << "Widening 5: added " << *vs
                             << " in the normalization queue\n";
              });
          break;
        }
      }
    }

    // for(vert_id s: r.verts()) {
    //   for(vert_id d : l.succs(s)) {
    //     if(!g.elem(s, d)) {
    //       unstable.push_back(s);
    //       break;
    //     }
    //   }
    // }
    return g;
  }

  bool need_normalization() const {
#ifdef TVPIDBM_NO_NORMALIZE
    return false;
#endif
    return unstable.size() > 0;
  }

  // dbm is already normalized
  linear_constraint_system_t
  to_linear_constraint_system(const DBM_t &dbm) const {
    linear_constraint_system_t csts;

    if (dbm.is_bottom()) {
      csts += linear_constraint_t::get_false();
      return csts;
    }

    // Extract all the edges
    SubGraph<graph_t> g_excl(const_cast<graph_t &>(dbm.g), 0);
    for (vert_id v : g_excl.verts()) {
      if (!dbm.rev_map[v])
        continue;
      if (dbm.g.elem(v, 0)) {
        variable_t vv = dbm.rev_map[v]->var();
        Wt c = dbm.g.edge_val(v, 0);
        csts += linear_constraint_t(linear_expression_t(vv) >= -number_t(c));
      }
      if (dbm.g.elem(0, v)) {
        variable_t vv = dbm.rev_map[v]->var();
        Wt c = dbm.g.edge_val(0, v);
        csts += linear_constraint_t(linear_expression_t(vv) <= number_t(c));
      }
    }

    for (vert_id s : g_excl.verts()) {
      if (!dbm.rev_map[s])
        continue;
      variable_t vs = dbm.rev_map[s]->var();
      for (vert_id d : g_excl.succs(s)) {
        if (!dbm.rev_map[d])
          continue;
        variable_t vd = dbm.rev_map[d]->var();
        csts += linear_constraint_t(vd - vs <= number_t(g_excl.edge_val(s, d)));
      }
    }
    return csts;
  }

  // dbm is already normalized
  void write(crab_os &o, const DBM_t &dbm) const {
#if 0
    o << "edges={";
    for(vert_id v : dbm.g.verts()) {
      for(vert_id d : dbm.g.succs(v)) {
	if(!dbm.rev_map[v] || !dbm.rev_map[d]) {
	  CRAB_WARN("Edge incident to un-mapped vertex.");
	  continue;
	}
	variable_t vv = *dbm.rev_map[v];
	variable_t vd = *dbm.rev_map[d];
	o << "(" << vv << "," << vd << ":"
	  << dbm.g.edge_val(v,d) << ")";
      }
    }
    o << "}";
    crab::outs() << "rev_map={";
    for(unsigned i=0, e = dbm.rev_map.size(); i!=e; i++) {
      if (dbm.rev_map[i]) {
	variable_t vi = *dbm.rev_map[i];
	crab::outs() << vi << "(" << i << ");";
      }
    }
    crab::outs() << "}\n";
#endif
    if (is_bottom()) {
      o << "_|_";
      return;
    } else if (is_top()) {
      o << "{}";
      return;
    } else {
      // Intervals
      bool first = true;
      o << "{";
      // Extract all the edges
      SubGraph<graph_t> g_excl(const_cast<graph_t &>(dbm.g), 0);
      for (vert_id v : g_excl.verts()) {
        if (!dbm.rev_map[v])
          continue;
        if (!dbm.g.elem(0, v) && !dbm.g.elem(v, 0))
          continue;
        interval_t v_out =
            interval_t(dbm.g.elem(v, 0) ? -number_t(dbm.g.edge_val(v, 0))
                                        : bound_t::minus_infinity(),
                       dbm.g.elem(0, v) ? number_t(dbm.g.edge_val(0, v))
                                        : bound_t::plus_infinity());
        if (first)
          first = false;
        else
          o << ", ";
        variable_t vv = dbm.rev_map[v]->var();
        unsigned cc = dbm.rev_map[v]->coeff();
        o << (cc != 1 ? std::to_string(cc) : "") << vv << " -> " << v_out;
      }

      for (vert_id s : g_excl.verts()) {
        if (!dbm.rev_map[s])
          continue;
        variable_t vs = dbm.rev_map[s]->var();
        unsigned cs = dbm.rev_map[s]->coeff();
        for (vert_id d : g_excl.succs(s)) {
          if (!dbm.rev_map[d])
            continue;
          variable_t vd = dbm.rev_map[d]->var();
          unsigned cd = dbm.rev_map[d]->coeff();
          if (first)
            first = false;
          else
            o << ", ";
          o << (cd != 1 ? std::to_string(cd) : "") << vd << "-"
            << (cs != 1 ? std::to_string(cs) : "") << vs
            << "<=" << dbm.g.edge_val(s, d);
        }
      }
      o << "}";
    }
  }

  tvpi_split_dbm_domain(vert_map_t &&_vert_map, rev_map_t &&_rev_map,
                        graph_t &&_g, std::vector<Wt> &&_potential,
                        vert_set_t &&_unstable)
      : vert_map(std::move(_vert_map)), rev_map(std::move(_rev_map)),
        g(std::move(_g)), potential(std::move(_potential)),
        unstable(std::move(_unstable)), _is_bottom(false) {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".copy");

    if (is_top()) {
      // Garbage collection from unconstrained variables in vert_map
      // and rev_map.
      set_to_top();
    }

    CRAB_LOG("tvpi-dbm-size", auto p = size();
             print_dbm_size(p.first, p.second));
  }

  void print_dbm_size(unsigned nodes, unsigned edges) {
    if (nodes > 1) {
      crab::outs() << "#nodes=" << nodes << " "
                   << "#edges=" << edges << " "
                   << "#max-edges=" << nodes * nodes << " "
                   << "(" << ((float)edges / (float)(nodes * nodes)) * 100
                   << ")"
                   << "\n";
    }
  }

#pragma endregion Helpers

  // ============================================================
  // public methods
  // ============================================================
public:
  /// split_dbm_domain implements only standard abstract operations
  /// of a numerical domain so it is intended to be used as a leaf
  /// domain in the hierarchy of domains.
  BOOL_OPERATIONS_NOT_IMPLEMENTED(DBM_t)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(DBM_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(DBM_t)

  // ============================================================
  // constructors
  // ============================================================
#pragma region Constructors
  tvpi_split_dbm_domain(bool is_bottom = false) : _is_bottom(is_bottom) {
    g.growTo(1); // Allocate the zero vector
    potential.push_back(Wt(0));
    rev_map.push_back(boost::none);
  }

  // FIXME: Rewrite to avoid copying if o is _|_
  tvpi_split_dbm_domain(const DBM_t &o)
      : vert_map(o.vert_map), rev_map(o.rev_map), g(o.g),
        potential(o.potential), unstable(o.unstable), _is_bottom(false) {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".copy");
    CRAB_LOG("tvpi-dbm-size", auto p = size();
             print_dbm_size(p.first, p.second));

    if (o._is_bottom)
      set_to_bottom();

    if (!_is_bottom)
      assert(g.size() > 0);
  }

  tvpi_split_dbm_domain(DBM_t &&o)
      : vert_map(std::move(o.vert_map)), rev_map(std::move(o.rev_map)),
        g(std::move(o.g)), potential(std::move(o.potential)),
        unstable(std::move(o.unstable)), _is_bottom(o._is_bottom) {}

  tvpi_split_dbm_domain &operator=(const tvpi_split_dbm_domain &o) {
    TVPI_TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS_ASSIGN_CTOR(".copy");
    if (this != &o) {
      if (o._is_bottom) {
        set_to_bottom();
      } else {
        _is_bottom = false;
        vert_map = o.vert_map;
        rev_map = o.rev_map;
        g = o.g;
        potential = o.potential;
        unstable = o.unstable;
        assert(g.size() > 0);
      }
    }

    CRAB_LOG("tvpi-dbm-size", auto p = size();
             print_dbm_size(p.first, p.second));

    return *this;
  }

  tvpi_split_dbm_domain &operator=(tvpi_split_dbm_domain &&o) {
    if (o._is_bottom) {
      set_to_bottom();
    } else {
      _is_bottom = false;
      vert_map = std::move(o.vert_map);
      rev_map = std::move(o.rev_map);
      g = std::move(o.g);
      potential = std::move(o.potential);
      unstable = std::move(o.unstable);
    }
    return *this;
  }
#pragma endregion Constructors

  // ============================================================
  // _|_ and top
  // ============================================================
#pragma region BotTop
  DBM_t make_top() const override { return DBM_t(false); }

  DBM_t make_bottom() const override { return DBM_t(true); }

  void set_to_top() override {
    tvpi_split_dbm_domain abs(false);
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    vert_map.clear();
    rev_map.clear();
    g.clear();
    potential.clear();
    unstable.clear();
    _is_bottom = true;
  }

  bool is_bottom() const override { return _is_bottom; }

  bool is_top() const override {
    if (_is_bottom)
      return false;
    return g.is_empty();
  }
#pragma endregion BotTop

  // ============================================================
  // main domain operations
  // ============================================================
#pragma region Leq
  bool operator<=(const DBM_t &o) const override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".leq");
    // TODO: check if operator<= cause termination issue

    // cover all trivial cases to avoid allocating a dbm matrix
    if (is_bottom())
      return true;
    else if (o.is_bottom())
      return false;
    else if (o.is_top())
      return true;
    else if (is_top())
      return false;
    else {
      CRAB_LOG("tvpi-dbm-<=", crab::outs() << "Before <=:\n"
                                           << "DBM 1\n"
                                           << *this << "\n"
                                           << "DBM 2\n"
                                           << o << "\n");

      auto leq_op = [](const DBM_t &left, const DBM_t &right) -> bool {
        // left is normalized but right doesn't need to.

        wt_ref_t wx, wy;

        if (left.vert_map.size() < right.vert_map.size()) {
          return false;
        }

        // Set up a mapping from o to this.
        std::vector<unsigned int> vert_renaming(right.g.size(), -1);
        vert_renaming[0] = 0;
        for (auto &p : right.vert_map) {
          if (right.g.succs(p.second).size() == 0 &&
              right.g.preds(p.second).size() == 0)
            continue;

          auto it = left.vert_map.find(key_t(p.first));
          // We can't have this <= o if we're missing some
          // vertex.
          if (it == left.vert_map.end())
            return false;
          vert_renaming[p.second] = (*it).second;
          // vert_renaming[(*it).second] = p.second;
        }

        assert(left.g.size() > 0);
        // GrPerm g_perm(vert_renaming, g);

        for (vert_id ox : right.g.verts()) {
          if (right.g.succs(ox).size() == 0)
            continue;

          assert(vert_renaming[ox] != -1);
          vert_id x = vert_renaming[ox];
          for (auto edge : right.g.e_succs(ox)) {
            vert_id oy = edge.vert;
            assert(vert_renaming[oy] != -1);
            vert_id y = vert_renaming[oy];
            Wt ow = edge.val;
            if (left.g.lookup(x, y, wx) && (wx.get() <= ow))
              continue;
            if (!left.g.lookup(x, 0, wx) || !left.g.lookup(0, y, wy))
              return false;
            if (!(wx.get() + wy.get() <= ow))
              return false;
          }
        }
        return true;
      };

      if (need_normalization()) {
        DBM_t left(*this);
        left.normalize();
        bool res = leq_op(left, o);
        CRAB_LOG("tvpi-dbm-<=", crab::outs()
                                    << "Result <=: " << (res ? "true" : "false")
                                    << "\n");
        return res;
      } else {
        bool res = leq_op(*this, o);
        CRAB_LOG("tvpi-dbm-<=", crab::outs()
                                    << "Result <=: " << (res ? "true" : "false")
                                    << "\n");
        return res;
      }
    }
  }
#pragma endregion Leq

#pragma region Join
  void operator|=(const DBM_t &o) override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".self_join");

    CRAB_LOG("tvpi-dbm", crab::outs() << "Before join:\n"
                                      << "DBM 1\n"
                                      << *this << "\n"
                                      << "DBM 2\n"
                                      << o << "\n");

    if (is_bottom()) {
      *this = o;
    } else if (o.is_top()) {
      set_to_top();
    } else if (is_top() || o.is_bottom()) {
      // do nothing
    } else {

      auto join_op = [](DBM_t &left, const DBM_t &right) {
        // Both left and right are normalized

        check_potential(left.g, left.potential, __LINE__);
        check_potential(right.g, right.potential, __LINE__);

        // Figure out the common renaming, initializing the
        // resulting potentials as we go.
        std::vector<vert_id> perm_x, perm_y;

        std::vector<Wt> pot_rx, pot_ry;
        vert_map_t out_vmap;
        rev_map_t out_revmap;
        // Add the zero vertex
        assert(left.potential.size() > 0);
        pot_rx.push_back(0);
        pot_ry.push_back(0);
        perm_x.push_back(0);
        perm_y.push_back(0);
        out_revmap.push_back(boost::none);

        for (auto &p : left.vert_map) {
          auto it = right.vert_map.find(p.first);
          // Variable exists in both
          if (it != right.vert_map.end()) { // find common vertices
            out_vmap.insert(vmap_elt_t(p.first, perm_x.size()));
            out_revmap.push_back(p.first);
            pot_rx.push_back(left.potential[p.second] - left.potential[0]);
            // XXX JNL: check this out
            // pot_ry.push_back(right.potential[p.second] - right.potential[0]);
            pot_ry.push_back(right.potential[(*it).second] -
                             right.potential[0]);
            perm_x.push_back(p.second);
            perm_y.push_back((*it).second);
          }
        }
        unsigned int sz = perm_x.size();

        // Build the permuted view of x and y.
        assert(left.g.size() > 0);
        GrPerm gx(perm_x, left.g);
        assert(right.g.size() > 0);
        GrPerm gy(perm_y, right.g);

        graph_t join_g = join(gx, gy, sz, pot_rx, pot_ry);
        // Conjecture: join_g remains closed.

        // Now garbage collect any unused vertices
        for (vert_id v : join_g.verts()) {
          if (v == 0)
            continue;
          if (join_g.succs(v).size() == 0 && join_g.preds(v).size() == 0) {
            join_g.forget(v);
            if (out_revmap[v]) {
              out_vmap.erase(*(out_revmap[v]));
              out_revmap[v] = boost::none;
            }
          }
        }

        left.vert_map = std::move(out_vmap);
        left.rev_map = std::move(out_revmap);
        left.g = std::move(join_g);
        left.potential = std::move(pot_rx);
        left.unstable.clear();
        left._is_bottom = false;
        CRAB_LOG("tvpi-dbm", crab::outs() << "Result join:\n" << left << "\n");
      };

      DBM_t &left = *this;
      left.normalize();
      if (o.need_normalization()) {
        DBM_t right(o);
        right.normalize();
        join_op(left, right);
      } else {
        join_op(left, o);
      }
    }
  }

  DBM_t operator|(const DBM_t &o) const override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".join");

    if (is_bottom()) {
      return o;
    } else if (o.is_top() || is_top()) {
      DBM_t res;
      return res;
    } else if (o.is_bottom()) {
      return *this;
    } else {
      CRAB_LOG("tvpi-dbm", crab::outs() << "Before join:\n"
                                        << "DBM 1\n"
                                        << *this << "\n"
                                        << "DBM 2\n"
                                        << o << "\n");

      auto join_op = [](const DBM_t &left, const DBM_t &right) -> DBM_t {
        // Both left and right are normalized

        check_potential(left.g, left.potential, __LINE__);
        check_potential(right.g, right.potential, __LINE__);

        // Figure out the common renaming, initializing the
        // resulting potentials as we go.
        std::vector<vert_id> perm_x, perm_y;
        std::vector<variable_t> perm_inv;

        std::vector<Wt> pot_rx, pot_ry;
        vert_map_t out_vmap;
        rev_map_t out_revmap;
        // Add the zero vertex
        assert(left.potential.size() > 0);
        pot_rx.push_back(0);
        pot_ry.push_back(0);
        perm_x.push_back(0);
        perm_y.push_back(0);
        out_revmap.push_back(boost::none);

        for (auto &p : left.vert_map) {
          auto it = right.vert_map.find(p.first);
          // Variable exists in both
          if (it != right.vert_map.end()) {
            out_vmap.insert(vmap_elt_t(p.first, perm_x.size()));
            out_revmap.push_back(p.first);

            pot_rx.push_back(left.potential[p.second] - left.potential[0]);
            // XXX JNL: check this out
            // pot_ry.push_back(right.potential[p.second] - right.potential[0]);
            pot_ry.push_back(right.potential[(*it).second] -
                             right.potential[0]);
            perm_inv.push_back(p.first.var());
            perm_x.push_back(p.second);
            perm_y.push_back((*it).second);
          }
        }
        unsigned int sz = perm_x.size();

        // Build the permuted view of x and y.
        assert(left.g.size() > 0);
        GrPerm gx(perm_x, left.g);
        assert(right.g.size() > 0);
        GrPerm gy(perm_y, right.g);

        graph_t join_g = join(gx, gy, sz, pot_rx, pot_ry);
        // Conjecture: join_g remains closed.

        // Now garbage collect any unused vertices
        for (vert_id v : join_g.verts()) {
          if (v == 0)
            continue;
          if (join_g.succs(v).size() == 0 && join_g.preds(v).size() == 0) {
            join_g.forget(v);
            if (out_revmap[v]) {
              out_vmap.erase(*(out_revmap[v]));
              out_revmap[v] = boost::none;
            }
          }
        }

        // DBM_t res(join_range, out_vmap, out_revmap, join_g, join_pot);
        DBM_t res(std::move(out_vmap), std::move(out_revmap), std::move(join_g),
                  std::move(pot_rx), vert_set_t());
        // join_g.check_adjs();
        CRAB_LOG("tvpi-dbm", crab::outs() << "Result join:\n" << res << "\n");
        return res;
      };

      if (need_normalization() && o.need_normalization()) {
        DBM_t left(*this);
        DBM_t right(o);
        left.normalize();
        right.normalize();
        return join_op(left, right);
      } else if (need_normalization()) {
        DBM_t left(*this);
        const DBM_t &right = o;
        left.normalize();
        return join_op(left, right);
      } else if (o.need_normalization()) {
        const DBM_t &left = *this;
        DBM_t right(o);
        right.normalize();
        return join_op(left, right);
      } else {
        return join_op(*this, o);
      }
    }
  }
#pragma endregion Join

#pragma region Widening
  DBM_t operator||(const DBM_t &o) const override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".widening");

    if (is_bottom())
      return o;
    else if (o.is_bottom())
      return *this;
    else {
      CRAB_LOG("tvpi-dbm",
               DBM_t left(*this); // to avoid closure on left operand
               crab::outs() << "Before widening:\n"
                            << "DBM 1\n"
                            << left << "\n"
                            << "DBM 2\n"
                            << o << "\n");

      auto widen_op = [](const DBM_t &left, const DBM_t &right) -> DBM_t {
        // Only right is normalized

        // Figure out the common renaming
        std::vector<vert_id> perm_x, perm_y;
        vert_map_t out_vmap;
        rev_map_t out_revmap;
        std::vector<Wt> widen_pot;
        vert_set_t widen_unstable(left.unstable);

        assert(left.potential.size() > 0);
        widen_pot.push_back(Wt(0));
        perm_x.push_back(0);
        perm_y.push_back(0);
        out_revmap.push_back(boost::none);
        for (auto &p : left.vert_map) {
          auto it = right.vert_map.find(p.first);
          // Variable exists in both
          if (it != right.vert_map.end()) {
            out_vmap.insert(vmap_elt_t(p.first, perm_x.size()));
            out_revmap.push_back(p.first);

            widen_pot.push_back(left.potential[p.second] - left.potential[0]);
            perm_x.push_back(p.second);
            perm_y.push_back((*it).second);
          }
        }

        // Build the permuted view of x and y.
        assert(left.g.size() > 0);
        GrPerm gx(perm_x, left.g);
        assert(right.g.size() > 0);
        GrPerm gy(perm_y, right.g);

        // Now perform the widening
        std::vector<vert_id> destabilized;
        graph_t widen_g(split_widen(gx, gy, destabilized, out_revmap));
        for (vert_id v : destabilized)
          widen_unstable.insert(v);

        DBM_t res(std::move(out_vmap), std::move(out_revmap),
                  std::move(widen_g), std::move(widen_pot),
                  std::move(widen_unstable));

        CRAB_LOG("tvpi-dbm", crab::outs() << "Result widening:\n"
                                          << res << "\n");
        return res;
      };

      // Do not normalize left operand
      const DBM_t &left = *this;
      if (o.need_normalization()) {
        DBM_t right(o);
        right.normalize();
        return widen_op(left, right);
      } else {
        return widen_op(left, o);
      }
    }
  }

  DBM_t widening_thresholds(const DBM_t &o,
                            const thresholds<number_t> &ts) const override {
    // TODO: use thresholds
    return (*this || o);
  }
#pragma endregion Widening

#pragma region Meet
  void operator&=(const DBM_t &o) override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".meet");

    if (is_bottom() || o.is_top()) {
      // do nothing
    } else if (is_top() || o.is_bottom()) {
      *this = o;
    } else {
      CRAB_LOG("tvpi-dbm", crab::outs() << "Before meet:\n"
                                        << "DBM 1\n"
                                        << *this << "\n"
                                        << "DBM 2\n"
                                        << o << "\n");

      auto meet_op = [](DBM_t &left, const DBM_t &right) {
        // Both left and right are normalized

        check_potential(left.g, left.potential, __LINE__);
        check_potential(right.g, right.potential, __LINE__);

        // We map vertices in the left operand onto a contiguous range.
        // This will often be the identity map, but there might be gaps.
        vert_map_t meet_verts;
        rev_map_t meet_rev;

        std::vector<vert_id> perm_x, perm_y;
        std::vector<Wt> meet_pi;
        perm_x.push_back(0);
        perm_y.push_back(0);
        meet_pi.push_back(Wt(0));
        meet_rev.push_back(boost::none);
        for (auto &p : left.vert_map) {
          vert_id vv = perm_x.size();
          meet_verts.insert(vmap_elt_t(p.first, vv));
          meet_rev.push_back(p.first);

          perm_x.push_back(p.second);
          perm_y.push_back(-1);
          meet_pi.push_back(left.potential[p.second] - left.potential[0]);
        }

        // Add missing mappings from the right operand.
        for (auto &p : right.vert_map) {
          auto it = meet_verts.find(p.first);

          if (it == meet_verts.end()) {
            vert_id vv = perm_y.size();
            meet_rev.push_back(p.first);

            perm_y.push_back(p.second);
            perm_x.push_back(-1);
            meet_pi.push_back(right.potential[p.second] - right.potential[0]);
            meet_verts.insert(vmap_elt_t(p.first, vv));
          } else {
            perm_y[(*it).second] = p.second;
          }
        }

        // Build the permuted view of x and y.
        assert(left.g.size() > 0);
        GrPerm gx(perm_x, left.g);
        assert(right.g.size() > 0);
        GrPerm gy(perm_y, right.g);

        // Compute the syntactic meet of the permuted graphs.
        bool is_closed;
        graph_t meet_g(GrOps::meet(gx, gy, is_closed));

        // Compute updated potentials on the zero-enriched graph
        // vector<Wt> meet_pi(meet_g.size());
        // We've warm-started pi with the operand potentials
        if (!GrOps::select_potentials(meet_g, meet_pi)) {
          // Potentials cannot be selected -- state is infeasible.
          left.set_to_bottom();
          return;
        }

        if (!is_closed) {
          edge_vector delta;
          SubGraph<graph_t> meet_g_excl(meet_g, 0);
          // GrOps::close_after_meet(meet_g_excl, meet_pi, gx, gy, delta);

          if (crab_domain_params_man::get().zones_chrome_dijkstra())
            GrOps::close_after_meet(meet_g_excl, meet_pi, gx, gy, delta);
          else
            GrOps::close_johnson(meet_g_excl, meet_pi, delta);

          GrOps::apply_delta(meet_g, delta);

          // Recover updated LBs and UBs.
          delta.clear();
          GrOps::close_after_assign(meet_g, meet_pi, 0, delta);
          GrOps::apply_delta(meet_g, delta);
        }

        check_potential(meet_g, meet_pi, __LINE__);

        left.vert_map = std::move(meet_verts);
        left.rev_map = std::move(meet_rev);
        left.g = std::move(meet_g);
        left.potential = std::move(meet_pi);
        left.unstable.clear();
        left._is_bottom = false;

        CRAB_LOG("tvpi-dbm", crab::outs() << "Result meet:\n" << left << "\n");
      };

      DBM_t &left = *this;
      left.normalize();

      if (o.need_normalization()) {
        DBM_t right(o);
        right.normalize();
        meet_op(left, right);
      } else {
        meet_op(left, o);
      }
    }
  }

  DBM_t operator&(const DBM_t &o) const override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".meet");

    if (is_bottom() || o.is_top())
      return *this;
    else if (is_top() || o.is_bottom()) {
      return o;
    } else {
      CRAB_LOG("tvpi-dbm", crab::outs() << "Before meet:\n"
                                        << "DBM 1\n"
                                        << *this << "\n"
                                        << "DBM 2\n"
                                        << o << "\n");

      auto meet_op = [](const DBM_t &left, const DBM_t &right) -> DBM_t {
        // Both left and right are normalized
        check_potential(left.g, left.potential, __LINE__);
        check_potential(right.g, right.potential, __LINE__);

        // We map vertices in the left operand onto a contiguous range.
        // This will often be the identity map, but there might be gaps.
        vert_map_t meet_verts;
        rev_map_t meet_rev;

        std::vector<vert_id> perm_x, perm_y;
        std::vector<Wt> meet_pi;
        perm_x.push_back(0);
        perm_y.push_back(0);
        meet_pi.push_back(Wt(0));
        meet_rev.push_back(boost::none);
        for (auto &p : left.vert_map) {
          vert_id vv = perm_x.size();
          meet_verts.insert(vmap_elt_t(p.first, vv));
          meet_rev.push_back(p.first);

          perm_x.push_back(p.second);
          perm_y.push_back(-1);
          meet_pi.push_back(left.potential[p.second] - left.potential[0]);
        }

        // Add missing mappings from the right operand.
        for (auto &p : right.vert_map) {
          auto it = meet_verts.find(p.first);

          if (it == meet_verts.end()) {
            vert_id vv = perm_y.size();
            meet_rev.push_back(p.first);

            perm_y.push_back(p.second);
            perm_x.push_back(-1);
            meet_pi.push_back(right.potential[p.second] - right.potential[0]);
            meet_verts.insert(vmap_elt_t(p.first, vv));
          } else {
            perm_y[(*it).second] = p.second;
          }
        }

        // Build the permuted view of x and y.
        assert(left.g.size() > 0);
        GrPerm gx(perm_x, left.g);
        assert(right.g.size() > 0);
        GrPerm gy(perm_y, right.g);

        // Compute the syntactic meet of the permuted graphs.
        bool is_closed;
        graph_t meet_g(GrOps::meet(gx, gy, is_closed));

        // Compute updated potentials on the zero-enriched graph
        // vector<Wt> meet_pi(meet_g.size());
        // We've warm-started pi with the operand potentials
        if (!GrOps::select_potentials(meet_g, meet_pi)) {
          // Potentials cannot be selected -- state is infeasible.
          DBM_t res;
          res.set_to_bottom();
          return res;
        }

        if (!is_closed) {
          edge_vector delta;
          SubGraph<graph_t> meet_g_excl(meet_g, 0);
          // GrOps::close_after_meet(meet_g_excl, meet_pi, gx, gy, delta);

          if (crab_domain_params_man::get().zones_chrome_dijkstra())
            GrOps::close_after_meet(meet_g_excl, meet_pi, gx, gy, delta);
          else
            GrOps::close_johnson(meet_g_excl, meet_pi, delta);

          GrOps::apply_delta(meet_g, delta);

          // Recover updated LBs and UBs.
          delta.clear();
          GrOps::close_after_assign(meet_g, meet_pi, 0, delta);
          GrOps::apply_delta(meet_g, delta);
        }
        check_potential(meet_g, meet_pi, __LINE__);
        DBM_t res(std::move(meet_verts), std::move(meet_rev), std::move(meet_g),
                  std::move(meet_pi), vert_set_t());
        CRAB_LOG("tvpi-dbm", crab::outs() << "Result meet:\n" << res << "\n");
        return res;
      };

      if (need_normalization() && o.need_normalization()) {
        DBM_t left(*this);
        DBM_t right(o);
        left.normalize();
        right.normalize();
        return meet_op(left, right);
      } else if (need_normalization()) {
        DBM_t left(*this);
        const DBM_t &right = o;
        left.normalize();
        return meet_op(left, right);
      } else if (o.need_normalization()) {
        const DBM_t &left = *this;
        DBM_t right(o);
        right.normalize();
        return meet_op(left, right);
      } else {
        return meet_op(*this, o);
      }
    }
  }
#pragma endregion Meet

#pragma region Narrowing
  DBM_t operator&&(const DBM_t &o) const override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".narrowing");

    if (is_bottom() || o.is_top())
      return *this;
    else if (is_top() || o.is_bottom()) {
      return o;
    } else {
      // TODO: Implement properly
#if 1
      // Narrowing implemented as meet might not terminate.
      // Make sure that there is always a maximum bound for narrowing
      // iterations.
      return *this & o;
#else
      // Narrowing as a no-op: sound and it will terminate
      CRAB_LOG("tvpi-dbm", crab::outs() << "Before narrowing:\n"
                                        << "DBM 1\n"
                                        << *this << "\n"
                                        << "DBM 2\n"
                                        << o << "\n");

      if (need_normalization()) {
        DBM_t res(*this);
        res.normalize();
        CRAB_LOG("tvpi-dbm", crab::outs() << "Result narrowing:\n"
                                          << res << "\n");
        return res;
      } else {
        CRAB_LOG("tvpi-dbm", crab::outs() << "Result narrowing:\n"
                                          << *this << "\n");
        return *this;
      }
#endif
    }
  }
#pragma endregion Narrowing

#pragma region Assignment
  void assign(const variable_t &x, const linear_expression_t &e) override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".assign");

    if (is_bottom()) {
      return;
    }

    CRAB_LOG("tvpi-dbm-assign", crab::outs() << "--- assigning " << x
                                             << ":=" << e << "\n--[[\n"
                                             << *this << "--]]\n");
    // request to be normalized
    normalize();

    check_potential(g, potential, __LINE__);

    interval_t x_int = eval_interval(e); // get upper and lower bounds for x
    vert_map_t tmp_vert_map;
    key_t x_key(x);

    boost::optional<Wt> lb_w, ub_w; // keep track of the bounds of x
    bool overflow;
    if (x_int.lb().is_finite()) { // -x <= lb
      lb_w = ntow::convert(-(*(x_int.lb().number())), overflow);
      if (overflow) {
        operator-=(x);
        CRAB_LOG("tvpi-dbm-assign", crab::outs()
                                        << "---" << x << ":=" << e << "\n"
                                        << *this << "\n");
        return;
      }
    }
    if (x_int.ub().is_finite()) { // x <= ub
      ub_w = ntow::convert(*(x_int.ub().number()), overflow);
      if (overflow) {
        operator-=(x);
        CRAB_LOG("tvpi-dbm-assign", crab::outs()
                                        << "---" << x << ":=" << e << "\n"
                                        << *this << "\n");
        return;
      }
    }

    bool is_rhs_constant = false;
    // If it's a constant, just assign the interval.
    if (boost::optional<number_t> x_n = x_int.singleton()) {
      set(x, *x_n);
      is_rhs_constant = true;
    }

    if (!is_rhs_constant) {
      std::vector<std::pair<key_t, Wt>> diffs_lb, diffs_ub;
      // Construct difference constraints from the assignment
      // IMPORTANT NOTE: each constraint contains assigned var x
      diffcsts_of_assign(x, e, diffs_lb, diffs_ub);
      if (diffs_lb.size() > 0 || diffs_ub.size() > 0) {
        // Assignment as a sequence of edge additions.
        vert_id v = g.new_vertex(); // graph tracks empty vertex id
        assert(v <= rev_map.size());
        if (v == rev_map.size()) {
          rev_map.push_back(x_key);
          potential.push_back(Wt(0));
        } else {
          potential[v] = Wt(0);
          rev_map[v] = x_key;
        }
        tmp_vert_map.insert(vmap_elt_t(x_key, v));
        Wt_min min_op;
        edge_vector cst_edges;

        // auto get_or_insert_vertex_for_newx = [&](key_t k) -> vert_id {
        //   auto it = tmp_vert_map.find(k);
        //   if (it == tmp_vert_map.end()) {
        //     vert_id v = g.new_vertex(); // graph tracks empty vertex id
        //     assert(v <= rev_map.size());
        //     if (v == rev_map.size()) {
        //       rev_map.push_back(k);
        //       potential.push_back(Wt(0));
        //     } else {
        //       potential[v] = Wt(0);
        //       rev_map[v] = k;
        //     }
        //     tmp_vert_map.insert(vmap_elt_t(k, v));
        //     assert(v != 0);
        //     return v;
        //   } else {
        //     return it->second;
        //   }
        // };

        for (auto diff : diffs_lb) {
          cst_edges.push_back({{v, get_vert(diff.first)}, -diff.second});
          // for (auto &c : crab_domain_params_man::get().coefficients()) {
          //   unsigned newc = diff.first.coeff() * c;
          //   if (auto vby = get_vert_opt(key_t(diff.first.var(), newc))) {
          //     vert_id vax = get_or_insert_vertex_for_newx(key_t(x, c));
          //     cst_edges.push_back({{vax, *vby}, -(diff.second * Wt(c))});
          //   }
          // }
        }

        for (auto diff : diffs_ub) {
          cst_edges.push_back({{get_vert(diff.first), v}, diff.second});
          // for (auto &c : crab_domain_params_man::get().coefficients()) {
          //   unsigned newc = diff.first.coeff() * c;
          //   if (auto vby = get_vert_opt(key_t(diff.first.var(), newc))) {
          //     vert_id vax = get_or_insert_vertex_for_newx(key_t(x, c));
          //     cst_edges.push_back({{*vby, vax}, diff.second * Wt(c)});
          //   }
          // }
        }

        for (auto diff : cst_edges) {
          vert_id src = diff.first.first;
          vert_id dest = diff.first.second;
          g.update_edge(src, diff.second, dest, min_op);
          if (!repair_potential(src, dest)) {
            assert(0 && "Unreachable");
            set_to_bottom();
          }
          check_potential(g, potential, __LINE__);
          close_over_edge(src, dest);
          reduce_tvpi_edge(src, dest, v, tmp_vert_map);
          check_potential(g, potential, __LINE__);
        }

        edge_vector delta;
        GrOps::close_after_assign(g, potential, 0, delta);
        GrOps::apply_delta(g, delta);

        if (lb_w) {
          g.update_edge(v, *lb_w, 0, min_op);
        }
        if (ub_w) {
          g.update_edge(0, *ub_w, v, min_op);
        }

        // Clear the old x vertex
        operator-=(x);
        vert_map.insert(tmp_vert_map.begin(), tmp_vert_map.end());
      } else {
        set(x, x_int);
      }
    }

    check_potential(g, potential, __LINE__);
    CRAB_LOG("tvpi-dbm-assign", crab::outs() << "--- assignment done for " << x
                                             << ":=" << e << "\n--[[\n"
                                             << *this << "--]]\n");
  }
#pragma endregion Assignment

#pragma region Apply1
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".apply");

    if (is_bottom()) {
      return;
    }

    normalize();

    switch (op) {
    case OP_ADDITION:
      assign(x, y + z);
      return;
    case OP_SUBTRACTION:
      assign(x, y - z);
      return;
    case OP_MULTIPLICATION: {
      // evaluate the interval of y and z, if either one is constant,
      // use assign to construct TVPI constraint
      interval_t y_int = get_interval(y);
      interval_t z_int = get_interval(z);
      if (boost::optional<number_t> y_n = y_int.singleton()) {
        assign(x, *y_n * z);
      }
      if (boost::optional<number_t> z_n = z_int.singleton()) {
        assign(x, y * (*z_n));
      }
      // set(x, get_interval(y) * get_interval(z));
      break;
    }
    // For the rest of operations, we fall back on intervals.
    case OP_SDIV:
      set(x, get_interval(y) / get_interval(z));
      break;
    case OP_UDIV:
      set(x, get_interval(y).UDiv(get_interval(z)));
      break;
    case OP_SREM:
      set(x, get_interval(y).SRem(get_interval(z)));
      break;
    default:
      // case OP_UREM:
      set(x, get_interval(y).URem(get_interval(z)));
      break;
    }

    CRAB_LOG("tvpi-dbm-apply", crab::outs()
                                   << "---" << x << ":=" << y << op << z << "\n"
                                   << *this << "\n");
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".apply");

    if (is_bottom()) {
      return;
    }

    normalize();

    switch (op) {
    case OP_ADDITION:
      assign(x, y + k);
      return;
    case OP_SUBTRACTION:
      assign(x, y - k);
      return;
    case OP_MULTIPLICATION:
      assign(x, k * y);
      return;
    // For the rest of operations, we fall back on intervals.
    case OP_SDIV:
      set(x, get_interval(y) / interval_t(k));
      break;
    case OP_UDIV:
      set(x, get_interval(y).UDiv(interval_t(k)));
      break;
    case OP_SREM:
      set(x, get_interval(y).SRem(interval_t(k)));
      break;
    default:
      // case OP_UREM:
      set(x, get_interval(y).URem(interval_t(k)));
      break;
    }

    CRAB_LOG("tvpi-dbm-apply", crab::outs()
                                   << "---" << x << ":=" << y << op << k << "\n"
                                   << *this << "\n");
  }
#pragma endregion Apply1

#pragma region AddEntails
  void operator+=(const linear_constraint_t &cst) {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".add_cst");
    CRAB_LOG("tvpi-dbm-+=", crab::outs() << "---" << cst << "\n"
                                         << *this << "\n");

    auto process_tvpi_expression = [&](const linear_expression_t &expr) {
      for (auto &c : crab_domain_params_man::get().coefficients()) {
        if (auto newexpr = try_rewrite_linear_expression(expr, c)) {
          CRAB_LOG("tvpi-dbm-+=3", crab::outs()
                                       << "processing (" << expr << ") * " << c
                                       << " rewritten " << *newexpr << "\n");
          if (!add_linear_leq(*newexpr)) {
            set_to_bottom();
          }
        }
      }
    };
    if (cst.is_tautology()) {
      return;
    }

    if (cst.is_contradiction()) {
      set_to_bottom();
      return;
    }

    if (is_bottom()) {
      return;
    }

    normalize();
    // g.check_adjs();
    if (cst.is_inequality()) { // e <= c
      if (!add_linear_leq(cst.expression())) {
        set_to_bottom();
      }
      process_tvpi_expression(cst.expression());
    } else if (cst.is_strict_inequality()) { // e < c
      // We try to convert a strict to non-strict.
      auto nc =
          ikos::linear_constraint_impl::strict_to_non_strict_inequality(cst);
      if (nc.is_inequality()) {
        // here we succeed
        if (!add_linear_leq(nc.expression())) {
          set_to_bottom();
        }
        process_tvpi_expression(nc.expression());
      }
    } else if (cst.is_equality()) {
      const linear_expression_t &exp = cst.expression();
      if (!add_linear_leq(exp) || !add_linear_leq(-exp)) {
        set_to_bottom();
      }
      process_tvpi_expression(exp);
      process_tvpi_expression(-exp);
      // g.check_adjs();
    } else if (cst.is_disequation()) {
      // We handle here the case x !=y by converting the disequation
      // into a strict inequality if possible.
      linear_constraint_system_t csts;
      constraint_simp_domain_traits<DBM_t>::lower_disequality(*this, cst, csts);
      for (auto const &c : csts) {
        // We try to convert a strict inequality into non-strict one
        auto nc =
            ikos::linear_constraint_impl::strict_to_non_strict_inequality(c);
        if (nc.is_inequality()) {
          // here we succeed
          if (!add_linear_leq(nc.expression())) {
            set_to_bottom();
          }
        }
      }

      if (!is_bottom()) {
        // We handle here the case x != c
        add_disequation(cst.expression());
      }
    }

    CRAB_LOG("tvpi-dbm-+=", crab::outs() << "---" << cst << "\n"
                                         << *this << "\n");
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    if (is_bottom())
      return;

    for (auto &cst : csts) {
      operator+=(cst);
    }
  }

  virtual bool entails(const linear_constraint_t &cst) const override {
    if (is_bottom()) {
      return true;
    }
    if (cst.is_tautology()) {
      return true;
    }
    if (cst.is_contradiction()) {
      return false;
    }

    bool res;
    if (cst.is_disequation()) {
      // |= c1.x1 + ... + cn.xn != k is iff
      // (1) |= c1.x1 + ... + cn.xn < k OR
      // (2) |= c1.x1 + ... + cn.xn > k
      linear_constraint_t pob1(cst.expression(),
                               linear_constraint_t::kind_t::STRICT_INEQUALITY);
      res = is_unsat(pob1.negate());
      if (!res) {
        linear_constraint_t pob2(
            cst.expression() * number_t(-1),
            linear_constraint_t::kind_t::STRICT_INEQUALITY);
        res = is_unsat(pob2.negate());
      }
    } else if (cst.is_equality()) {
      // |= c1.x1 + ... + cn.xn == k is iff
      // (1) |= c1.x1 + ... + cn.xn <= k AND
      // (2) |= c1.x1 + ... + cn.xn >= k
      linear_constraint_t pob1(cst.expression(),
                               linear_constraint_t::kind_t::INEQUALITY);
      res = is_unsat(pob1.negate());
      if (res) {
        linear_constraint_t pob2(cst.expression() * number_t(-1),
                                 linear_constraint_t::kind_t::INEQUALITY);
        res = is_unsat(pob2.negate());
      }
    } else {
      // cst is an inequality
      res = is_unsat(cst.negate());
    }

    return res;
  }
#pragma endregion AddEntails

#pragma region Intervals
  interval_t operator[](const variable_t &x) override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".to_intervals");
    normalize();
    return (is_bottom() ? interval_t::bottom() : get_interval(vert_map, g, x));
  }

  interval_t at(const variable_t &x) const override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".to_intervals");
    return (is_bottom() ? interval_t::bottom() : get_interval(vert_map, g, x));
  }

  void set(const variable_t &x, interval_t intv) {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".assign");

    if (is_bottom())
      return;

    if (intv.is_bottom()) {
      set_to_bottom();
      return;
    }

    this->operator-=(x);

    if (intv.is_top()) {
      return;
    }

    vert_id v = get_vert(x);
    bool overflow;
    if (intv.ub().is_finite()) {
      Wt ub = ntow::convert(*(intv.ub().number()), overflow);
      if (overflow) {
        return;
      }
      potential[v] = potential[0] + ub;
      g.set_edge(0, ub, v);
    }
    if (intv.lb().is_finite()) {
      Wt lb = ntow::convert(*(intv.lb().number()), overflow);
      if (overflow) {
        return;
      }
      potential[v] = potential[0] + lb;
      g.set_edge(v, -lb, 0);
    }
  }
#pragma endregion Intervals

#pragma region Apply2
  // int_cast_operators_api
  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    int_cast_domain_traits<DBM_t>::apply(*this, op, dst, src);
  }

  // bitwise_operators_api
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".apply");

    if (is_bottom())
      return;
    normalize();

    // Convert to intervals and perform the operation
    interval_t yi = operator[](y);
    interval_t zi = operator[](z);
    interval_t xi = interval_t::bottom();
    switch (op) {
    case OP_AND: {
      xi = yi.And(zi);
      break;
    }
    case OP_OR: {
      xi = yi.Or(zi);
      break;
    }
    case OP_XOR: {
      xi = yi.Xor(zi);
      break;
    }
    case OP_SHL: {
      xi = yi.Shl(zi);
      break;
    }
    case OP_LSHR: {
      xi = yi.LShr(zi);
      break;
    }
    default:
      // case OP_ASHR:
      xi = yi.AShr(zi);
      break;
    }
    set(x, xi);
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".apply");

    if (is_bottom())
      return;
    normalize();

    // Convert to intervals and perform the operation
    interval_t yi = operator[](y);
    interval_t zi(k);
    interval_t xi = interval_t::bottom();

    switch (op) {
    case OP_AND: {
      xi = yi.And(zi);
      break;
    }
    case OP_OR: {
      xi = yi.Or(zi);
      break;
    }
    case OP_XOR: {
      xi = yi.Xor(zi);
      break;
    }
    case OP_SHL: {
      xi = yi.Shl(zi);
      break;
    }
    case OP_LSHR: {
      xi = yi.LShr(zi);
      break;
    }
    default:
      // case OP_ASHR:
      xi = yi.AShr(zi);
      break;
    }
    set(x, xi);
  }
#pragma endregion Apply2

  DEFAULT_SELECT(DBM_t)
  DEFAULT_WEAK_ASSIGN(DBM_t)

#pragma region Projection
  void project(const variable_vector_t &variables) override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".project");

    if (is_bottom() || is_top()) {
      return;
    }

    if (variables.empty()) {
      set_to_top();
      return;
    }

    normalize();

    CRAB_LOG("tvpi-dbm-project", crab::outs() << "Before projecting (";
             tvpi_utils::print_vector(crab::outs(), variables);
             crab::outs() << ")=" << *this << "\n");

    std::vector<bool> save(rev_map.size(), false);
    for (auto &x : variables) {
      for (auto &p : vert_map) {
        if (p.first.var() == x) {
          save[p.second] = true;
        }
      }
    }

    for (vert_id v = 0; v < rev_map.size(); v++) {
      if (!save[v] && rev_map[v]) {
        variable_t vv = rev_map[v]->var();
        operator-=(vv);
      }
    }
    CRAB_LOG("tvpi-dbm-project", crab::outs()
                                     << "After Projection:" << *this << "\n";);
  }

  void operator-=(const variable_t &v) override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".forget");

    if (is_bottom())
      return;
    normalize();

    CRAB_LOG("tvpi-dbm--=", crab::outs()
                                << "Before forget " << v << ": " << g << "\n");
    for (auto it = vert_map.begin(); it != vert_map.end();) {
      if (it->first.var() == v) {
        CRAB_LOG("tvpi-dbm--=", crab::outs()
                                    << "forgetting " << it->second << "\n");
        g.forget(it->second);
        rev_map[it->second] = boost::none;
        it = vert_map.erase(it);
      } else {
        ++it;
      }
    }
    CRAB_LOG("tvpi-dbm--=", crab::outs() << "After: " << g << "\n");
  }

  void forget(const variable_vector_t &variables) override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".forget");

    if (is_bottom() || is_top()) {
      return;
    }

    for (auto &v : variables) {
      operator-=(v);
    }
  }
#pragma endregion Projection

#pragma region Refactoring
  void expand(const variable_t &x, const variable_t &y) override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".expand");

    if (is_bottom() || is_top()) {
      return;
    }

    CRAB_LOG("tvpi-dbm-expand", crab::outs() << "Before expand " << x
                                             << " into " << y << ":\n"
                                             << *this << "\n");

    for (auto &p : vert_map) {
      if (p.first.var() == y) {
        CRAB_ERROR(
            "split_dbm expand operation failed because y already exists");
      }
    }

    std::vector<std::pair<vert_id, vert_id>> to_expand;
    std::vector<unsigned> coeffs;
    for (auto &p : vert_map) {
      if (p.first.var() == x) {
        coeffs.push_back(p.first.coeff());
      }
    }
    for (auto &c : coeffs) {
      vert_id ii = get_vert(key_t(x, c));
      vert_id jj = get_vert(key_t(y, c));
      to_expand.push_back({ii, jj});
    }

    for (auto &p : to_expand) {
      vert_id ii = p.first;
      vert_id jj = p.second;
      edge_vector delta;
      for (auto edge : g.e_preds(ii)) {
        delta.push_back({{edge.vert, jj}, edge.val});
      }

      for (auto edge : g.e_succs(ii)) {
        delta.push_back({{jj, edge.vert}, edge.val});
      }
      GrOps::apply_delta(g, delta);

      potential[jj] = potential[ii];
    }

    CRAB_LOG("tvpi-dbm-expand", crab::outs() << "After expand " << x << " into "
                                             << y << ":\n"
                                             << *this << "\n");
  }

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".rename");

    if (is_top() || is_bottom())
      return;

    CRAB_LOG("tvpi-dbm-rename", crab::outs() << "Renaming {";
             for (auto &v
                  : from) crab::outs()
             << v << ";";
             crab::outs() << "} with "; for (auto &v
                                             : to) crab::outs()
                                        << v << ";";
             crab::outs() << "}:\n"; crab::outs() << *this << "\n";);

    for (unsigned i = 0, sz = from.size(); i < sz; ++i) {
      const variable_t &v = from[i];
      const variable_t &new_v = to[i];
      if (v == new_v) { // nothing to rename
        continue;
      }

      for (auto &p : vert_map) {
        if (p.first.var() == new_v) {
          // We do garbage collection of unconstrained variables only
          // after joins so it's possible to find new_v but we are ok as
          // long as it's unconstrained.
          vert_id dim = p.second;
          if (g.succs(dim).size() != 0 || g.preds(dim).size() != 0) {
            CRAB_ERROR(domain_name() + "::rename assumes that ", new_v,
                       " does not exist");
          }
        }
      }

      vert_map_t tmp_map;

      for (auto it = vert_map.begin(); it != vert_map.end();) {
        if (it->first.var() == v) {
          vert_id dim = it->second;
          key_t nk(new_v, it->first.coeff());
          it = vert_map.erase(it);
          tmp_map.insert(vmap_elt_t(nk, dim));
          rev_map[dim] = nk;
        } else {
          ++it;
        }
      }

      for (auto &kv : tmp_map) {
        vert_map.insert(kv);
      }
    }

    CRAB_LOG("tvpi-dbm-rename", crab::outs() << "RESULT=" << *this << "\n");
  }

  void extract(const variable_t &x, linear_constraint_system_t &csts,
               bool only_equalities) {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".extract");

    normalize();
    if (is_bottom()) {
      return;
    }

    auto it = vert_map.find(key_t(x));
    if (it != vert_map.end()) {
      vert_id s = (*it).second;
      if (rev_map[s]) {
        variable_t vs = rev_map[s]->var();
        SubGraph<graph_t> g_excl(g, 0);
        for (vert_id d : g_excl.verts()) {
          if (rev_map[d]) {
            variable_t vd = rev_map[d]->var();
            // We give priority to equalities since some domains
            // might not understand inequalities
            // FIXME: this is wrong for TVPI constraints
            if (g_excl.elem(s, d) && g_excl.elem(d, s) &&
                g_excl.edge_val(s, d) == Wt(0) &&
                g_excl.edge_val(d, s) == Wt(0)) {
              linear_constraint_t cst(linear_expression_t(vs) == vd);
              csts += cst;
            } else {
              if (!only_equalities && g_excl.elem(s, d)) {
                linear_constraint_t cst(vd - vs <=
                                        number_t(g_excl.edge_val(s, d)));
                csts += cst;
              }
              if (!only_equalities && g_excl.elem(d, s)) {
                linear_constraint_t cst(vs - vd <=
                                        number_t(g_excl.edge_val(d, s)));
                csts += cst;
              }
            }
          }
        }
      }
    }
  }
#pragma endregion Refactoring

  void normalize() override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".normalize");

    // Always maintained in normal form, except for widening
    if (!need_normalization()) {
      return;
    }

    SubGraph<graph_t> g_excl(g, 0);
    edge_vector delta;
    // GrOps::close_after_widen(g, potential, vert_set_wrap_t(unstable), delta);
    // GKG: Check
    if (crab_domain_params_man::get().zones_widen_restabilize())
      GrOps::close_after_widen(g_excl, potential, vert_set_wrap_t(unstable),
                               delta);
    else
      GrOps::close_johnson(g_excl, potential, delta);
    // Retrive variable bounds
    GrOps::close_after_assign(g, potential, 0, delta);

    GrOps::apply_delta(g, delta);

    unstable.clear();
  }

#pragma region Printing
  // Output function
  void write(crab_os &o) const override {
    // linear_constraint_system_t inv = to_linear_constraint_system();
    // o << inv;

    if (need_normalization()) {
      DBM_t tmp(*this);
      tmp.normalize();
      write(o, tmp);
    } else {
      write(o, *this);
    }
  }

  void print_potentials(crab_os &o) const {
    o << "\tPotentials={";
    for (vert_id v = 0; v < rev_map.size(); v++) {
      if (v == 0) {
        o << "v0=" << potential[v] << ";";
      }
      if (rev_map[v]) {
        o << (*rev_map[v]) << "=" << potential[v] << ";";
      }
    }
    o << "}\n";
  }

  void print_graph(crab_os &o) const {
    o << "\tGraphs=";
    g.write(o);
    o << "\n";
  }

  void print_details(crab_os &o) const {
    o << "DBM:\n";
    write(o);
    o << "\n";
    o << "Details:\n";
    print_potentials(o);
    print_graph(o);
  }
#pragma endregion Printing

  linear_constraint_system_t to_linear_constraint_system() const override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".to_linear_constraints");

    if (need_normalization()) {
      DBM_t tmp(*this);
      tmp.normalize();
      return to_linear_constraint_system(tmp);
    } else {
      return to_linear_constraint_system(*this);
    }
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    auto lin_csts = to_linear_constraint_system();
    if (lin_csts.is_false()) {
      return disjunctive_linear_constraint_system_t(true /*is_false*/);
    } else if (lin_csts.is_true()) {
      return disjunctive_linear_constraint_system_t(false /*is_false*/);
    } else {
      return disjunctive_linear_constraint_system_t(lin_csts);
    }
  }

  // return number of vertices and edges
  std::pair<std::size_t, std::size_t> size() const {
    return {g.size(), g.num_edges()};
  }

  std::vector<variable_t> vars() const {
    std::vector<variable_t> res;
    std::unordered_set<variable_t> var_set;
    var_set.reserve(vert_map.size());
    for (const auto &p : vert_map) {
      var_set.insert(p.first);
    }
    return std::vector<variable_t>(res.begin(), res.end());
  }

  bool exists(const variable_t &x) const {
    return vert_map.find(key_t(x)) != vert_map.end();
  }

#pragma region NotUsed
  // intrinsics operations
  void intrinsic(std::string name, const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    // CRAB_WARN("intrinsic ", name, " not implemented by ", domain_name());
  }

  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const DBM_t &invariant) override {
    // CRAB_WARN("backward_intrinsic", name, " not implemented by ",
    // domain_name());
  }

  void minimize() override {
    // CRAB_WARN("minimize ", name, " not implemented by ", domain_name());
  }

  void callee_entry(const callsite_info<variable_t> &callsite,
                    const DBM_t &caller) override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".callee_entry");
    inter_abstract_operations<
        DBM_t,
        DomainParams::implement_inter_transformers>::callee_entry(callsite,
                                                                  caller,
                                                                  *this);
  }

  void caller_continuation(const callsite_info<variable_t> &callsite,
                           const DBM_t &callee) override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".caller_cont");
    inter_abstract_operations<DBM_t,
                              DomainParams::implement_inter_transformers>::
        caller_continuation(callsite, callee, *this);
  }

  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const DBM_t &inv) override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".backward_assign");

    crab::domains::BackwardAssignOps<DBM_t>::assign(*this, x, e, inv);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const DBM_t &inv) override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".backward_apply");

    crab::domains::BackwardAssignOps<DBM_t>::apply(*this, op, x, y, z, inv);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const DBM_t &inv) override {
    TVPI_SPLIT_DBM_DOMAIN_SCOPED_STATS(".backward_apply");

    crab::domains::BackwardAssignOps<DBM_t>::apply(*this, op, x, y, z, inv);
  }

  std::string domain_name() const override { return "TVPISplitDBM"; }
#pragma endregion NotUsed
}; // class tvpi_split_dbm_domain

template <typename Number, typename VariableName, typename DBMParams,
          typename DomainParams>
struct abstract_domain_traits<
    tvpi_split_dbm_domain<Number, VariableName, DBMParams, DomainParams>> {
  using number_t = Number;
  using varname_t = VariableName;
};

template <typename Number, typename VariableName, typename DBMParams,
          typename DomainParams>
class reduced_domain_traits<
    tvpi_split_dbm_domain<Number, VariableName, DBMParams, DomainParams>> {
public:
  using tvpi_sdbm_domain_t =
      tvpi_split_dbm_domain<Number, VariableName, DBMParams, DomainParams>;
  using variable_t = typename tvpi_sdbm_domain_t::variable_t;
  using linear_constraint_system_t =
      typename tvpi_sdbm_domain_t::linear_constraint_system_t;

  static void extract(tvpi_sdbm_domain_t &dom, const variable_t &x,
                      linear_constraint_system_t &csts, bool only_equalities) {
    dom.extract(x, csts, only_equalities);
  }
};
} // namespace domains
} // namespace crab

#pragma GCC diagnostic pop
