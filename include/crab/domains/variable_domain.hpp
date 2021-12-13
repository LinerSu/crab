#pragma once

/* Classical variable propagation domain */

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/backward_assign_operations.hpp>
#include <crab/domains/interval.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/support/stats.hpp>

#include <boost/optional.hpp>

namespace crab {
namespace domains {
/**
 *  Each address variable is mapped to an element in this lattice:
 *
 *                 top
 *                  |
 * ...,Var_1,Var_2,...,Var_{n-1},Var_n,...
 *                  |
 *                bottom
 *
 * top means that it might refer to multiple objects
 **/
template <typename VariableName> class variable_domain {
  boost::optional<VariableName> m_variable;
  bool m_is_bottom;
  using variable_t = variable_domain<VariableName>;
  variable_domain(bool is_bottom)
      : m_variable(boost::none), m_is_bottom(is_bottom) {}

public:
  variable_domain(VariableName c) : m_variable(c), m_is_bottom(false) {}

  static variable_t bottom() { return variable_domain(true); }

  static variable_t top() { return variable_domain(false); }

  // static variable_t singleton() { return variable_domain(VariableName(0)); }

  bool is_bottom() const { return m_is_bottom; }

  bool is_top() const { return (!is_bottom() && !m_variable); }

  bool is_variable() const { return m_variable != boost::none; }

  VariableName get_variable() const {
    assert(is_variable());
    return *m_variable;
  }

  /** Begin Lattice Operations **/

  /**
    Indicates either this <= v property based on lattice.

    @param v a variable object has the same type as current.
    @return true if <= holds; else false.
  */
  bool operator<=(const variable_t &v) const {
    if (is_bottom() || v.is_top()) {
      return true;
    } else if (v.is_bottom() || is_top()) {
      return false;
    } else {
      assert(is_variable());
      assert(v.is_variable());
      return get_variable() == v.get_variable();
    }
  }

  /**
    Indicates either this == v property based on lattice.
    Note top == top?

    @param v a variable object has the same type as current.
    @return true if == holds; else false.
  */
  bool operator==(const variable_t &v) const {
    return (m_is_bottom == v.m_is_bottom && m_variable == v.m_variable);
  }

  /**
    Indicates either this | v property based on lattice.
    bot | v = v; this | bot = this
    if this == v then this | v = this
    otherwise return top

    @param v a variable object has the same type as current.
    @return what we described above.
  */
  variable_t operator|(const variable_t &v) const {
    if (is_bottom() || v.is_top())
      return v;
    else if (is_top() || v.is_bottom())
      return *this;
    else {
      assert(is_variable());
      assert(v.is_variable());
      if (get_variable() == v.get_variable()) {
        return *this;
      } else {
        return variable_t::top();
      }
    }
  }

  variable_t operator||(const variable_t &v) const { return *this | v; }

  template <typename Thresholds>
  variable_t widening_thresholds(const variable_t &v,
                                 const Thresholds &ts /*unused*/) const {
    return *this | v;
    // Q: how we handle widening with threshold
  }

  /**
    Indicates either this & v property based on lattice.
    this & top = this; top & v = v
    if this == v then this & v = this
    otherwise return bot

    @param v a variable object has the same type as current.
    @return what we described above.
  */
  variable_t operator&(const variable_t &v) const {
    if (is_bottom() || v.is_top())
      return *this;
    else if (is_top() || v.is_bottom()) {
      return v;
    } else {
      assert(is_variable());
      assert(v.is_variable());
      if (get_variable() == v.get_variable()) {
        return *this;
      } else {
        return variable_t::bottom();
      }
    }
  }

  variable_t operator&&(const variable_t &v) const { return *this & v; }
  /** End  Lattice Operations **/

  /** Begin arithmetic operations **/
  // The following operations -- op
  // Var_i op Var_j = top
  variable_t Add(const variable_t &v) const {
    if (is_bottom() || v.is_bottom()) {
      return variable_t::bottom();
    } else {
      return variable_t::top();
    }
  }
  variable_t Sub(const variable_t &v) const {
    if (is_bottom() || v.is_bottom()) {
      return variable_t::bottom();
    } else {
      return variable_t::top();
    }
  }
  variable_t Mul(const variable_t &v) const {
    if (is_bottom() || v.is_bottom()) {
      return variable_t::bottom();
    } else {
      return variable_t::top();
    }
  }
  variable_t SDiv(const variable_t &v) const {
    if (is_bottom() || v.is_bottom()) {
      return variable_t::bottom();
    } else {
      return variable_t::top();
    }
  }
  variable_t SRem(const variable_t &v) const {
    if (is_bottom() || v.is_bottom()) {
      return variable_t::bottom();
    } else {
      return variable_t::top();
    }
  }
  variable_t UDiv(const variable_t &v) const {
    if (is_bottom() || v.is_bottom()) {
      return variable_t::bottom();
    } else {
      return variable_t::top();
    }
  }
  variable_t URem(const variable_t &v) const {
    if (is_bottom() || v.is_bottom()) {
      return variable_t::bottom();
    } else {
      return variable_t::top();
    }
  }
  variable_t Arith(crab::domains::arith_operation_t op,
                   const variable_t &v) const {
    if (is_bottom() || v.is_bottom()) {
      return variable_t::bottom();
    } else {
      return variable_t::top();
    }
  }
  /** End arithmetic operations **/

  /** Begin bitwise operations **/
  // These operations depend on the type of Number
  variable_t BitwiseAnd(const variable_t &v) const {
    if (is_bottom() || v.is_bottom()) {
      return variable_t::bottom();
    } else {
      return variable_t::top();
    }
  }
  variable_t BitwiseOr(const variable_t &v) const {
    if (is_bottom() || v.is_bottom()) {
      return variable_t::bottom();
    } else {
      return variable_t::top();
    }
  }
  variable_t BitwiseXor(const variable_t &v) const {
    if (is_bottom() || v.is_bottom()) {
      return variable_t::bottom();
    } else {
      return variable_t::top();
    }
  }
  variable_t BitwiseShl(const variable_t &v) const {
    if (is_bottom() || v.is_bottom()) {
      return variable_t::bottom();
    } else {
      return variable_t::top();
    }
  }
  variable_t BitwiseLShr(const variable_t &v) const {
    if (is_bottom() || v.is_bottom()) {
      return variable_t::bottom();
    } else {
      return variable_t::top();
    }
  }
  variable_t BitwiseAShr(const variable_t &v) const {
    if (is_bottom() || v.is_bottom()) {
      return variable_t::bottom();
    } else {
      return variable_t::top();
    }
  }
  variable_t BitwiseOp(const variable_t &v) const {
    if (is_bottom() || v.is_bottom()) {
      return variable_t::bottom();
    } else {
      return variable_t::top();
    }
  }
  /** End bitwise operations **/

  void write(crab::crab_os &o) const {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "top";
    } else {
      assert(is_variable());
      o << get_variable();
    }
  }

  // void rename(const variable_vector_t &from,
  //             const variable_vector_t &to) {
  //   if (!is_variable()) {
  //     return;
  //   }
  //   for (int i = 0; i < from.size(); ++i) {
  //     auto tmp = from[i].name();
  //     if (tmp == m_variable) {
  //       m_variable = to[i].name();
  //       break;
  //     }
  //   }
  // }

  // void forget(const variable_vector_t &variables) {
  //   if (is_bottom() || is_top()) {
  //     return;
  //   }
  //   for (auto const &var : variables) {
  //     auto tmp = from[i].name();
  //     if (m_variable == tmp) {
  //       m_variable = boost::none;
  //       m_is_bottom = true;
  //       break;
  //     }
  //   }
  // }

  friend inline crab_os &operator<<(crab_os &o, const variable_t &v) {
    v.write(o);
    return o;
  }
};

// template <typename VariableName>
// class variable_domain final
//     : public
//     crab::domains::abstract_domain_api<variable_domain<VariableName>> {
// public:
//   using variable_domain_t = variable_domain<VariableName>;
//   using abstract_domain_t =
//       crab::domains::abstract_domain_api<variable_domain_t>;
//   using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
//   using typename abstract_domain_t::interval_t;
//   using typename abstract_domain_t::linear_constraint_system_t;
//   using typename abstract_domain_t::linear_constraint_t;
//   using typename abstract_domain_t::linear_expression_t;
//   using typename abstract_domain_t::reference_constraint_t;
//   using typename abstract_domain_t::variable_or_constant_t;
//   using typename abstract_domain_t::variable_or_constant_vector_t;
//   using typename abstract_domain_t::variable_t;
//   using typename abstract_domain_t::variable_vector_t;
//   using variable_value_t = variable<VariableName>;
//   using varname_t = VariableName;

// private:
//   using separate_domain_t = ikos::separate_domain<variable_t,
//   variable_value_t>;

// private:
//   separate_domain_t m_env;

//   variable_domain(separate_domain_t &&env) : m_env(std::move(env)) {}

//   variable_value_t eval(const linear_expression_t &expr) const {
//     assert(!is_bottom());
//     variable_value_t r(expr.constant());
//     for (auto const &kv : expr) {
//       variable_value_t c(kv.first);
//       r = r.Add(c.Mul(m_env[kv.second]));
//       if (r.is_top()) {
//         break;
//       }
//     }
//     return r;
//   }

//   constant_t compute_residual(const linear_constraint_t &cst,
//                               const variable_value_t &pivot) const {
//     constant_t residual(cst.constant());
//     for (auto const &kv : cst) {
//       constant_t c(kv.first);
//       const variable_value_t &v = kv.second;
//       if (!(v == pivot)) {
//         residual = residual.Sub(c.Mul(m_env[v]));
//         if (residual.is_top()) {
//           break;
//         }
//       }
//     }
//     return residual;
//   }

//   void propagate(const linear_constraint_t &cst) {
//     if (is_bottom()) {
//       return;
//     }

//     if (cst.is_inequality() || cst.is_strict_inequality() ||
//         cst.is_disequation()) {
//       constant_t e = eval(cst.expression());
//       if (e.is_constant()) {
//         if (cst.is_inequality()) {
//           if (!(e.get_constant() <= number_t(0))) {
//             set_to_bottom();
//             return;
//           }
//         } else if (cst.is_disequation()) {
//           if (!(e.get_constant() != number_t(0))) {
//             set_to_bottom();
//             return;
//           }
//         } else if (cst.is_strict_inequality()) {
//           if (!(e.get_constant() < number_t(0))) {
//             set_to_bottom();
//             return;
//           }
//         }
//       }
//     } else if (cst.is_equality()) {
//       for (auto kv : cst) {
//         number_t c = kv.first;
//         const variable_t &pivot = kv.second;
//         constant_t new_c = compute_residual(cst, pivot).SDiv(c);
//         if (!new_c.is_top()) {
//           m_env.set(pivot, m_env[pivot] & new_c);
//         }
//       }
//     }
//   }

//   void solve_constraints(const linear_constraint_system_t &csts) {
//     for (auto const &c : csts) {
//       if (is_bottom()) {
//         return;
//       }
//       if (c.is_inequality() && c.is_unsigned()) {
//         // we don't handle unsigned constraints
//         continue;
//       }
//       if (c.is_tautology()) {
//         continue;
//       }
//       if (c.is_contradiction()) {
//         set_to_bottom();
//         return;
//       }
//       propagate(c);
//     }
//   }

// public:
//   variable_domain_t make_top() const override {
//     return variable_domain_t(separate_domain_t::top());
//   }

//   variable_domain_t make_bottom() const override {
//     return variable_domain_t(separate_domain_t::bottom());
//   }

//   void set_to_top() override {
//     variable_domain abs(separate_domain_t::top());
//     std::swap(*this, abs);
//   }

//   void set_to_bottom() override {
//     variable_domain abs(separate_domain_t::bottom());
//     std::swap(*this, abs);
//   }

//   variable_domain() : m_env(separate_domain_t::top()) {}

//   variable_domain(const variable_domain_t &e) : m_env(e.m_env) {
//     crab::CrabStats::count(domain_name() + ".count.copy");
//     crab::ScopedCrabStats __st__(domain_name() + ".copy");
//   }

//   variable_domain(variable_domain_t &&e) : m_env(std::move(e.m_env)) {}

//   variable_domain_t &operator=(const variable_domain_t &v) {
//     crab::CrabStats::count(domain_name() + ".count.copy");
//     crab::ScopedCrabStats __st__(domain_name() + ".copy");
//     if (this != &v) {
//       m_env = v.m_env;
//     }
//     return *this;
//   }

//   variable_domain_t &operator=(variable_domain_t &&v) {
//     if (this != &v) {
//       m_env = std::move(v.m_env);
//     }
//     return *this;
//   }

//   variable_value_t get_variable_value(const variable_t &v) const {
//     return m_env[v];
//   }

//   void set_variable_value(const variable_t &v, variable_value_t v_val) {
//     m_env.set(v, v_val);
//   }

//   bool is_bottom() const override { return m_env.is_bottom(); }

//   bool is_top() const override { return m_env.is_top(); }

//   bool operator<=(const variable_domain_t &v) const override {
//     crab::CrabStats::count(domain_name() + ".count.leq");
//     crab::ScopedCrabStats __st__(domain_name() + ".leq");
//     return (m_env <= v.m_env);
//   }

//   void operator|=(const variable_domain_t &v) override {
//     crab::CrabStats::count(domain_name() + ".count.join");
//     crab::ScopedCrabStats __st__(domain_name() + ".join");
//     CRAB_LOG("constant-domain",
//              crab::outs() << "Join " << m_env << " and " << v.m_env <<
//              "\n";);
//     m_env = m_env | v.m_env;
//     CRAB_LOG("constant-domain", crab::outs() << "Res=" << m_env << "\n";);
//   }

//   variable_domain_t operator|(const variable_domain_t &v) const override {
//     crab::CrabStats::count(domain_name() + ".count.join");
//     crab::ScopedCrabStats __st__(domain_name() + ".join");
//     return (m_env | v.m_env);
//   }

//   variable_domain_t operator&(const variable_domain_t &v) const override {
//     crab::CrabStats::count(domain_name() + ".count.meet");
//     crab::ScopedCrabStats __st__(domain_name() + ".meet");
//     return (m_env & v.m_env);
//   }

//   variable_domain_t operator||(const variable_domain_t &v) const override {
//     crab::CrabStats::count(domain_name() + ".count.widening");
//     crab::ScopedCrabStats __st__(domain_name() + ".widening");
//     return (m_env || v.m_env);
//   }

//   variable_domain_t widening_thresholds(
//       const variable_domain_t &v,
//       const crab::iterators::thresholds<number_t> &ts) const override {
//     crab::CrabStats::count(domain_name() + ".count.widening");
//     crab::ScopedCrabStats __st__(domain_name() + ".widening");
//     return m_env.widening_thresholds(v.m_env, ts);
//   }

//   variable_domain_t operator&&(const variable_domain_t &v) const override {
//     crab::CrabStats::count(domain_name() + ".count.narrowing");
//     crab::ScopedCrabStats __st__(domain_name() + ".narrowing");
//     return (m_env && v.m_env);
//   }

//   void operator-=(const variable_t &v) override {
//     crab::CrabStats::count(domain_name() + ".count.forget");
//     crab::ScopedCrabStats __st__(domain_name() + ".forget");
//     m_env -= v;
//   }

//   interval_t operator[](const variable_t &v) override {
//     variable_value_t val = m_env[v];
//     if (val.is_bottom()) {
//       return interval_t::bottom();
//     } else if (val.is_top()) {
//       return interval_t::top();
//     } else {
//       assert(val.is_variable());
//       return interval_t(val.get_variable());
//     }
//   }

//   // // TODO: update this
//   // void operator+=(const linear_constraint_system_t &csts) override {
//   //   crab::CrabStats::count(domain_name() + ".count.add_constraints");
//   //   crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");
//   //   solve_constraints(csts);
//   // }

//   // // TODO: update this
//   // void assign(const variable_t &x, const linear_expression_t &e) override
//   {
//   //   crab::CrabStats::count(domain_name() + ".count.assign");
//   //   crab::ScopedCrabStats __st__(domain_name() + ".assign");
//   //   if (boost::optional<variable_t> v = e.get_variable()) {
//   //     m_env.set(x, m_env[(*v)]);
//   //   } else {
//   //     m_env.set(x, eval(e));
//   //   }
//   // }

//   // void apply(crab::domains::arith_operation_t op, const variable_t &x,
//   //            const variable_t &y, const variable_t &z) override {
//   //   crab::CrabStats::count(domain_name() + ".count.apply");
//   //   crab::ScopedCrabStats __st__(domain_name() + ".apply");

//   //   if (!is_bottom()) {
//   //     variable_value_t yc = m_env[y];
//   //     variable_value_t zc = m_env[z];
//   //     variable_value_t xc = variable_value_t::top();
//   //     xc = yc.Arith(op, zc);
//   //     m_env.set(x, xc);
//   //   }
//   // }

//   // void apply(crab::domains::arith_operation_t op, const variable_t &x,
//   //            const variable_t &y, number_t k) override {
//   //   crab::CrabStats::count(domain_name() + ".count.apply");
//   //   crab::ScopedCrabStats __st__(domain_name() + ".apply");

//   //   // Var_i op num(k) should not be applied
//   //   CRAB_ERROR("Operation ", op, " not supported for variable domain");

//   //   if (!is_bottom()) {
//   //     variable_value_t xc = variable_value_t::top();
//   //     m_env.set(x, xc);
//   //   }
//   // }

//   // intrinsics operations
//   void intrinsic(std::string name, const variable_or_constant_vector_t
//   &inputs,
//                  const variable_vector_t &outputs) override {
//     CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
//   }

//   // void backward_intrinsic(std::string name,
//   //                         const variable_or_constant_vector_t &inputs,
//   //                         const variable_vector_t &outputs,
//   //                         const constant_domain_t &invariant) override {
//   //   CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
//   // }

//   // // backward arithmetic operations
//   // void backward_assign(const variable_t &x, const linear_expression_t &e,
//   //                      const constant_domain_t &inv) override {
//   //   crab::CrabStats::count(domain_name() + ".count.backward_assign");
//   //   crab::ScopedCrabStats __st__(domain_name() + ".backward_assign");
//   //   // TODO
//   // }

//   // void backward_apply(crab::domains::arith_operation_t op, const
//   variable_t &x,
//   //                     const variable_t &y, number_t z,
//   //                     const constant_domain_t &inv) override {
//   //   crab::CrabStats::count(domain_name() + ".count.backward_apply");
//   //   crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");
//   //   // TODO
//   // }

//   // void backward_apply(crab::domains::arith_operation_t op, const
//   variable_t &x,
//   //                     const variable_t &y, const variable_t &z,
//   //                     const constant_domain_t &inv) override {
//   //   crab::CrabStats::count(domain_name() + ".count.backward_apply");
//   //   crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");
//   //   // TODO
//   // }

//   // cast operations
//   // void apply(crab::domains::int_conv_operation_t /*op*/, const variable_t
//   &dst,
//   //            const variable_t &src) override {
//   //   // ignore the widths
//   //   assign(dst, src);
//   // }

//   // bitwise operations
//   // void apply(crab::domains::bitwise_operation_t op, const variable_t &x,
//   //            const variable_t &y, const variable_t &z) override {
//   //   crab::CrabStats::count(domain_name() + ".count.apply");
//   //   crab::ScopedCrabStats __st__(domain_name() + ".apply");

//   //   if (!is_bottom()) {
//   //     variable_value_t yc = m_env[y];
//   //     variable_value_t zc = m_env[z];
//   //     variable_value_t xc = variable_value_t::top();
//   //     xc = yc.BitwiseOp(zc);
//   //     m_env.set(x, xc);
//   //   }
//   // }

//   // void apply(crab::domains::bitwise_operation_t op, const variable_t &x,
//   //            const variable_t &y, number_t k) override {
//   //   crab::CrabStats::count(domain_name() + ".count.apply");
//   //   crab::ScopedCrabStats __st__(domain_name() + ".apply");

//   //   // Var_i op num(k) should not be applied
//   //   CRAB_ERROR("Operation ", op, " not supported for variable domain");

//   //   if (!is_bottom()) {
//   //     variable_value_t xc = variable_value_t::top();
//   //     m_env.set(x, xc);
//   //   }
//   // }

//   // virtual void select(const variable_t &lhs, const linear_constraint_t
//   &cond,
//   //                     const linear_expression_t &e1,
//   //                     const linear_expression_t &e2) override {
//   //   crab::CrabStats::count(domain_name() + ".count.select");
//   //   crab::ScopedCrabStats __st__(domain_name() + ".select");

//   //   if (!is_bottom()) {
//   //     constant_domain_t inv1(*this);
//   //     inv1 += cond;
//   //     if (inv1.is_bottom()) {
//   //       assign(lhs, e2);
//   //       return;
//   //     }
//   //     constant_domain_t inv2(*this);
//   //     inv2 += cond.negate();
//   //     if (inv2.is_bottom()) {
//   //       assign(lhs, e1);
//   //       return;
//   //     }
//   //     m_env.set(lhs, eval(e1) | eval(e2));
//   //   }
//   // }

//   /// variable_domain_t implements only standard abstract operations of
//   /// a equality domain so it is intended to be used as a leaf domain
//   /// in the hierarchy of domains.
//   NUMERICAL_OPERATIONS_NOT_IMPLEMENTED(variable_domain_t)
//   BOOL_OPERATIONS_NOT_IMPLEMENTED(variable_domain_t)
//   ARRAY_OPERATIONS_NOT_IMPLEMENTED(variable_domain_t)
//   REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(variable_domain_t)

//   void forget(const variable_vector_t &variables) override {
//     if (is_bottom() || is_top()) {
//       return;
//     }
//     for (auto const &var : variables) {
//       this->operator-=(var);
//     }
//   }

//   void project(const variable_vector_t &variables) override {
//     crab::CrabStats::count(domain_name() + ".count.project");
//     crab::ScopedCrabStats __st__(domain_name() + ".project");

//     m_env.project(variables);
//   }

//   void rename(const variable_vector_t &from,
//               const variable_vector_t &to) override {
//     crab::CrabStats::count(domain_name() + ".count.rename");
//     crab::ScopedCrabStats __st__(domain_name() + ".rename");

//     m_env.rename(from, to);
//   }

//   void expand(const variable_t &x, const variable_t &new_x) override {
//     crab::CrabStats::count(domain_name() + ".count.expand");
//     crab::ScopedCrabStats __st__(domain_name() + ".expand");

//     if (is_bottom() || is_top()) {
//       return;
//     }

//     m_env.set(new_x, m_env[x]);
//   }

//   void normalize() override {}

//   void minimize() override {}

//   void write(crab::crab_os &o) const override {
//     crab::CrabStats::count(domain_name() + ".count.write");
//     crab::ScopedCrabStats __st__(domain_name() + ".write");

//     m_env.write(o);
//   }

//   linear_constraint_system_t to_linear_constraint_system() const override {
//     crab::CrabStats::count(domain_name() +
//                            ".count.to_linear_constraint_system");
//     crab::ScopedCrabStats __st__(domain_name() +
//                                  ".to_linear_constraint_system");

//     linear_constraint_system_t csts;

//     if (this->is_bottom()) {
//       csts += linear_constraint_t::get_false();
//       return csts;
//     }

//     for (auto it = m_env.begin(); it != m_env.end(); ++it) {
//       const variable_t &v = it->first;
//       const variable_value_t &val = it->second;
//       if (val.is_variable()) {
//         csts += linear_constraint_t(v == val.get_variable());
//       }
//     }
//     return csts;
//   }

//   disjunctive_linear_constraint_system_t
//   to_disjunctive_linear_constraint_system() const override {
//     auto lin_csts = to_linear_constraint_system();
//     if (lin_csts.is_false()) {
//       return disjunctive_linear_constraint_system_t(true /*is_false*/);
//     } else if (lin_csts.is_true()) {
//       return disjunctive_linear_constraint_system_t(false /*is_false*/);
//     } else {
//       return disjunctive_linear_constraint_system_t(lin_csts);
//     }
//   }

//   std::string domain_name() const override { return "VariableDomain"; }

// }; // class variable_domain
} // namespace domains
} // namespace crab

// namespace crab {
// namespace domains {
// template <typename VariableName>
// struct abstract_domain_traits<variable_domain<VariableName>> {
//   using number_t = Number;
//   using varname_t = VariableName;
// };

// } // namespace domains
// } // namespace crab
