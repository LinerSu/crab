#pragma once

/*
   Implementation of the abstract transfer functions by reducing them
   to abstract domain operations.

   These are the main Crab statements for which we define their abstract
   transfer functions:

   ARITHMETIC and BOOLEAN
     x := y bin_op z;
     x := y;
     assume(cst)
     assert(cst);
     x := select(cond, y, z);

   ARRAYS
     a[l...u] := v (a,b are arrays and v can be bool or integer)
     a[i] := v;
     v := a[i];
     a := b

   REFERENCES
     TODO

   FUNCTIONS
     x := foo(arg1,...,argn);
     return r;

   havoc(x);

 */

#include <crab/cfg/cfg.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>
#include <crab/domains/abstract_domain_operators.hpp>

namespace crab {  
namespace analyzer {

/**
 * API abstract transformer
 **/
template <typename Number, typename VariableName>
class abs_transformer_api
    : public crab::cfg::statement_visitor<Number, VariableName> {
public:
  typedef Number number_t;
  typedef VariableName varname_t;

  typedef variable<number_t, VariableName> var_t;
  typedef ikos::linear_expression<number_t, VariableName> lin_exp_t;
  typedef ikos::linear_constraint<number_t, VariableName> lin_cst_t;
  typedef ikos::linear_constraint_system<number_t, VariableName> lin_cst_sys_t;

  typedef crab::cfg::havoc_stmt<number_t, VariableName> havoc_t;
  typedef crab::cfg::unreachable_stmt<number_t, VariableName> unreach_t;

  typedef crab::cfg::binary_op<number_t, VariableName> bin_op_t;
  typedef crab::cfg::assignment<number_t, VariableName> assign_t;
  typedef crab::cfg::assume_stmt<number_t, VariableName> assume_t;
  typedef crab::cfg::select_stmt<number_t, VariableName> select_t;
  typedef crab::cfg::assert_stmt<number_t, VariableName> assert_t;
  typedef crab::cfg::int_cast_stmt<number_t, VariableName> int_cast_t;
  
  typedef crab::cfg::callsite_stmt<number_t, VariableName> callsite_t;
  typedef crab::cfg::return_stmt<number_t, VariableName> return_t;
  typedef crab::cfg::intrinsic_stmt<number_t, VariableName> intrinsic_t;  

  typedef crab::cfg::array_init_stmt<number_t, VariableName> arr_init_t;
  typedef crab::cfg::array_store_stmt<number_t, VariableName> arr_store_t;
  typedef crab::cfg::array_load_stmt<number_t, VariableName> arr_load_t;
  typedef crab::cfg::array_assign_stmt<number_t, VariableName> arr_assign_t;

  typedef crab::cfg::region_init_stmt<number_t, varname_t> region_init_t;  
  typedef crab::cfg::make_ref_stmt<number_t, varname_t>  make_ref_t;
  typedef crab::cfg::load_from_ref_stmt<number_t, varname_t> load_from_ref_t;
  typedef crab::cfg::store_to_ref_stmt<number_t, varname_t> store_to_ref_t;
  typedef crab::cfg::gep_ref_stmt<number_t, varname_t> gep_ref_t;
  typedef crab::cfg::load_from_arr_ref_stmt<number_t, varname_t> load_from_arr_ref_t;
  typedef crab::cfg::store_to_arr_ref_stmt<number_t, varname_t> store_to_arr_ref_t;    
  typedef crab::cfg::assume_ref_stmt<number_t, VariableName> assume_ref_t;
  typedef crab::cfg::assert_ref_stmt<number_t, VariableName> assert_ref_t;

  typedef crab::cfg::bool_binary_op<number_t, VariableName> bool_bin_op_t;
  typedef crab::cfg::bool_assign_cst<number_t, VariableName> bool_assign_cst_t;
  typedef crab::cfg::bool_assign_var<number_t, VariableName> bool_assign_var_t;
  typedef crab::cfg::bool_assume_stmt<number_t, VariableName> bool_assume_t;
  typedef crab::cfg::bool_select_stmt<number_t, VariableName> bool_select_t;
  typedef crab::cfg::bool_assert_stmt<number_t, VariableName> bool_assert_t;

protected:
  virtual void exec(havoc_t &) {}
  virtual void exec(unreach_t &) {}
  virtual void exec(bin_op_t &) {}
  virtual void exec(assign_t &) {}
  virtual void exec(assume_t &) {}
  virtual void exec(select_t &) {}
  virtual void exec(assert_t &) {}
  virtual void exec(int_cast_t &) {}
  virtual void exec(callsite_t &) {}
  virtual void exec(return_t &) {}
  virtual void exec(intrinsic_t &) {}  
  virtual void exec(arr_init_t &) {}
  virtual void exec(arr_store_t &) {}
  virtual void exec(arr_load_t &) {}
  virtual void exec(arr_assign_t &) {}
  virtual void exec(region_init_t &) {}  
  virtual void exec(make_ref_t &) {}
  virtual void exec(load_from_ref_t &) {}
  virtual void exec(store_to_ref_t &) {}
  virtual void exec(gep_ref_t &) {}
  virtual void exec(load_from_arr_ref_t &) {}
  virtual void exec(store_to_arr_ref_t &) {}  
  virtual void exec(assume_ref_t &) {}
  virtual void exec(assert_ref_t &) {}
  virtual void exec(bool_bin_op_t &) {}
  virtual void exec(bool_assign_cst_t &) {}
  virtual void exec(bool_assign_var_t &) {}
  virtual void exec(bool_assume_t &) {}
  virtual void exec(bool_select_t &) {}
  virtual void exec(bool_assert_t &) {}

public: /* visitor api */
  void visit(havoc_t &s) { exec(s); }
  void visit(unreach_t &s) { exec(s); }
  void visit(bin_op_t &s) { exec(s); }
  void visit(assign_t &s) { exec(s); }
  void visit(assume_t &s) { exec(s); }
  void visit(select_t &s) { exec(s); }
  void visit(assert_t &s) { exec(s); }
  void visit(int_cast_t &s) { exec(s); }
  void visit(callsite_t &s) { exec(s); }
  void visit(return_t &s) { exec(s); }
  void visit(intrinsic_t &s) { exec(s); }  
  void visit(arr_init_t &s) { exec(s); }
  void visit(arr_store_t &s) { exec(s); }
  void visit(arr_load_t &s) { exec(s); }
  void visit(arr_assign_t &s) { exec(s); }
  void visit(region_init_t &s) { exec(s); }  
  void visit(make_ref_t &s) { exec(s); }
  void visit(load_from_ref_t &s) { exec(s); }
  void visit(store_to_ref_t &s) { exec(s); }
  void visit(gep_ref_t &s) { exec(s); }
  void visit(load_from_arr_ref_t &s) { exec(s); }
  void visit(store_to_arr_ref_t &s) { exec(s); }  
  void visit(assume_ref_t &s) { exec(s); }
  void visit(assert_ref_t &s) { exec(s); }
  void visit(bool_bin_op_t &s) { exec(s); }
  void visit(bool_assign_cst_t &s) { exec(s); }
  void visit(bool_assign_var_t &s) { exec(s); }
  void visit(bool_assume_t &s) { exec(s); }
  void visit(bool_select_t &s) { exec(s); }
  void visit(bool_assert_t &s) { exec(s); }
};


/**
 * Convert CFG operations into abstract domain operations
 **/

template <typename T>
inline boost::optional<T> conv_op(cfg::binary_operation_t op);
template <typename T>
inline boost::optional<T> conv_op(cfg::bool_binary_operation_t op);
template <typename T>
inline boost::optional<T> conv_op(cfg::cast_operation_t op);
  
template <>
inline boost::optional<domains::arith_operation_t>
conv_op(cfg::binary_operation_t op) {
  switch (op) {
  case cfg::BINOP_ADD: return domains::OP_ADDITION;
  case cfg::BINOP_SUB: return domains::OP_SUBTRACTION;
  case cfg::BINOP_MUL: return domains::OP_MULTIPLICATION;
  case cfg::BINOP_SDIV:return domains::OP_SDIV; 
  case cfg::BINOP_UDIV:return domains::OP_UDIV; 
  case cfg::BINOP_SREM:return domains::OP_SREM;
  case cfg::BINOP_UREM:return domains::OP_UREM;
  default:return boost::optional<domains::arith_operation_t>();
  }
}

template <>
inline boost::optional<domains::bitwise_operation_t>
conv_op(cfg::binary_operation_t op) {
  switch (op) {
  case cfg::BINOP_AND: return domains::OP_AND;
  case cfg::BINOP_OR:  return domains::OP_OR;
  case cfg::BINOP_XOR: return domains::OP_XOR;
  case cfg::BINOP_SHL: return domains::OP_SHL;
  case cfg::BINOP_LSHR:return domains::OP_LSHR; 
  case cfg::BINOP_ASHR:return domains::OP_ASHR;
  default: return boost::optional<domains::bitwise_operation_t>();
  }
}

template <>
inline boost::optional<domains::int_conv_operation_t>
conv_op(cfg::cast_operation_t op) {
  switch (op) {
  case cfg::CAST_TRUNC: return domains::OP_TRUNC;
  case cfg::CAST_SEXT:  return domains::OP_SEXT;
  case cfg::CAST_ZEXT:  return domains::OP_ZEXT;
  default: return boost::optional<domains::int_conv_operation_t>();         
  }
}

template <>
inline boost::optional<domains::bool_operation_t>
conv_op(cfg::bool_binary_operation_t op) {
  switch (op) {
  case cfg::BINOP_BAND: return domains::OP_BAND;
  case cfg::BINOP_BOR:  return domains::OP_BOR;
  case cfg::BINOP_BXOR: return domains::OP_BXOR;
  default: return boost::optional<domains::bool_operation_t>();
  }
}

  
/**
 * Abstract forward transformer for all statements. Function calls
 * can be redefined by derived classes. By default, all function
 * calls are ignored in a sound manner (by havoc'ing all outputs).
 **/
template <class AbsD>
class intra_abs_transformer
    : public abs_transformer_api<typename AbsD::number_t, typename AbsD::varname_t> {
public:
  typedef AbsD abs_dom_t;
  typedef typename abs_dom_t::number_t number_t;
  typedef typename abs_dom_t::varname_t varname_t;
  typedef typename abs_dom_t::variable_t variable_t;

public:
  typedef abs_transformer_api<number_t, varname_t> abs_transform_api_t;
  using typename abs_transform_api_t::arr_assign_t;
  using typename abs_transform_api_t::arr_init_t;
  using typename abs_transform_api_t::arr_load_t;
  using typename abs_transform_api_t::arr_store_t;
  using typename abs_transform_api_t::assert_t;
  using typename abs_transform_api_t::assign_t;
  using typename abs_transform_api_t::assume_t;
  using typename abs_transform_api_t::bin_op_t;
  using typename abs_transform_api_t::bool_assert_t;
  using typename abs_transform_api_t::bool_assign_cst_t;
  using typename abs_transform_api_t::bool_assign_var_t;
  using typename abs_transform_api_t::bool_assume_t;
  using typename abs_transform_api_t::bool_bin_op_t;
  using typename abs_transform_api_t::bool_select_t;
  using typename abs_transform_api_t::callsite_t;
  using typename abs_transform_api_t::intrinsic_t;  
  using typename abs_transform_api_t::havoc_t;
  using typename abs_transform_api_t::int_cast_t;
  using typename abs_transform_api_t::lin_cst_sys_t;
  using typename abs_transform_api_t::lin_cst_t;
  using typename abs_transform_api_t::lin_exp_t;
  using typename abs_transform_api_t::region_init_t;  
  using typename abs_transform_api_t::make_ref_t;
  using typename abs_transform_api_t::load_from_ref_t;
  using typename abs_transform_api_t::store_to_ref_t;
  using typename abs_transform_api_t::gep_ref_t;
  using typename abs_transform_api_t::load_from_arr_ref_t;
  using typename abs_transform_api_t::store_to_arr_ref_t;
  using typename abs_transform_api_t::assume_ref_t;    
  using typename abs_transform_api_t::assert_ref_t;
  using typename abs_transform_api_t::return_t;
  using typename abs_transform_api_t::select_t;
  using typename abs_transform_api_t::unreach_t;
  using typename abs_transform_api_t::var_t;

protected:
  abs_dom_t m_inv;
  bool m_ignore_assert;

private:
  template <typename NumOrVar>
  void apply(abs_dom_t &inv, cfg::binary_operation_t op,
	     const variable_t &x, const variable_t &y, NumOrVar z) {
    if (auto top = conv_op<domains::arith_operation_t>(op)) {
      inv.apply(*top, x, y, z);
    } else if (auto top = conv_op<domains::bitwise_operation_t>(op)) {
      inv.apply(*top, x, y, z);
    } else {
      CRAB_ERROR("unsupported binary operator", op);
    }
  }

public:
  intra_abs_transformer(abs_dom_t inv, bool ignore_assert = false)
      : m_inv(inv), m_ignore_assert(ignore_assert) {}

  virtual ~intra_abs_transformer() {}

  void set_abs_value(abs_dom_t &&inv) { m_inv = std::move(inv); }

  abs_dom_t get_abs_value() const { return m_inv;}
  abs_dom_t &get_abs_value() { return m_inv; }

  void exec(bin_op_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag &&
        (!(stmt.op() >= cfg::BINOP_SDIV && stmt.op() <= cfg::BINOP_UREM))) {
      pre_bot = m_inv.is_bottom();
    }

    const lin_exp_t &op1 = stmt.left();
    const lin_exp_t &op2 = stmt.right();
    if (op1.get_variable() && op2.get_variable()) {
      apply(m_inv, stmt.op(), stmt.lhs(), (*op1.get_variable()),
            (*op2.get_variable()));
    } else {
      assert(op1.get_variable() && op2.is_constant());
      apply(m_inv, stmt.op(), stmt.lhs(), (*op1.get_variable()),
            op2.constant());
    }

    if (::crab::CrabSanityCheckFlag &&
        (!(stmt.op() >= cfg::BINOP_SDIV && stmt.op() <= cfg::BINOP_UREM))) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
  }

  void exec(select_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    abs_dom_t inv1(m_inv);
    abs_dom_t inv2(m_inv);

    inv1 += stmt.cond();
    inv2 += stmt.cond().negate();

    if (::crab::CrabSanityCheckFlag) {
      if (!pre_bot && (inv1.is_bottom() && inv2.is_bottom())) {
        CRAB_ERROR(
            "select condition and its negation cannot be false simultaneously ",
            stmt);
      }
    }

    if (inv2.is_bottom()) {
      inv1.assign(stmt.lhs(), stmt.left());
      m_inv = inv1;
    } else if (inv1.is_bottom()) {
      inv2.assign(stmt.lhs(), stmt.right());
      m_inv = inv2;
    } else {
      inv1.assign(stmt.lhs(), stmt.left());
      inv2.assign(stmt.lhs(), stmt.right());
      m_inv = inv1 | inv2;
    }

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
  }

  void exec(assign_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.assign(stmt.lhs(), stmt.rhs());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
  }

  void exec(assume_t &stmt) { m_inv.operator+=(stmt.constraint()); }

  void exec(assert_t &stmt) {
    if (m_ignore_assert)
      return;

    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.operator+=(stmt.constraint());

    if (::crab::CrabSanityCheckFlag) {
      if (!stmt.constraint().is_contradiction()) {
        bool post_bot = m_inv.is_bottom();
        if (!(pre_bot || !post_bot)) {
          CRAB_WARN("Invariant became bottom after ", stmt, ".",
                    " This might indicate that the assertion is violated");
        }
      }
    }
  }

  void exec(int_cast_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }
    if (auto op = conv_op<crab::domains::int_conv_operation_t>(stmt.op())) {
      m_inv.apply(*op, stmt.dst(), stmt.src());
    } else {
      CRAB_ERROR("unsupported cast operator ", stmt.op());
    }

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
  }

  void exec(bool_assign_cst_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.assign_bool_cst(stmt.lhs(), stmt.rhs());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
  }

  void exec(bool_assign_var_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }
    m_inv.assign_bool_var(stmt.lhs(), stmt.rhs(), stmt.is_rhs_negated());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
  }

  void exec(bool_bin_op_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    if (auto op = conv_op<domains::bool_operation_t>(stmt.op())) {
      m_inv.apply_binary_bool(*op, stmt.lhs(), stmt.left(), stmt.right());
    } else {
      CRAB_WARN("Unsupported statement ", stmt);
    }

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
  }

  void exec(bool_assume_t &stmt) {
    m_inv.assume_bool(stmt.cond(), stmt.is_negated());
  }

  void exec(bool_select_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    abs_dom_t inv1(m_inv);
    abs_dom_t inv2(m_inv);
    const bool negate = true;
    inv1.assume_bool(stmt.cond(), !negate);
    inv2.assume_bool(stmt.cond(), negate);
    if (inv2.is_bottom()) {
      inv1.assign_bool_var(stmt.lhs(), stmt.left(), !negate);
      m_inv = inv1;
    } else if (inv1.is_bottom()) {
      inv2.assign_bool_var(stmt.lhs(), stmt.right(), !negate);
      m_inv = inv2;
    } else {
      inv1.assign_bool_var(stmt.lhs(), stmt.left(), !negate);
      inv2.assign_bool_var(stmt.lhs(), stmt.right(), !negate);
      m_inv = inv1 | inv2;
    }

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
  }

  void exec(bool_assert_t &stmt) {
    if (m_ignore_assert)
      return;

    m_inv.assume_bool(stmt.cond(), false);
  }

  void exec(havoc_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.operator-=(stmt.get_variable());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
  }

  void exec(unreach_t &stmt) { m_inv.set_to_bottom(); }

  void exec(arr_init_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.array_init(stmt.array(), stmt.elem_size(), stmt.lb_index(),
                     stmt.ub_index(), stmt.val());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
  }

  void exec(arr_store_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    if (stmt.lb_index().equal(stmt.ub_index())) {
      m_inv.array_store(stmt.array(), stmt.elem_size(), stmt.lb_index(),
			stmt.value(), stmt.is_strong_update());
    } else {
      m_inv.array_store_range(stmt.array(), stmt.elem_size(), stmt.lb_index(),
			      stmt.ub_index(), stmt.value());
    }

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
  }

  void exec(arr_load_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.array_load(stmt.lhs(), stmt.array(), stmt.elem_size(), stmt.index());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
  }

  void exec(arr_assign_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.array_assign(stmt.lhs(), stmt.rhs());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
  }

  void exec(make_ref_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }
    
    m_inv.ref_make(stmt.lhs(), stmt.region());
    
    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }    
  }

  void exec(region_init_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }
    
    m_inv.region_init(stmt.region());
    
    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }    
  }

  void exec(load_from_ref_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }
    
    m_inv.ref_load(stmt.ref(), stmt.region(), stmt.lhs());
    
    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }    
  }
  
  void exec(store_to_ref_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }
    
    m_inv.ref_store(stmt.ref(), stmt.region(), stmt.val());
    
    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }        
  }
  
  void exec(gep_ref_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }
    
    m_inv.ref_gep(stmt.rhs(), stmt.rhs_region(), stmt.lhs(), stmt.lhs_region(), stmt.offset());
    
    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }        
  }
  
  void exec(load_from_arr_ref_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }
    
    m_inv.ref_load_from_array(stmt.lhs(), stmt.ref(), stmt.region(),
			      stmt.index(), stmt.elem_size());
    
    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }        
  }
  
  void exec(store_to_arr_ref_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    if (stmt.lb_index().equal(stmt.ub_index())) {
      m_inv.ref_store_to_array(stmt.ref(), stmt.region(),
			       stmt.lb_index(), stmt.elem_size(),
			       stmt.value());
    } else {
      CRAB_ERROR("TODO store_to_array_ref for ranges");
    }
    
    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }        
  }  
  
  void exec(assume_ref_t &stmt) {
    m_inv.ref_assume(stmt.constraint());
  }

  void exec(assert_ref_t &stmt) {
    if (m_ignore_assert)
      return;
    m_inv.ref_assume(stmt.constraint());
  }

  void exec(intrinsic_t &cs) {
    m_inv.intrinsic(cs.get_intrinsic_name(), cs.get_args(), cs.get_lhs());
  }
  
  virtual void exec(callsite_t &cs) {
    for (const variable_t &vt: cs.get_lhs()) {
      m_inv.operator-=(vt); // havoc
    }
  }
  
  virtual void exec(return_t &ret) {}
  
};

///////////////////////////////////////
/// For inter-procedural analysis
///////////////////////////////////////

template <typename AbsDom> class inter_transformer_helpers {
public:
  typedef typename AbsDom::linear_expression_t linear_expression_t;
  typedef typename AbsDom::reference_constraint_t reference_constraint_t;  
  typedef typename AbsDom::variable_t variable_t;
  typedef typename AbsDom::number_t number_t;

  static void unify(AbsDom &inv, const variable_t &lhs, const variable_t &rhs) {
    assert(lhs.get_type() == rhs.get_type());
    switch (lhs.get_type()) {
    case BOOL_TYPE:
      inv.assign_bool_var(lhs, rhs, false);
      break;
    case INT_TYPE:
    case REAL_TYPE:
      inv.assign(lhs, rhs);
      break;
    case REF_TYPE:
      inv -= lhs;
      inv.ref_assume(reference_constraint_t::mk_eq(lhs, rhs, number_t(0)));
      break;
    case ARR_BOOL_TYPE:
    case ARR_INT_TYPE:
    case ARR_REAL_TYPE:
      inv.array_assign(lhs, rhs);
      break;
    default:
      CRAB_ERROR("unsuported type");
    }
  }
};

/////////////////////////////////
/// For backward analysis
/////////////////////////////////

/**
 * Abstract transformer to compute necessary preconditions.
 **/
template <class AbsD, class InvT>
class intra_necessary_preconditions_abs_transformer final
    : public abs_transformer_api<typename AbsD::number_t,
                                 typename AbsD::varname_t> {
public:
  typedef AbsD abs_dom_t;
  typedef typename abs_dom_t::number_t number_t;
  typedef typename abs_dom_t::varname_t varname_t;
  typedef typename abs_dom_t::variable_t variable_t;
  typedef crab::cfg::statement<number_t, varname_t> statement_t;
  typedef abs_transformer_api<number_t, varname_t> abs_transform_api_t;
  using typename abs_transform_api_t::arr_assign_t;
  using typename abs_transform_api_t::arr_init_t;
  using typename abs_transform_api_t::arr_load_t;
  using typename abs_transform_api_t::arr_store_t;
  using typename abs_transform_api_t::assert_t;
  using typename abs_transform_api_t::assign_t;
  using typename abs_transform_api_t::assume_t;
  using typename abs_transform_api_t::bin_op_t;
  using typename abs_transform_api_t::bool_assert_t;
  using typename abs_transform_api_t::bool_assign_cst_t;
  using typename abs_transform_api_t::bool_assign_var_t;
  using typename abs_transform_api_t::bool_assume_t;
  using typename abs_transform_api_t::bool_bin_op_t;
  using typename abs_transform_api_t::bool_select_t;
  using typename abs_transform_api_t::callsite_t;
  using typename abs_transform_api_t::intrinsic_t;  
  using typename abs_transform_api_t::havoc_t;
  using typename abs_transform_api_t::int_cast_t;
  using typename abs_transform_api_t::lin_cst_sys_t;
  using typename abs_transform_api_t::lin_cst_t;
  using typename abs_transform_api_t::lin_exp_t;
  using typename abs_transform_api_t::make_ref_t;
  using typename abs_transform_api_t::region_init_t;  
  using typename abs_transform_api_t::load_from_ref_t;
  using typename abs_transform_api_t::store_to_ref_t;
  using typename abs_transform_api_t::gep_ref_t;
  using typename abs_transform_api_t::load_from_arr_ref_t;
  using typename abs_transform_api_t::store_to_arr_ref_t;
  using typename abs_transform_api_t::assume_ref_t;  
  using typename abs_transform_api_t::assert_ref_t;
  using typename abs_transform_api_t::return_t;
  using typename abs_transform_api_t::select_t;
  using typename abs_transform_api_t::unreach_t;
  using typename abs_transform_api_t::var_t;

private:
  // used to compute the (necessary) preconditions
  abs_dom_t m_pre;
  // used to refine the preconditions: map from statement_t to abs_dom_t.
  InvT *m_invariants;
  // ignore assertions
  bool m_ignore_assert;
  // if m_ignore_assert is false then enable compute preconditions
  // from good states, otherwise from bad states (by negating the
  // conditions of the assert statements).
  bool m_good_states;

  abs_dom_t make_top() const {
    return m_pre.make_top();
  }
  
  abs_dom_t get_forward_invariant(const statement_t *stmt) const{
    assert(m_invariants);
    assert(stmt);
    
    auto it = m_invariants->find(stmt);
    if (it != m_invariants->end()) {
      return it->second;
    } else {
      return make_top();
    }
  }
  
  
public:
  intra_necessary_preconditions_abs_transformer(abs_dom_t post, InvT *invars,
                                                bool good_states,
                                                bool ignore_assert = false)
      : m_pre(post), m_invariants(invars), m_ignore_assert(ignore_assert),
        m_good_states(good_states) {

    if (!m_invariants) {
      CRAB_ERROR("Invariant table cannot be null");
    }
  }

  ~intra_necessary_preconditions_abs_transformer() = default;

  abs_dom_t preconditions() { return m_pre; }

  void exec(bin_op_t &stmt) {
    auto op = conv_op<domains::arith_operation_t>(stmt.op());
    if (!op || op >= domains::OP_UDIV) {
      // ignore UDIV, SREM, UREM
      // CRAB_WARN("backward operation ", stmt.op(), " not implemented");
      m_pre -= stmt.lhs();
      return;
    }

    const lin_exp_t &op1 = stmt.left();
    const lin_exp_t &op2 = stmt.right();
    abs_dom_t invariant = get_forward_invariant(&stmt);

    CRAB_LOG("backward-tr", crab::outs()
                                << "** " << stmt.lhs() << " := " << op1 << " "
                                << *op << " " << op2 << "\n"
                                << "\tFORWARD INV=" << invariant << "\n"
                                << "\tPOST=" << m_pre << "\n");

    if (op1.get_variable() && op2.get_variable()) {
      m_pre.backward_apply(*op, stmt.lhs(), (*op1.get_variable()),
                           (*op2.get_variable()), std::move(invariant));
    } else {
      assert(op1.get_variable() && op2.is_constant());
      m_pre.backward_apply(*op, stmt.lhs(), (*op1.get_variable()),
                           op2.constant(), std::move(invariant));
    }
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  // select(x := cond ? e1: e2, post) can be reduced to
  //   pre: goto b_then;
  //   pre: goto b_else;
  //   b_then:
  //     assume(cond);
  //     x := e1;
  //     goto post;
  //   b_else:
  //     assume(not(cond));
  //     x := e2;
  //     goto post;
  //   post: ....
  void exec(select_t &stmt) {
    abs_dom_t old_pre = get_forward_invariant(&stmt);

    // -- one of the two branches is false
    abs_dom_t then_inv(old_pre);
    then_inv += stmt.cond();
    if (then_inv.is_bottom()) {
      m_pre.backward_assign(stmt.lhs(), stmt.right(), std::move(old_pre));
      m_pre += stmt.cond().negate();
      return;
    }

    abs_dom_t else_inv(old_pre);
    else_inv += stmt.cond().negate();
    if (else_inv.is_bottom()) {
      m_pre.backward_assign(stmt.lhs(), stmt.left(), std::move(old_pre));
      m_pre += stmt.cond();
      return;
    }

    // -- both branches can be possible so we join them
    abs_dom_t pre_then(m_pre);
    pre_then.backward_assign(stmt.lhs(), stmt.left(), old_pre);
    pre_then += stmt.cond();

    abs_dom_t pre_else(m_pre);
    pre_else.backward_assign(stmt.lhs(), stmt.right(), old_pre);
    pre_else += stmt.cond().negate();

    m_pre = pre_then | pre_else;
  }

  // x := e
  void exec(assign_t &stmt) {
    abs_dom_t invariant = get_forward_invariant(&stmt);

    CRAB_LOG("backward-tr", auto rhs = stmt.rhs();
             crab::outs() << "** " << stmt.lhs() << " := " << rhs << "\n"
                          << "\tFORWARD INV=" << invariant << "\n"
                          << "\tPOST=" << m_pre << "\n");

    m_pre.backward_assign(stmt.lhs(), stmt.rhs(), std::move(invariant));
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  // assume(c)
  // the precondition must contain c so forward and backward are the same.
  void exec(assume_t &stmt) {
    CRAB_LOG("backward-tr", crab::outs() << "** " << stmt << "\n"
                                         << "\tPOST=" << m_pre << "\n");
    m_pre += stmt.constraint();
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  // assert(c)
  void exec(assert_t &stmt) {
    if (!m_ignore_assert) {
      CRAB_LOG("backward-tr", crab::outs() << "** " << stmt << "\n"
                                           << "\tPOST=" << m_pre << "\n");
      if (m_good_states) {
        // similar to assume(c)
        m_pre += stmt.constraint();
      } else {
        // here we are interested in computing preconditions of the
        // error states. Thus, we propagate backwards "not c" which
        // represents the error states.
        abs_dom_t error = m_pre.make_top();
        error += stmt.constraint().negate();
        // we need to join to consider all possible preconditions to
        // error. Otherwise, if we would have two assertions
        // "assert(x >= -2); assert(x <= 2);" we would have
        // incorrectly contradictory constraints.
        m_pre |= error;
      }

      CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
    }
  }

  // similar to assume(false)
  void exec(unreach_t &stmt) { m_pre.set_to_bottom(); }

  // x := *
  // x can be anything before the assignment
  void exec(havoc_t &stmt) { m_pre -= stmt.get_variable(); }

  void exec(int_cast_t &stmt) {
    abs_dom_t invariant = get_forward_invariant(&stmt);
    CRAB_LOG("backward-tr", crab::outs() << "** " << stmt << "\n"
                                         << "\tPOST=" << m_pre << "\n");
    m_pre.backward_assign(stmt.dst(), stmt.src(), std::move(invariant));
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  void exec(bool_assign_cst_t &stmt) { m_pre -= stmt.lhs(); }
  void exec(bool_assign_var_t &stmt) { m_pre -= stmt.lhs(); }
  void exec(bool_bin_op_t &stmt) { m_pre -= stmt.lhs(); }
  void exec(bool_select_t &stmt) { m_pre -= stmt.lhs(); }

  void exec(bool_assume_t &stmt) {
    // same as forward
    CRAB_LOG("backward-tr", crab::outs() << "** " << stmt << "\n"
                                         << "\tPOST=" << m_pre << "\n");
    m_pre.assume_bool(stmt.cond(), stmt.is_negated());
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  void exec(bool_assert_t &stmt) {
    if (!m_ignore_assert) {
      CRAB_LOG("backward-tr", crab::outs() << "** " << stmt << "\n"
                                           << "\tPOST=" << m_pre << "\n");
      if (m_good_states) {
        // similar to bool_assume(c)
        m_pre.assume_bool(stmt.cond(), false /*non-negated*/);
      } else {
        abs_dom_t error = m_pre.make_top();
        error.assume_bool(stmt.cond(), true /*negated*/);
        m_pre |= error;
      }
      CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
    }
  }

  void exec(arr_init_t &stmt) {
    abs_dom_t invariant = get_forward_invariant(&stmt);

    CRAB_LOG("backward-tr", crab::outs()
                                << "** " << stmt << "\n"
                                << "\tFORWARD INV=" << invariant << "\n"
                                << "\tPOST=" << m_pre << "\n");
    m_pre.backward_array_init(stmt.array(), stmt.elem_size(), stmt.lb_index(),
                              stmt.ub_index(), stmt.val(),
                              std::move(invariant));
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  void exec(arr_load_t &stmt) {
    abs_dom_t invariant = get_forward_invariant(&stmt);

    CRAB_LOG("backward-tr", crab::outs()
                                << "** " << stmt << "\n"
                                << "\tFORWARD INV=" << invariant << "\n"
                                << "\tPOST=" << m_pre << "\n");
    m_pre.backward_array_load(stmt.lhs(), stmt.array(), stmt.elem_size(),
                              stmt.index(), std::move(invariant));
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  void exec(arr_store_t &stmt) {
    abs_dom_t invariant = get_forward_invariant(&stmt);
    CRAB_LOG("backward-tr", crab::outs()
                                << "** " << stmt << "\n"
                                << "\tFORWARD INV=" << invariant << "\n"
                                << "\tPOST=" << m_pre << "\n");

    if (stmt.lb_index().equal(stmt.ub_index())) {
      m_pre.backward_array_store(
            stmt.array(), stmt.elem_size(), stmt.lb_index(), stmt.value(),
            stmt.is_strong_update(), std::move(invariant));
    } else {
      m_pre.backward_array_store_range(stmt.array(), stmt.elem_size(),
				       stmt.lb_index(), stmt.ub_index(),
				       stmt.value(), std::move(invariant));
    }
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  void exec(arr_assign_t &stmt) {
    abs_dom_t invariant = get_forward_invariant(&stmt);
    CRAB_LOG("backward-tr", crab::outs()
                                << "** " << stmt << "\n"
                                << "\tFORWARD INV=" << invariant << "\n"
                                << "\tPOST=" << m_pre << "\n");
    m_pre.backward_array_assign(stmt.lhs(), stmt.rhs(), std::move(invariant));
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  // NOT IMPLEMENTED
  void exec(region_init_t &stmt) {}  
  void exec(make_ref_t &stmt) {}
  void exec(load_from_ref_t &stmt) {}
  void exec(store_to_ref_t &stmt) {}
  void exec(gep_ref_t &stmt) {}
  void exec(load_from_arr_ref_t &stmt) {}
  void exec(store_to_arr_ref_t &stmt) {}
  void exec(assume_ref_t &stmt) {}
  void exec(assert_ref_t &stmt) {}
  
  /// -- Call and return can be redefined by derived classes

  virtual void exec(callsite_t &cs) {
    for (const variable_t &vt : cs.get_lhs()) {
      m_pre -= vt;
    }
  }
  virtual void exec(return_t &stmt) {}

  void exec(intrinsic_t &cs) {
    abs_dom_t invariant = get_forward_invariant(&cs);    
    m_pre.backward_intrinsic(cs.get_intrinsic_name(), cs.get_args(), cs.get_lhs(),
			     std::move(invariant));    
  }
  
};

} // namespace analyzer
} // namespace crab
