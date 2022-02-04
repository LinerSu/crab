#pragma once

#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/flat_boolean_domain.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/domains/small_range.hpp>
#include <crab/domains/types.hpp>
#include <crab/domains/uf_domain.hpp>
#include <crab/domains/union_find_domain.hpp>
#include <crab/types/tag.hpp>
#include <crab/types/varname_factory.hpp>

#include "region/ghost_variables.hpp"
#include "region/tags.hpp"

#include <unordered_map>

namespace crab {
namespace domains {
namespace mru_region_domain_impl {
template <class Number, class VariableName, class BaseAbsDom> class Params {
public:
  using number_t = Number;
  using varname_t = VariableName;
  using varname_allocator_t = crab::var_factory_impl::str_var_alloc_col;
  using base_abstract_domain_t = BaseAbsDom;
  using base_varname_t = typename BaseAbsDom::varname_t;

  static_assert(std::is_same<Number, typename BaseAbsDom::number_t>::value,
                "Number type and BaseAbsDom::number_t must be the same");
  // This is a strong requirement
  static_assert(
      std::is_same<base_varname_t,
                   typename varname_allocator_t::varname_t>::value,
      "BaseAbsDom::varname_t and allocator_varname_t must be the same");
};
} // namespace mru_region_domain_impl

template <typename Params>
class mru_region_domain
    : public abstract_domain_api<mru_region_domain<Params>> {
  using mru_region_domain_t = mru_region_domain<Params>;
  using abstract_domain_t = abstract_domain_api<mru_region_domain_t>;

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
  using number_t = typename Params::number_t;
  using varname_t = typename Params::varname_t;

private:
  using base_varname_t = typename Params::base_varname_t;
  using base_uf_domain_t = uf_domain<number_t, base_varname_t>;
  using base_abstract_domain_t = typename Params::base_abstract_domain_t;
  using base_variable_vector_t =
      typename base_abstract_domain_t::variable_vector_t;
  using base_variable_t = typename base_abstract_domain_t::variable_t;
  using base_variable_or_constant_t =
      typename base_abstract_domain_t::variable_or_constant_t;
  using base_linear_expression_t =
      typename base_abstract_domain_t::linear_expression_t;
  using base_linear_constraint_t =
      typename base_abstract_domain_t::linear_constraint_t;
  using base_linear_constraint_system_t =
      typename base_abstract_domain_t::linear_constraint_system_t;
  using base_varname_allocator_t = typename Params::varname_allocator_t;

  /**------------------ Begin type definitions ------------------**/
  using ghost_variables_t = region_domain_impl::ghost_variables<
      abstract_domain_t, base_abstract_domain_t, base_varname_allocator_t>;
  using tag_t = region_domain_impl::tag<number_t>;
  // Environment domains: map regions to finite domain
  using rgn_counting_env_t = ikos::separate_domain<variable_t, small_range>;
  using rgn_bool_env_t = ikos::separate_domain<variable_t, boolean_value>;
  using rgn_type_env_t = ikos::separate_domain<variable_t, type_value>;
  // Union-find where equivalence classes are attached to boolean values
  using rgn_dealloc_t = union_find_domain<variable_t, boolean_value>;
  // Map from variable to its ghost variables
  using var_map_t = std::unordered_map<variable_t, ghost_variables_t>;
  // Reverse map used only for pretty printing: a CrabIR variable has
  // associated a set of ghost variables. This reverse map maps back
  // each ghost variable to its CrabIR variable.
  using ghost_var_id = unsigned;
  using rev_var_map_t =
      std::unordered_map<base_variable_t, std::pair<variable_t, ghost_var_id>>;
  // Map from reference variable to a representation of base address
  // So far the uf domain is used for address domain.
  // The base address is an uniterpreted symbol in uf domain.
  // TODO: we need to change if we introduce ghost variable in address dom
  using uninterpreted_symbol_t = term::term_operator_t;
  using uf_var_map_t = union_find_domain<variable_t, uninterpreted_symbol_t>;
  using base_dom_binop_t = std::function<base_abstract_domain_t(
      base_abstract_domain_t, base_abstract_domain_t)>;
  using base_uf_dom_binop_t =
      std::function<base_uf_domain_t(base_uf_domain_t, base_uf_domain_t)>;
  // Union-find where equivalence classes are attached to boolean values
  // Each equivalence class represents a set of regions belonging to an object
  using rgn_equiv_t = union_find_domain<variable_t, boolean_value>;
  // Map each reference or region variable to a set of allocation sites
  using alloc_site_env_t =
      separate_discrete_domain<variable_t, allocation_site>;
  using allocation_sites = typename alloc_site_env_t::value_type;
  // Map variables to sets of tags
  using tag_env_t = separate_discrete_domain<variable_t, tag_t>;
  using c_base_addr_t =
      boost::optional<std::pair<ghost_variables_t, uninterpreted_symbol_t>>;
  using tag_set = typename tag_env_t::value_type;
  /**------------------ End type definitions ------------------**/

  /**------------------ Begin field definitions ------------------**/
  /*
      abs_domain mem;
      abs_domain cache_lines;
      uf_domain regs;
      uf_domain addrs;
    The implementation follows the most recently used memory model including
    cache and memory. The memory in our abstraction is represented as a domain
    and cache is a structure includeing cache lines, registers, and addresses.
    See the following comments for details:
  */

  // a special symbol to represent bottom state of the mru region domain
  // mru region domain is bottom iff mem is bottom.
  bool m_is_bottom;

  // To create ghost variables for each abstract subdomain.
  base_varname_allocator_t m_alloc;
  // Map a variable_t to its ghost variables:
  var_map_t m_var_map;
  // Reverse map from ghost variables to variable_t
  rev_var_map_t m_rev_var_map;
  // Map a region variable to its ghost variables
  // This map is used to track regions within the mru object
  var_map_t m_cache_var_map;

public:
  // The base abstract domains (over ghost variables): the domains only
  // contain variables, there is no variable of region / reference type.
  // NOTE: For each subdomain, all base variables remian as the same.
  // Any other operations do not require renaming, such as:
  //  - cache lines and regs reduction;
  //  - join / meet over cache lines and memory

  // an abstract domain that contains all the memory and register state.
  base_abstract_domain_t m_mem;
  // an abstract domain that contains all properties for the cache lines.
  // limitation:
  //    - cache is used only for one object
  base_abstract_domain_t m_cache_lines;
  // an abstract domain that models equalities between scalar variables
  // in m_mem and variables in m_cache_lines
  base_uf_domain_t m_regs;
  // an abstract domain that keeps equalities of base address between references
  // we assume each object is allocated by some address called base address.
  // each reference refers to an object so each reference has additional info
  // about base address. This domain tracks all base address alloacted so far.
  // Thus, two variables are mapped to the same uninterpreted symbol
  // if they belong to the same memory object.
  base_uf_domain_t m_addrs;

  // a ghost variable models the base address showing the cache refers to.
  // If cache is not used, this variable is not meaningful;
  // This one is updated once cache is used.
  // This variable is not presented in the Crab IR.
  // The name (C_base) is used for printing.
  c_base_addr_t m_base;

  // a 3-valued boolean. if true then cache is used.
  // if false then cache is not used. Otherwise, we do not know
  boolean_value m_used;
  // a 3-valued boolean. if true then cache has been updated.
  // if false then cache is not updated. Otherwise, we do not know
  boolean_value m_dirty;

private:
  // Abstract domain to count how many addresses are in a region.
  // This allows us to decide when strong update is sound: only if
  // one address per region (i.e., singleton).
  rgn_counting_env_t m_rgn_counting_dom;
  // Whether some data might have been written to any address within
  // the region.
  rgn_bool_env_t m_rgn_init_dom;
  // For each unknown region we keep track of the (dynamic) type of
  // its last written value.
  rgn_type_env_t m_rgn_type_dom;
  // Keep track of whether some memory within a region has been
  // deallocated.
  //
  // To reason about deallocation, we need to know which regions might
  // belong to the same allocated memory object.  Each call to
  // ref_make models a new allocation returning a reference to the
  // base address of the allocated block. Then, ref_gep is used to
  // perform pointer arithmetic within an allocated block. Very
  // importantly, ref_gep can switch between regions although we
  // assume that those regions always belong to the same allocated
  // block.
  //
  // We partition region variables into equivalence classes attached
  // to a boolean value. Two region variables are in the same
  // equivalence class if they might belong to the same allocated
  // block. The boolean flag indicates whether the allocated block
  // might have been deallocated.
  //
  // The partitioning is done as follows:
  //
  //   region_init(rgn): create a singleton equivalence class with rgn.
  //   ref_gep(reg1, rgn1, ref2, rgn2, o): join together the
  //                                       equivalence classes of rgn1
  //                                       and rgn2.
  rgn_dealloc_t m_rgn_dealloc_dom;
  // Tag analysis: map each (any type) variable to a set of tags
  tag_env_t m_tag_env;

public:
  // To reason about MRU region-based memory model. We need to know which
  // regions belong to which object. This info is passed by sea-dsa pointer
  // analysis. This map tracks all objects' regions used in the CrabIR.
  rgn_equiv_t m_obj_rgn_map;

private:
  uf_var_map_t m_addr_var_map;
  // rev_var_map_t m_rev_addr_var_map;
  int m_uninterpreted_symbol_index =
      term::term_operator_t::first_nonreserved_value();
  /**------------------ End field definitions ------------------**/

  /**------------------ Begin helper method definitions ------------------**/
  mru_region_domain(
      base_varname_allocator_t &&alloc, var_map_t &&var_map,
      rev_var_map_t &&rev_var_map, var_map_t &&cache_var_map,
      base_abstract_domain_t &&mem, base_abstract_domain_t &&cache_lines,
      base_uf_domain_t &&regs, base_uf_domain_t &&addrs, c_base_addr_t &&base,
      boolean_value &&used, boolean_value &&dirty,
      rgn_counting_env_t &&rgn_counting_dom, rgn_bool_env_t &&rgn_init_dom,
      rgn_type_env_t &&rgn_type_dom, rgn_dealloc_t &&rgn_dealloc_dom,
      tag_env_t &&tag_env, rgn_equiv_t &&obj_rgn_map,
      uf_var_map_t &&addr_var_map, int uninterpreted_symbol_index)
      : m_is_bottom(mem.is_bottom()), m_alloc(std::move(alloc)),
        m_var_map(std::move(var_map)), m_rev_var_map(std::move(rev_var_map)),
        m_cache_var_map(std::move(cache_var_map)), m_mem(std::move(mem)),
        m_cache_lines(std::move(cache_lines)), m_regs(std::move(regs)),
        m_addrs(std::move(addrs)), m_base(std::move(base)),
        m_used(std::move(used)), m_dirty(std::move(dirty)),
        m_rgn_counting_dom(std::move(rgn_counting_dom)),
        m_rgn_init_dom(std::move(rgn_init_dom)),
        m_rgn_type_dom(std::move(rgn_type_dom)),
        m_rgn_dealloc_dom(std::move(rgn_dealloc_dom)),
        m_tag_env(std::move(tag_env)), m_obj_rgn_map(std::move(obj_rgn_map)),
        m_addr_var_map(std::move(addr_var_map)),
        m_uninterpreted_symbol_index(m_uninterpreted_symbol_index) {}

  // Perform *this = join(*this, right)
  void do_join(const mru_region_domain_t &right_const) {

    mru_region_domain_t right = mru_region_domain_t(right_const);
    // Check whether need to perform cache reduction
    if (!(m_used == boolean_value::get_false()) &&
        !(right.m_used == boolean_value::get_false())) { // if both dirty
      // left.m_base == right.m_base if refer to the same uniterpreted symbol
      if (!((*m_base).second == (*right.m_base).second)) {
        commit_cache();
        right.commit_cache();
      }
    } else if (!(m_used == boolean_value::get_false())) {
      // left: perform reduction before join
      commit_cache();
    } else if (!(right.m_used == boolean_value::get_false())) {
      // right: perform reduction before join
      right.commit_cache();
    }

    // At this point, the left and right caches are
    // (1) empty, (2) dirty but refer to the same object

    // The following domains do not require common renaming
    // (i.e. no ghost variables). The domains are finite.
    rgn_counting_env_t out_rgn_counting_dom(m_rgn_counting_dom |
                                            right.m_rgn_counting_dom);
    rgn_bool_env_t out_rgn_init_dom(m_rgn_init_dom | right.m_rgn_init_dom);
    rgn_type_env_t out_rgn_type_dom(m_rgn_type_dom | right.m_rgn_type_dom);
    rgn_dealloc_t out_rgn_dealloc_dom(m_rgn_dealloc_dom |
                                      right.m_rgn_dealloc_dom);
    tag_env_t out_tag_env(m_tag_env | right.m_tag_env);
    rgn_equiv_t out_obj_rgn_map(m_obj_rgn_map | right.m_obj_rgn_map);
    uf_var_map_t out_addr_var_map(m_addr_var_map | right.m_addr_var_map);

    // Merge allocator
    base_varname_allocator_t out_alloc(m_alloc, right.m_alloc);

    boolean_value out_used(m_used | right.m_used);
    boolean_value out_dirty(m_dirty | right.m_dirty);

    // Perform renaming before join
    base_abstract_domain_t right_mem(right.m_mem);
    base_abstract_domain_t right_cache_lines(right.m_cache_lines);
    base_uf_domain_t right_regs(right.m_regs);
    base_uf_domain_t right_addrs(right.m_addrs);

    var_map_t out_var_map, out_cache_var_map;
    rev_var_map_t out_rev_var_map;
    // -- Compute common renamings
    base_variable_vector_t left_vars, right_vars, out_vars;
    // reserve vec space by upper bound to avoid reallocations
    size_t num_renamings = m_var_map.size();
    left_vars.reserve(num_renamings);
    right_vars.reserve(num_renamings);
    out_vars.reserve(num_renamings);
    // perform the common renamings
    for (auto &kv : m_var_map) {
      const variable_t &v = kv.first;
      auto it = right.m_var_map.find(v);
      if (it != right.m_var_map.end()) {
        // if we ignore unknown regions, types are the same by construction
        if (!crab_domain_params_man::get().region_skip_unknown_regions() &&
            !kv.second.same_type(it->second)) {
          crab::CrabStats::count(
              domain_name() + ".count.join.skipped.inconsistent_dynamic_type");
          continue;
        }
        ghost_variables_t out_gvars(
            ghost_variables_t::create(out_alloc, kv.second));
        // Store all ghost variables into each vector
        kv.second.add(left_vars);
        it->second.add(right_vars);
        out_gvars.add(out_vars);
        out_var_map.insert({v, out_gvars});
        out_gvars.update_rev_varmap(out_rev_var_map, v);

        // update cache used var map
        if (m_cache_var_map.find(v) != m_cache_var_map.end()) {
          if (right.m_cache_var_map.find(v) != right.m_cache_var_map.end()) {
            out_cache_var_map.insert({v, out_gvars});
          }
        }
      }
    }

    if (!(m_used == boolean_value::get_false()) &&
        !(right.m_used == boolean_value::get_false())) {
      assert(m_base);
      assert(right.m_base);
      if (!(*m_base).first.same_type((*right.m_base).first)) {
        // the type of ghost variable for base address should be the same
        CRAB_ERROR(
            domain_name(),
            "::do_join current and right should contain the same m_base type");
      }
      // base address refers to the same object
      (*m_base).first.add(left_vars);
      (*right.m_base).first.add(right_vars);
      (*m_base).first.add(out_vars);
    }

    // might need project the variables used in common because common variables
    // are subset of the variables used in domain
    m_mem.project(left_vars);
    m_mem.rename(left_vars, out_vars);
    right_mem.project(right_vars);
    right_mem.rename(right_vars, out_vars);

    m_cache_lines.project(left_vars);
    m_cache_lines.rename(left_vars, out_vars);
    right_cache_lines.project(right_vars);
    right_cache_lines.rename(right_vars, out_vars);

    m_regs.project(left_vars);
    m_regs.rename(left_vars, out_vars);
    right_regs.project(right_vars);
    right_regs.rename(right_vars, out_vars);

    m_addrs.project(left_vars);
    m_addrs.rename(left_vars, out_vars);
    right_addrs.project(right_vars);
    right_addrs.rename(right_vars, out_vars);

    m_mem |= right_mem;
    m_cache_lines |= right_cache_lines;
    m_regs |= right_regs;
    m_addrs |= right_addrs;
    m_is_bottom = m_mem.is_bottom();
    std::swap(m_alloc, out_alloc);
    std::swap(m_var_map, out_var_map);
    std::swap(m_rev_var_map, out_rev_var_map);
    std::swap(m_cache_var_map, out_cache_var_map);
    std::swap(m_used, out_used);
    std::swap(m_dirty, out_dirty);
    std::swap(m_rgn_counting_dom, out_rgn_counting_dom);
    std::swap(m_rgn_init_dom, out_rgn_init_dom);
    std::swap(m_rgn_type_dom, out_rgn_type_dom);
    std::swap(m_rgn_dealloc_dom, out_rgn_dealloc_dom);
    std::swap(m_tag_env, out_tag_env);
    std::swap(m_obj_rgn_map, out_obj_rgn_map);
    std::swap(m_addr_var_map, out_addr_var_map);
    m_uninterpreted_symbol_index = std::max(m_uninterpreted_symbol_index,
                                            right.m_uninterpreted_symbol_index);

    CRAB_LOG("mru-region", crab::outs() << *this << "\n");
  }

  mru_region_domain_t do_join_or_widening(
      const mru_region_domain_t &left_const,
      const mru_region_domain_t &right_const, const bool is_join /*unused*/,
      base_dom_binop_t base_dom_op, base_uf_dom_binop_t base_uf_dom_op) const {
    mru_region_domain_t left = mru_region_domain_t(left_const);
    mru_region_domain_t right = mru_region_domain_t(right_const);
    // Check whether need to perform cache reduction
    if (!(left.m_used == boolean_value::get_false()) &&
        !(right.m_used == boolean_value::get_false())) { // if both dirty
      // left.m_base == right.m_base if refer to the same uniterpreted symbol
      if (!((*left.m_base).second == (*right.m_base).second)) {
        left.commit_cache();
        right.commit_cache();
      }
    } else if (!(left.m_used == boolean_value::get_false())) {
      // left: perform reduction before join
      left.commit_cache();
    } else if (!(right.m_used == boolean_value::get_false())) {
      // right: perform reduction before join
      right.commit_cache();
    }

    // crab::outs() << "does left commit?" << left <<"\n";

    // At this point, the left and right caches are
    // (1) empty, (2) dirty but refer to the same object

    // The following domains do not require common renaming
    // (i.e. no ghost variables). The domains are finite.
    rgn_counting_env_t out_rgn_counting_dom(left.m_rgn_counting_dom |
                                            right.m_rgn_counting_dom);
    rgn_bool_env_t out_rgn_init_dom(left.m_rgn_init_dom | right.m_rgn_init_dom);
    rgn_type_env_t out_rgn_type_dom(left.m_rgn_type_dom | right.m_rgn_type_dom);
    rgn_dealloc_t out_rgn_dealloc_dom(left.m_rgn_dealloc_dom |
                                      right.m_rgn_dealloc_dom);
    tag_env_t out_tag_env(left.m_tag_env | right.m_tag_env);
    rgn_equiv_t out_obj_rgn_map(left.m_obj_rgn_map | right.m_obj_rgn_map);
    uf_var_map_t out_addr_var_map(left.m_addr_var_map | right.m_addr_var_map);

    // Merge allocator
    base_varname_allocator_t out_alloc(left.m_alloc, right.m_alloc);

    // Merge cache stuff
    boolean_value out_used(left.m_used | right.m_used);
    boolean_value out_dirty(left.m_dirty | right.m_dirty);
    c_base_addr_t out_base;

    // Perform renaming before join
    base_abstract_domain_t left_mem(left.m_mem);
    base_abstract_domain_t right_mem(right.m_mem);
    base_abstract_domain_t left_cache_lines(left.m_cache_lines);
    base_abstract_domain_t right_cache_lines(right.m_cache_lines);
    base_uf_domain_t left_regs(left.m_regs);
    base_uf_domain_t right_regs(right.m_regs);
    base_uf_domain_t left_addrs(left.m_addrs);
    base_uf_domain_t right_addrs(right.m_addrs);
    var_map_t out_var_map, out_cache_var_map;
    rev_var_map_t out_rev_var_map;
    // -- Compute common renamings
    base_variable_vector_t left_vars, right_vars, out_vars;
    // reserve vec space by upper bound to avoid reallocations
    size_t num_renamings = left.m_var_map.size();
    out_cache_var_map = std::move(left.m_cache_var_map);
    left_vars.reserve(num_renamings);
    right_vars.reserve(num_renamings);
    out_vars.reserve(num_renamings);

    for (auto &kv : left.m_var_map) {
      const variable_t &v = kv.first;
      auto it = right.m_var_map.find(v);
      if (it != right.m_var_map.end()) {
        // if we ignore unknown regions, types are the same by construction
        if (!crab_domain_params_man::get().region_skip_unknown_regions() &&
            !kv.second.same_type(it->second)) {
          crab::CrabStats::count(
              domain_name() +
              ".count.join_or_widen.skipped.inconsistent_dynamic_type");
          continue;
        }
        ghost_variables_t out_gvars(
            ghost_variables_t::create(out_alloc, kv.second));
        kv.second.add(left_vars);
        it->second.add(right_vars);
        out_gvars.add(out_vars);
        out_var_map.insert({v, out_gvars});
        out_gvars.update_rev_varmap(out_rev_var_map, v);

        // update cache used var map
        if (left.m_cache_var_map.find(v) != left.m_cache_var_map.end()) {
          if (right.m_cache_var_map.find(v) != right.m_cache_var_map.end()) {
            out_cache_var_map.insert({v, out_gvars});
          }
        }
      }
    }

    if (!(left.m_used == boolean_value::get_false()) &&
        !(right.m_used == boolean_value::get_false())) {
      assert(left.m_base);
      assert(right.m_base);
      if (!(*left.m_base).first.same_type((*right.m_base).first)) {
        // the type of ghost variable for base address should be the same
        CRAB_ERROR(domain_name(), "::do_join_or_widening left and right should "
                                  "contain the same m_base type");
      }
      // base address refers to the same object
      (*left.m_base).first.add(left_vars);
      (*right.m_base).first.add(right_vars);
      (*left.m_base).first.add(out_vars);
    }

    // JN: project might be necessary to avoid keeping variables that
    // exist in the base domain but they doen't exist on m_var_map.
    //
    // If such a variable exists only on either left_dom or right_dom
    // the join removes it.  However, if we have the same variable in
    // both left_dom and righ_dom then the join will preserve it.
    //
    left_mem.project(left_vars);
    left_mem.rename(left_vars, out_vars);
    right_mem.project(right_vars);
    right_mem.rename(right_vars, out_vars);

    left_cache_lines.project(left_vars);
    left_cache_lines.rename(left_vars, out_vars);
    right_cache_lines.project(right_vars);
    right_cache_lines.rename(right_vars, out_vars);

    left_regs.project(left_vars);
    left_regs.rename(left_vars, out_vars);
    right_regs.project(right_vars);
    right_regs.rename(right_vars, out_vars);

    left_addrs.project(left_vars);
    left_addrs.rename(left_vars, out_vars);
    right_addrs.project(right_vars);
    right_addrs.rename(right_vars, out_vars);

    // Final join or widening
    base_abstract_domain_t out_mem(base_dom_op(left_mem, right_mem));
    base_abstract_domain_t out_cache_lines(
        base_dom_op(left_cache_lines, right_cache_lines));
    base_uf_domain_t out_regs(base_uf_dom_op(left_regs, right_regs));
    base_uf_domain_t out_addrs(base_uf_dom_op(left_addrs, right_addrs));
    int out_uninterpreted_symbol_index = std::max(
        left.m_uninterpreted_symbol_index, right.m_uninterpreted_symbol_index);

    mru_region_domain_t res(
        std::move(out_alloc), std::move(out_var_map),
        std::move(out_rev_var_map), std::move(out_cache_var_map),
        std::move(out_mem), std::move(out_cache_lines), std::move(out_regs),
        std::move(out_addrs), std::move(out_base), std::move(out_used),
        std::move(out_dirty), std::move(out_rgn_counting_dom),
        std::move(out_rgn_init_dom), std::move(out_rgn_type_dom),
        std::move(out_rgn_dealloc_dom), std::move(out_tag_env),
        std::move(out_obj_rgn_map), std::move(out_addr_var_map),
        out_uninterpreted_symbol_index);
    res.m_is_bottom = res.m_mem.is_bottom();
    return res;
  }

  mru_region_domain_t do_meet_or_narrowing(
      const mru_region_domain_t &left_const,
      const mru_region_domain_t &right_const, const bool is_meet /*unused*/,
      base_dom_binop_t base_dom_op, base_uf_dom_binop_t base_uf_dom_op) const {

    mru_region_domain_t left = mru_region_domain_t(left_const);
    mru_region_domain_t right = mru_region_domain_t(right_const);
    // Check whether need to perform cache reduction
    if (!(left.m_used == boolean_value::get_false()) &&
        !(right.m_used == boolean_value::get_false())) { // if both dirty
      // left.m_base == right.m_base if refer to the same uniterpreted symbol
      if (!((*left.m_base).second == (*right.m_base).second)) {
        left.commit_cache();
        right.commit_cache();
      }
    } else if (!(left.m_used == boolean_value::get_false())) {
      // left: perform reduction before join
      left.commit_cache();
    } else if (!(right.m_used == boolean_value::get_false())) {
      // right: perform reduction before join
      right.commit_cache();
    }

    // At this point, the left and right caches are
    // (1) empty, (2) dirty but refer to the same object

    // The following domains do not require common renaming
    // (i.e. no ghost variables). The domains are finite.
    rgn_counting_env_t out_rgn_counting_dom(left.m_rgn_counting_dom &
                                            right.m_rgn_counting_dom);
    rgn_bool_env_t out_rgn_init_dom(left.m_rgn_init_dom & right.m_rgn_init_dom);
    rgn_type_env_t out_rgn_type_dom(left.m_rgn_type_dom & right.m_rgn_type_dom);
    rgn_dealloc_t out_rgn_dealloc_dom(left.m_rgn_dealloc_dom &
                                      right.m_rgn_dealloc_dom);
    tag_env_t out_tag_env(left.m_tag_env & right.m_tag_env);
    rgn_equiv_t out_obj_rgn_map(left.m_obj_rgn_map & right.m_obj_rgn_map);
    uf_var_map_t out_addr_var_map(left.m_addr_var_map & right.m_addr_var_map);

    // Merge allocator
    base_varname_allocator_t out_alloc(left.m_alloc, right.m_alloc);

    // Merge cache stuff
    boolean_value out_used(left.m_used & right.m_used);
    boolean_value out_dirty(left.m_dirty & right.m_dirty);
    c_base_addr_t out_base;

    // Perform renaming before meet
    base_abstract_domain_t left_mem(left.m_mem);
    base_abstract_domain_t right_mem(right.m_mem);
    base_abstract_domain_t left_cache_lines(left.m_cache_lines);
    base_abstract_domain_t right_cache_lines(right.m_cache_lines);
    base_uf_domain_t left_regs(left.m_regs);
    base_uf_domain_t right_regs(right.m_regs);
    base_uf_domain_t left_addrs(left.m_addrs);
    base_uf_domain_t right_addrs(right.m_addrs);
    var_map_t out_var_map, out_cache_var_map;
    rev_var_map_t out_rev_var_map;
    // -- Compute common renamings
    // out vars are used for output domain
    base_variable_vector_t left_vars, right_vars, out_vars;
    // Only left (right) track variables used in left (right) dom only (not in
    // common)
    base_variable_vector_t only_left_vars, only_left_out_vars;
    base_variable_vector_t only_right_vars, only_right_out_vars;
    // reserve vec space by upper bound to avoid reallocations
    size_t left_renamings = left.m_var_map.size();
    size_t right_renamings = right.m_var_map.size();
    size_t num_renamings = left_renamings + right_renamings;
    left_vars.reserve(num_renamings);
    right_vars.reserve(num_renamings);
    out_vars.reserve(num_renamings);
    only_left_vars.reserve(left_renamings);
    only_left_out_vars.reserve(left_renamings);
    only_right_vars.reserve(right_renamings);
    only_right_out_vars.reserve(right_renamings);

    for (auto &kv : left.m_var_map) {
      const variable_t &v = kv.first;
      auto it = right.m_var_map.find(v);
      if (it != right.m_var_map.end()) {
        // if we ignore unknown regions, types are the same by construction
        if (!crab_domain_params_man::get().region_skip_unknown_regions() &&
            !kv.second.same_type(it->second)) {
          crab::CrabStats::count(
              domain_name() +
              ".count.meet_or_narrowing.skipped.inconsistent_dynamic_type");
          continue;
        }
        // common renaming
        ghost_variables_t out_gvars(
            ghost_variables_t::create(out_alloc, kv.second));
        kv.second.add(left_vars);
        it->second.add(right_vars);
        out_gvars.add(out_vars);
        out_var_map.insert({v, out_gvars});
        out_gvars.update_rev_varmap(out_rev_var_map, v);
        // update cache used var map
        if (left.m_cache_var_map.find(v) != left.m_cache_var_map.end()) {
          if (right.m_cache_var_map.find(v) != right.m_cache_var_map.end()) {
            out_cache_var_map.insert({v, out_gvars});
          }
        }
      } else {
        // only on left
        ghost_variables_t out_gvars(
            ghost_variables_t::create(out_alloc, kv.second));
        kv.second.add(only_left_vars);
        out_gvars.add(only_left_out_vars);
        out_var_map.insert({kv.first, out_gvars});
        out_gvars.update_rev_varmap(out_rev_var_map, kv.first);
        // update cache used var map
        if (left.m_cache_var_map.find(v) != left.m_cache_var_map.end()) {
          out_cache_var_map.insert({v, out_gvars});
        }
      }
    }

    // add missing maps from right operand
    for (auto &kv : right.m_var_map) {
      const variable_t &v = kv.first;
      auto it = left.m_var_map.find(v);
      if (it == left.m_var_map.end()) {
        // only on right
        ghost_variables_t out_gvars(
            ghost_variables_t::create(out_alloc, kv.second));
        kv.second.add(only_right_vars);
        out_gvars.add(only_right_out_vars);
        out_var_map.insert({kv.first, out_gvars});
        out_gvars.update_rev_varmap(out_rev_var_map, kv.first);
        // update cache used var map
        if (right.m_cache_var_map.find(v) != right.m_cache_var_map.end()) {
          out_cache_var_map.insert({v, out_gvars});
        }
      }
    }

    if (!(left.m_used == boolean_value::get_false()) &&
        !(right.m_used == boolean_value::get_false())) {
      assert(left.m_base);
      assert(right.m_base);
      if (!(*left.m_base).first.same_type((*right.m_base).first)) {
        // the type of ghost variable for base address should be the same
        CRAB_ERROR(domain_name(), "::do_join_or_widening left and right should "
                                  "contain the same m_base type");
      }
      // base address refers to the same object
      (*left.m_base).first.add(left_vars);
      (*right.m_base).first.add(right_vars);
      (*left.m_base).first.add(out_vars);
    }

    // append common into only left
    only_left_vars.insert(only_left_vars.end(), left_vars.begin(),
                          left_vars.end());
    // append common into only left out
    only_left_out_vars.insert(only_left_out_vars.end(), out_vars.begin(),
                              out_vars.end());
    // need to project before renaming
    left_mem.project(only_left_vars);
    left_mem.rename(only_left_vars, only_left_out_vars);
    left_cache_lines.project(only_left_vars);
    left_cache_lines.rename(only_left_vars, only_left_out_vars);
    left_regs.project(only_left_vars);
    left_regs.rename(only_left_vars, only_left_out_vars);
    left_addrs.project(only_left_vars);
    left_addrs.rename(only_left_vars, only_left_out_vars);

    // append common into only right
    only_right_vars.insert(only_right_vars.end(), right_vars.begin(),
                           right_vars.end());
    // append common into only right out
    only_right_out_vars.insert(only_right_out_vars.end(), out_vars.begin(),
                               out_vars.end());
    // need to project before renaming
    right_mem.project(only_left_vars);
    right_mem.rename(only_left_vars, only_left_out_vars);
    right_cache_lines.project(only_right_vars);
    right_cache_lines.rename(only_right_vars, only_right_out_vars);
    right_regs.project(only_right_vars);
    right_regs.rename(only_right_vars, only_right_out_vars);
    right_addrs.project(only_left_vars);
    right_addrs.rename(only_left_vars, only_left_out_vars);

    // Final meet or narrowing
    base_abstract_domain_t out_mem(base_dom_op(left_mem, right_mem));
    base_abstract_domain_t out_cache_lines(
        base_dom_op(left_cache_lines, right_cache_lines));
    base_uf_domain_t out_regs(base_uf_dom_op(left_regs, right_regs));
    base_uf_domain_t out_addrs(base_uf_dom_op(left_addrs, right_addrs));

    int out_uninterpreted_symbol_index = std::max(
        left.m_uninterpreted_symbol_index, right.m_uninterpreted_symbol_index);

    mru_region_domain_t res(
        std::move(out_alloc), std::move(out_var_map),
        std::move(out_rev_var_map), std::move(out_cache_var_map),
        std::move(out_mem), std::move(out_cache_lines), std::move(out_regs),
        std::move(out_addrs), std::move(out_base), std::move(out_used),
        std::move(out_dirty), std::move(out_rgn_counting_dom),
        std::move(out_rgn_init_dom), std::move(out_rgn_type_dom),
        std::move(out_rgn_dealloc_dom), std::move(out_tag_env),
        std::move(out_obj_rgn_map), std::move(out_addr_var_map),
        out_uninterpreted_symbol_index);
    res.m_is_bottom = res.m_mem.is_bottom();
    return res;
  }

  const ghost_variables_t &get_or_insert_gvars(const variable_t &v) {
    return get_or_insert_gvars(v, m_var_map, m_rev_var_map, m_alloc);
  }

  const ghost_variables_t &
  get_or_insert_gvars(const variable_t &v, var_map_t &varmap,
                      rev_var_map_t &rev_varmap,
                      base_varname_allocator_t &alloc) const {
    auto it = varmap.find(v);
    if (it != varmap.end()) {
      return it->second;
    } else {
      variable_type vty = get_dynamic_type_or_fail(v);
      auto gvars = ghost_variables_t::create(alloc, vty, __LINE__);
      auto res = varmap.insert({v, gvars});
      if (res.second) {
        gvars.update_rev_varmap(rev_varmap, v);
      }
      return res.first->second;
    }
  }

  ghost_variables_t get_gvars_or_fail(const variable_t &v) const {
    auto it = m_var_map.find(v);
    if (it != m_var_map.end()) {
      return it->second;
    } else {
      CRAB_ERROR("get_gvars_or_fail failed");
    }
  }

  boost::optional<ghost_variables_t> get_gvars(const variable_t &v) const {
    auto it = m_var_map.find(v);
    if (it != m_var_map.end()) {
      return it->second;
    } else {
      return boost::none;
    }
  }

  base_variable_or_constant_t
  rename_variable_or_constant(const variable_or_constant_t &v) {
    if (v.is_constant()) {
      return base_variable_or_constant_t(v.get_constant(), v.get_type());
    } else {
      base_variable_t bv = get_or_insert_gvars(v.get_variable()).get_var();
      return base_variable_or_constant_t(bv);
    }
  }

  // Rename linear expression
  base_linear_expression_t rename_linear_expr(const linear_expression_t &e,
                                              bool is_cache_rename = false) {
    base_linear_expression_t out(e.constant());
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const variable_t &v = (*it).second;
      const number_t &coef = (*it).first;
      const ghost_variables_t &new_gvars = get_or_insert_gvars(v);
      base_variable_t nv = new_gvars.get_var();
      m_cache_var_map.insert({v, new_gvars});
      out = out + (coef * nv);
    }
    return out;
  }

  // Rename linear constraint
  base_linear_constraint_t rename_linear_cst(const linear_constraint_t &cst,
                                             bool is_cache_rename = false) {
    if (cst.is_inequality() || cst.is_strict_inequality()) {
      return base_linear_constraint_t(
          rename_linear_expr(cst.expression(), is_cache_rename),
          (typename base_linear_constraint_t::kind_t)cst.kind(),
          cst.is_signed());

    } else {
      return base_linear_constraint_t(
          rename_linear_expr(cst.expression(), is_cache_rename),
          (typename base_linear_constraint_t::kind_t)cst.kind());
    }
  }

  // used only for to_linear_constraint_system()
  boost::optional<variable_t> rev_rename_var(const base_variable_t &v,
                                             bool ignore_references) const {
    auto it = m_rev_var_map.find(v);
    if (it != m_rev_var_map.end()) {
      variable_t v = it->second.first;
      ghost_var_id ghost_id = it->second.second;
      if (!ignore_references || !v.get_type().is_reference()) {
        if (ghost_id == 1) { // we only get the first ghost variable
                             // which represents the address of the
                             // reference. Skipping other ghost
                             // variables such as offsets and sizes.
          return v;
        }
      }
    }
    return boost::none;
  }

  // used only for to_linear_constraint_system()
  boost::optional<linear_expression_t>
  rev_rename_linear_expr(const base_linear_expression_t &e,
                         bool ignore_references) const {
    linear_expression_t out(e.constant());
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const base_variable_t &v = (*it).second;
      const number_t &coef = (*it).first;
      if (boost::optional<variable_t> v_opt =
              rev_rename_var(v, ignore_references)) {
        out = out + (coef * (*v_opt));
      } else {
        return boost::optional<linear_expression_t>();
      }
    }
    return out;
  }

  // used only for to_linear_constraint_system()
  boost::optional<linear_constraint_t>
  rev_rename_linear_cst(const base_linear_constraint_t &cst,
                        bool ignore_references) const {
    if (boost::optional<linear_expression_t> e =
            rev_rename_linear_expr(cst.expression(), ignore_references)) {

      if (cst.is_inequality() || cst.is_strict_inequality()) {
        return linear_constraint_t(
            *e, (typename linear_constraint_t::kind_t)cst.kind(),
            cst.is_signed());

      } else {
        return linear_constraint_t(
            *e, (typename linear_constraint_t::kind_t)cst.kind());
      }
    } else {
      return boost::optional<linear_constraint_t>();
    }
  }

  // Rename linear constraints
  base_linear_constraint_system_t
  rename_linear_cst_sys(const linear_constraint_system_t &csts,
                        bool is_cache_rename = false) {
    base_linear_constraint_system_t out;
    for (auto const &cst : csts) {
      out += rename_linear_cst(cst, is_cache_rename);
    }
    return out;
  }

  base_linear_constraint_t
  convert_ref_cst_to_linear_cst(const reference_constraint_t &ref_cst) {
    if (ref_cst.is_tautology()) {
      return base_linear_constraint_t::get_true();
    } else if (ref_cst.is_contradiction()) {
      return base_linear_constraint_t::get_false();
    } else {
      if (ref_cst.is_unary()) {
        assert(ref_cst.lhs().get_type().is_reference());
        base_variable_t x = get_or_insert_gvars(ref_cst.lhs()).get_var();
        if (ref_cst.is_equality()) {
          return base_linear_constraint_t(x == number_t(0));
        } else if (ref_cst.is_disequality()) {
          return base_linear_constraint_t(x != number_t(0));
        } else if (ref_cst.is_less_or_equal_than()) {
          return base_linear_constraint_t(x <= number_t(0));
        } else if (ref_cst.is_less_than()) {
          return base_linear_constraint_t(x < number_t(0));
        } else if (ref_cst.is_greater_or_equal_than()) {
          return base_linear_constraint_t(x >= number_t(0));
        } else if (ref_cst.is_greater_than()) {
          return base_linear_constraint_t(x > number_t(0));
        }
      } else {
        assert(ref_cst.lhs().get_type().is_reference());
        assert(ref_cst.rhs().get_type().is_reference());
        base_variable_t x = get_or_insert_gvars(ref_cst.lhs()).get_var();
        base_variable_t y = get_or_insert_gvars(ref_cst.rhs()).get_var();
        number_t offset = ref_cst.offset();
        if (ref_cst.is_equality()) {
          return base_linear_constraint_t(x == y + offset);
        } else if (ref_cst.is_disequality()) {
          return base_linear_constraint_t(x != y + offset);
        } else if (ref_cst.is_less_or_equal_than()) {
          return base_linear_constraint_t(x <= y + offset);
        } else if (ref_cst.is_less_than()) {
          return base_linear_constraint_t(x < y + offset);
        } else if (ref_cst.is_greater_or_equal_than()) {
          return base_linear_constraint_t(x >= y + offset);
        } else if (ref_cst.is_greater_than()) {
          return base_linear_constraint_t(x > y + offset);
        }
      }
    }
    CRAB_ERROR("unexpected reference constraint");
  }

  template <class RangeVars>
  void merge_tags(const variable_t &x, RangeVars vars) {
    tag_set tags = tag_set::bottom();
    for (auto const &v : vars) {
      tags = tags | m_tag_env[v];
    }
    m_tag_env.set(x, tags);
  }

  bool is_tracked_region(const variable_t &v) const {
    auto ty = v.get_type();
    if (!ty.is_region()) {
      return false;
    }
    if (!ty.is_unknown_region()) {
      return true;
    }
    return (crab_domain_params_man::get().region_skip_unknown_regions()
                ? false
                : has_dynamic_type(v));
  }

  static bool is_tracked_unknown_region(const variable_t &v) {
    return (!crab_domain_params_man::get().region_skip_unknown_regions() &&
            v.get_type().is_unknown_region());
  }

  // if the variable has a static type (different from unknown region)
  // then returns true.  If the static type is an unknown region then
  // the actual variable's type must be inferred dynamically (as part
  // of the abstract domain). In this case, it will return true if the
  // type inference done by the abstract domain succeeds (type is
  // inferred different from unknown region).
  bool has_dynamic_type(const variable_t &v) const {
    if (!v.get_type().is_unknown_region()) {
      return true;
    }

    if (crab_domain_params_man::get().region_skip_unknown_regions()) {
      return false;
    } else {
      type_value dyn_ty = m_rgn_type_dom[v];
      if (dyn_ty.is_top() || dyn_ty.is_bottom()) {
        return false;
      }
      return (!dyn_ty.get().is_unknown_region());
    }
  }

  // if has_dynamic_type(v) returns true then it returns the type,
  // otherwise it will fail.
  variable_type get_dynamic_type_or_fail(const variable_t &v) const {
    assert(has_dynamic_type(v));

    if (!v.get_type().is_unknown_region()) {
      return v.get_type();
    }

    type_value dyn_ty = m_rgn_type_dom[v];
    if (dyn_ty.is_top() || dyn_ty.is_bottom()) {
      CRAB_ERROR("get_dynamic_type_or_fail cannot be called on top or bottom");
    }
    if (dyn_ty.get().is_unknown_region()) {
      CRAB_ERROR("get_dynamic_type_or_fail should not return unknown region");
    }
    return dyn_ty.get();
  }

  static void ERROR_IF_NOT_REGION(const variable_t &v, unsigned line) {
    if (!v.get_type().is_region()) {
      CRAB_ERROR(v, ":", v.get_type(), " is not a region at line ", line);
    }
  }

  static void ERROR_IF_ARRAY_REGION(const variable_t &v, unsigned line) {
    if (v.get_type().is_array_region()) {
      CRAB_ERROR(v, ":", v.get_type(), " cannot contain an array at line ",
                 line);
    }
  }

  static void ERROR_IF_NOT_REF(const variable_t &v, unsigned line) {
    if (!v.get_type().is_reference()) {
      CRAB_ERROR(v, ":", v.get_type(), " is not a reference at line ", line);
    }
  }

  static void ERROR_IF_NOT_INT(const variable_t &v, unsigned line) {
    if (!v.get_type().is_integer()) {
      CRAB_ERROR(v, ":", v.get_type(), " is not an integer at line ", line);
    }
  }

  /**------------------ End helper method definitions ------------------**/

  /**------------- Begin MRU helper method definitions -------------------**/
  // perform a similiar work for cache.lines.constraint(cache.regs)
  // This one requires to be renaming because mem base dom used different base
  // vars
  base_abstract_domain_t cache_lines_regs_reduce() {
    base_abstract_domain_t res(m_cache_lines);
    auto csts = m_regs.to_linear_constraint_system();
    res += csts;
    return res;
  }

  // get all used Crab IR vars by region dom
  void vars(variable_vector_t &used_vars) const {
    used_vars.reserve(m_var_map.size());
    for (auto &kv : m_var_map) {
      used_vars.push_back(kv.first);
    }
  }

  void region_vars_in_mem(variable_vector_t &used_vars) const {
    used_vars.reserve(m_var_map.size());
    for (auto &kv : m_var_map) {
      if (kv.first.get_type().is_region()) {
        used_vars.push_back(kv.first);
      }
    }
  }

  void region_vars_in_cache(variable_vector_t &used_vars) const {
    used_vars.reserve(m_cache_var_map.size());
    for (auto &kv : m_cache_var_map) {
      if (kv.first.get_type().is_region()) {
        used_vars.push_back(kv.first);
      }
    }
  }

  void vars_in_cache(variable_vector_t &used_vars) const {
    used_vars.reserve(m_cache_var_map.size());
    for (auto &kv : m_cache_var_map) {
      used_vars.push_back(kv.first);
    }
  }

  variable_vector_t vec_minus(variable_vector_t left, variable_vector_t right) {
    using variable_map_t = ikos::patricia_tree<variable_t, variable_t>;
    variable_map_t vars;
    for (auto it = left.begin(); it != left.end(); ++it) {
      vars.insert(*it, *it);
    }
    for (auto it = right.begin(); it != right.end(); ++it) {
      if (vars.find(*it) != nullptr)
        vars.remove(*it);
    }
    variable_vector_t res;
    for (auto it = vars.begin(); it != vars.end(); ++it) {
      res.push_back(it->first);
    }
    return res;
  }

  // NOTE: For debug use
  void print_base_vec(const base_variable_vector_t &vec) {
    for (auto elem : vec)
      crab::outs() << elem << " ";
    crab::outs() << "\n";
  }

  // NOTE: For debug use
  void print_var_map() {
    crab::outs() << "Cache varmap\t\n";
    for (auto &kv : m_var_map)
      crab::outs() << kv.first << " ";
    crab::outs() << "\n";
  }

  // NOTE: For debug use
  void print_vec(const variable_vector_t &vec) {
    for (auto elem : vec)
      crab::outs() << elem << " ";
    crab::outs() << "\n";
  }

  void invalidate_cache_if_miss(const variable_t &rgn, const variable_t &base) {
    crab::CrabStats::count(domain_name() + ".count.invalidate_cache_if_miss");
    crab::ScopedCrabStats __st__(domain_name() + ".invalidate_cache_if_miss");

    if (m_base == boost::none) { // cache is empty
      update_cache(rgn, base);
      return;
    }

    base_variable_t c_base = (*m_base).first.get_var();
    base_variable_t base_in_base = get_or_insert_gvars(base).get_var();
    base_linear_constraint_t addrs_cons =
        base_linear_constraint_t(base_linear_expression_t(c_base) !=
                                 base_linear_expression_t(base_in_base));

    base_linear_constraint_t tmp =
        base_linear_constraint_t(base_linear_expression_t(c_base) ==
                                 base_linear_expression_t(base_in_base));
    // crab::outs() << "check: "<< "C_base == " << base << ": " <<
    //         c_base << "==" << base_in_base << " Result: " <<
    //         checker_domain_traits<base_uf_domain_t>::entail(
    //         m_addrs, tmp) << "\n";
    // crab::outs() << m_addrs << "\n";
    uninterpreted_symbol_t cbase_symb = (*m_base).second;
    uninterpreted_symbol_t base_symb = get_or_insert_uninterpreted_symbol(base);
    if (m_used == boolean_value::get_true() && !(cbase_symb == base_symb)) { // C_used == true && C_base != base
      CRAB_LOG("mru-region-cache", crab::outs() << "cache missed\n";);
      // Cache missed:
      // Step1: commit cache if cache is dirty
      // if cache is used, C_dirty can be true or top
      if (!(m_dirty == boolean_value::get_false())) { // C_dirty != false
        commit_cache();
      } else {
        // if cache is not dirty, clear the cache
        CRAB_LOG("mru-region-cache", crab::outs() << "clear cache\n";);
        make_fresh_cache();
      }
      // Step2: update cache for new MRU object
      update_cache(rgn, base);
    }
  }

  // fresh cache as: { cache lines: bot, regs: top, addrs: addrs -= C_base }
  // for addrs, we keep all equalities of base address for each reference
  // also need to assign cache used flag C_used as false
  void make_fresh_cache() {
    m_cache_lines.set_to_top();
    m_regs.set_to_top();
    if (m_base) {
      m_addrs.forget({(*m_base).first.get_var()});
    }
    m_used = boolean_value::get_false();  // C_used = false
    m_dirty = boolean_value::get_false(); // C_dirty = false
  }

  // commit cache contents into memory
  // require to fresh cache after commit.
  void commit_cache() {
    fold();
    make_fresh_cache();
  }

  void update_cache(const variable_t &rgn, const variable_t &base) {
    // clear cache var map
    m_cache_var_map.clear();
    expand(rgn); // perform cache_line expand

    if (!m_addr_var_map.contains(base)) {
      CRAB_ERROR(domain_name(),
                 "::update_cache should found the uninterpreted symbol for "
                 "the ref's base address");
    }

    if (m_base == boost::none) { // cache is empty
      // need to construct a ghost variable for address domain
      ghost_variables_t ghost_base =
          ghost_variables_t::create(m_alloc, base.get_type(), __LINE__);
      m_base = std::make_pair(ghost_base, m_addr_var_map.get(base));
    } else {
      m_base = std::make_pair((*m_base).first, m_addr_var_map.get(base));
    }

    base_variable_t c_base = (*m_base).first.get_var();
    base_variable_t base_in_base = get_or_insert_gvars(base).get_var();
    m_addrs.assign(c_base, base_in_base); // C_base == base
    m_used = boolean_value::get_true();   // C_used = true
    m_dirty = boolean_value::get_false(); // C_dirty = false
  }

  /**------------- End MRU helper method definitions -------------------**/

public:
  /**------------------ Begin domain API definitions ------------------**/
  mru_region_domain_t make_top() const override {
    return mru_region_domain_t(true);
  }

  mru_region_domain_t make_bottom() const override {
    return mru_region_domain_t(false);
  }

  void set_to_top() override {
    mru_region_domain_t abs(true);
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    mru_region_domain_t abs(false);
    std::swap(*this, abs);
  }

  mru_region_domain(bool is_top = true)
      : m_is_bottom(!is_top), m_base(boost::none) {
    make_fresh_cache();
  }

  mru_region_domain(const mru_region_domain_t &o)
      : m_is_bottom(o.m_is_bottom), m_alloc(o.m_alloc), m_var_map(o.m_var_map),
        m_rev_var_map(o.m_rev_var_map), m_cache_var_map(o.m_cache_var_map),
        m_mem(o.m_mem), m_cache_lines(o.m_cache_lines), m_regs(o.m_regs),
        m_addrs(o.m_addrs), m_base(o.m_base), m_used(o.m_used),
        m_dirty(o.m_dirty), m_rgn_counting_dom(o.m_rgn_counting_dom),
        m_rgn_init_dom(o.m_rgn_init_dom), m_rgn_type_dom(o.m_rgn_type_dom),
        m_rgn_dealloc_dom(o.m_rgn_dealloc_dom), m_tag_env(o.m_tag_env),
        m_obj_rgn_map(o.m_obj_rgn_map), m_addr_var_map(o.m_addr_var_map),
        m_uninterpreted_symbol_index(o.m_uninterpreted_symbol_index) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }

  mru_region_domain(mru_region_domain_t &&o)
      : m_is_bottom(o.m_is_bottom), m_alloc(std::move(o.m_alloc)),
        m_var_map(std::move(o.m_var_map)),
        m_rev_var_map(std::move(o.m_rev_var_map)),
        m_cache_var_map(std::move(o.m_cache_var_map)),
        m_mem(std::move(o.m_mem)), m_cache_lines(std::move(o.m_cache_lines)),
        m_regs(std::move(o.m_regs)), m_addrs(std::move(o.m_addrs)),
        m_base(std::move(o.m_base)), m_used(std::move(o.m_used)),
        m_dirty(std::move(o.m_dirty)),
        m_rgn_counting_dom(std::move(o.m_rgn_counting_dom)),
        m_rgn_init_dom(std::move(o.m_rgn_init_dom)),
        m_rgn_type_dom(std::move(o.m_rgn_type_dom)),
        m_rgn_dealloc_dom(std::move(o.m_rgn_dealloc_dom)),
        m_tag_env(std::move(o.m_tag_env)),
        m_obj_rgn_map(std::move(o.m_obj_rgn_map)),
        m_addr_var_map(std::move(o.m_addr_var_map)),
        m_uninterpreted_symbol_index(o.m_uninterpreted_symbol_index) {}

  mru_region_domain_t &operator=(const mru_region_domain_t &o) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &o) {
      m_is_bottom = o.m_is_bottom;
      m_alloc = o.m_alloc;
      m_var_map = o.m_var_map;
      m_rev_var_map = o.m_rev_var_map;
      m_cache_var_map = o.m_cache_var_map;
      m_mem = o.m_mem;
      m_cache_lines = o.m_cache_lines;
      m_regs = o.m_regs;
      m_addrs = o.m_addrs;
      m_base = o.m_base;
      m_used = o.m_used;
      m_dirty = o.m_dirty;
      m_rgn_counting_dom = o.m_rgn_counting_dom;
      m_rgn_init_dom = o.m_rgn_init_dom;
      m_rgn_type_dom = o.m_rgn_type_dom;
      m_rgn_dealloc_dom = o.m_rgn_dealloc_dom;
      m_tag_env = o.m_tag_env;
      m_obj_rgn_map = o.m_obj_rgn_map;
      m_addr_var_map = o.m_addr_var_map;
      m_uninterpreted_symbol_index = o.m_uninterpreted_symbol_index;
    }
    return *this;
  }

  mru_region_domain_t &operator=(mru_region_domain_t &&o) {
    if (this != &o) {
      m_is_bottom = o.m_is_bottom;
      m_alloc = std::move(o.m_alloc);
      m_var_map = std::move(o.m_var_map);
      m_rev_var_map = std::move(o.m_rev_var_map);
      m_cache_var_map = std::move(o.m_cache_var_map);
      m_mem = std::move(o.m_mem);
      m_cache_lines = std::move(o.m_cache_lines);
      m_regs = std::move(o.m_regs);
      m_addrs = std::move(o.m_addrs);
      m_base = std::move(o.m_base);
      m_used = std::move(o.m_used);
      m_dirty = std::move(o.m_dirty);
      m_rgn_counting_dom = std::move(o.m_rgn_counting_dom);
      m_rgn_init_dom = std::move(o.m_rgn_init_dom);
      m_rgn_type_dom = std::move(o.m_rgn_type_dom);
      m_rgn_dealloc_dom = std::move(o.m_rgn_dealloc_dom);
      m_tag_env = std::move(o.m_tag_env);
      m_obj_rgn_map = std::move(o.m_obj_rgn_map);
      m_addr_var_map = std::move(o.m_addr_var_map);
      m_uninterpreted_symbol_index = o.m_uninterpreted_symbol_index;
    }
    return *this;
  }

  bool is_bottom() const override { return m_is_bottom; }

  bool is_top() const override {
    bool res = !is_bottom() && m_mem.is_top() && m_cache_lines.is_top() &&
               m_regs.is_top() && m_addrs.is_top() &&
               m_rgn_counting_dom.is_top();
    if (crab_domain_params_man::get().region_deallocation()) {
      res = res && m_rgn_dealloc_dom.is_top();
    }
    if (crab_domain_params_man::get().region_tag_analysis()) {
      res = res && m_tag_env.is_top();
    }
    return res;
  }

  bool operator<=(const mru_region_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");

    CRAB_LOG("mru-region-leq", crab::outs() << "Inclusion test:\n\t" << *this
                                            << "\n\t" << o << "\n";);
    if (is_bottom() || o.is_top()) { // this indeed <= o
      CRAB_LOG("mru-region-leq", crab::outs() << "Result1=1\n";);
      return true;
    } else if (is_top() || o.is_bottom()) { // this indeed > o
      CRAB_LOG("mru-region-leq", crab::outs() << "Result2=0\n";);
      return false;
    }

    if (!(m_rgn_counting_dom <= o.m_rgn_counting_dom)) {
      CRAB_LOG("mru-region-leq", crab::outs() << "Result3=0\n";);
      return false;
    }

    if (!(m_rgn_init_dom <= o.m_rgn_init_dom)) {
      CRAB_LOG("mru-region-leq", crab::outs() << "Result4=0\n";);
      return false;
    }

    if (!(m_rgn_type_dom <= o.m_rgn_type_dom)) {
      CRAB_LOG("mru-region-leq", crab::outs() << "Result5=0\n";);
      return false;
    }

    if (crab_domain_params_man::get().region_deallocation()) {
      if (!(m_rgn_dealloc_dom <= o.m_rgn_dealloc_dom)) {
        CRAB_LOG("mru-region-leq", crab::outs() << "Result6=0\n";);
        return false;
      }
    }

    if (crab_domain_params_man::get().region_tag_analysis()) {
      if (!(m_tag_env <= o.m_tag_env)) {
        CRAB_LOG("mru-region-leq", crab::outs() << "Result7=0\n";);
        return false;
      }
    }

    if (!(m_obj_rgn_map <= o.m_obj_rgn_map)) {
      CRAB_LOG("mru-region-leq", crab::outs() << "Result8=0\n";);
      return false;
    }

    base_varname_allocator_t out_alloc(m_alloc, o.m_alloc);

    base_abstract_domain_t left_mem(m_mem);
    base_abstract_domain_t right_mem(o.m_mem);
    base_abstract_domain_t left_cache_lines(m_cache_lines);
    base_abstract_domain_t right_cache_lines(o.m_cache_lines);
    base_uf_domain_t left_regs(m_regs);
    base_uf_domain_t right_regs(o.m_regs);
    base_uf_domain_t left_addrs(m_addrs);
    base_uf_domain_t right_addrs(o.m_addrs);

    // perform the common renaming
    base_variable_vector_t left_vars, right_vars, out_vars;
    // reserve vec space by upper bound to avoid reallocations
    size_t num_renamings = m_var_map.size();
    left_vars.reserve(num_renamings);
    right_vars.reserve(num_renamings);
    out_vars.reserve(num_renamings);

    for (auto &kv : m_var_map) {
      const variable_t &v = kv.first;
      auto it = o.m_var_map.find(v);
      if (it != o.m_var_map.end()) {
        // if we ignore unknown regions, types are the same by construction
        if (!crab_domain_params_man::get().region_skip_unknown_regions() &&
            !kv.second.same_type(it->second)) {
          CRAB_ERROR(domain_name(), "::operator<= ", " dynamic type of ", v,
                     " must be the same in both left and right operands");
        }
        ghost_variables_t out_gvars(
            ghost_variables_t::create(out_alloc, kv.second));
        kv.second.add(left_vars);
        it->second.add(right_vars);
        out_gvars.add(out_vars);
      }
    }

    left_mem.project(left_vars);
    left_mem.rename(left_vars, out_vars);
    right_mem.project(right_vars);
    right_mem.rename(right_vars, out_vars);

    left_cache_lines.project(left_vars);
    left_cache_lines.rename(left_vars, out_vars);
    right_cache_lines.project(right_vars);
    right_cache_lines.rename(right_vars, out_vars);

    left_regs.project(left_vars);
    left_regs.rename(left_vars, out_vars);
    right_regs.project(right_vars);
    right_regs.rename(right_vars, out_vars);

    left_addrs.project(left_vars);
    left_addrs.rename(left_vars, out_vars);
    right_addrs.project(right_vars);
    right_addrs.rename(right_vars, out_vars);

    CRAB_LOG("mru-region-leq", crab::outs()
                                   << "Inclusion test (after renaming):\n\t"
                                   << "cache line:" << left_cache_lines
                                   << "\n\t" << right_cache_lines << "\n"
                                   << "reg" << left_regs << "\n\t" << right_regs
                                   << "\n";);

    bool res = left_mem <= right_mem && left_cache_lines <= right_cache_lines &&
               left_regs <= right_regs && left_addrs <= right_addrs;
    CRAB_LOG("mru-region-leq", crab::outs() << "Result9=" << res << "\n";);
    return res;
  }

  void operator|=(const mru_region_domain_t &o) override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    if (is_bottom()) { // bot | o = o
      *this = o;
      return;
    } else if (o.is_bottom()) { // this | bot = this
      return;
    } else if (is_top() || o.is_top()) { // top | o or this | top, return top
      set_to_top();
      return;
    }

    CRAB_LOG("mru-region", crab::outs()
                               << "Join " << *this << " and " << o << "=");

    variable_vector_t left_regions, right_regions;
    std::vector<ghost_variables_t> regions_left_base_vars,
        regions_right_base_vars;

    do_join(o);
  }

  mru_region_domain_t operator|(const mru_region_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    if (is_bottom()) { // bot | o = o
      return o;
    } else if (o.is_bottom()) { // this | bot = this
      return *this;
    } else if (is_top() || o.is_top()) { // top | o or this | top, return top
      mru_region_domain_t abs;
      abs.set_to_top();
      return abs;
    }

    CRAB_LOG("mru-region", crab::outs()
                               << "Join " << *this << " and " << o << "=");

    auto base_dom_op = [](const base_abstract_domain_t &v1,
                          const base_abstract_domain_t &v2) { return v1 | v2; };
    auto base_uf_dom_op = [](const base_uf_domain_t &v1,
                             const base_uf_domain_t &v2) { return v1 | v2; };

    variable_vector_t left_regions, right_regions;
    std::vector<ghost_variables_t> regions_left_base_vars,
        regions_right_base_vars;

    mru_region_domain_t res(std::move(do_join_or_widening(
        *this, o, true /*is join*/, base_dom_op, base_uf_dom_op)));

    CRAB_LOG("region", crab::outs() << res << "\n");
    return res;
  }

  mru_region_domain_t operator&(const mru_region_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    if (is_bottom() || o.is_top()) { // bot & o or this & top, return this
      return *this;
    } else if (o.is_bottom() || is_top()) { // this & bot or top & o, return o
      return o;
    }

    CRAB_LOG("mru-region", crab::outs()
                               << "Meet " << *this << " and " << o << "=");

    auto base_dom_op = [](const base_abstract_domain_t &v1,
                          const base_abstract_domain_t &v2) { return v1 & v2; };
    auto base_uf_dom_op = [](const base_uf_domain_t &v1,
                             const base_uf_domain_t &v2) { return v1 & v2; };

    mru_region_domain_t res(std::move(do_meet_or_narrowing(
        *this, o, true /*is meet*/, base_dom_op, base_uf_dom_op)));

    CRAB_LOG("mru-region", crab::outs() << res << "\n");
    return res;
  }

  mru_region_domain_t operator||(const mru_region_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");

    // Trivial cases first: we don't cover cases where one operand is
    // top because is_top() calls the base domain which we don't know
    // whether it will perform some normalization or not.
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    }

    CRAB_LOG("mru-region", crab::outs()
                               << "Widening " << *this << " and " << o << "=");

    auto base_dom_op = [](const base_abstract_domain_t &v1,
                          const base_abstract_domain_t &v2) {
      return v1 || v2;
    };
    auto base_uf_dom_op = [](const base_uf_domain_t &v1,
                             const base_uf_domain_t &v2) { return v1 || v2; };

    mru_region_domain_t res(std::move(do_join_or_widening(
        *this, o, false /*is widen*/, base_dom_op, base_uf_dom_op)));

    CRAB_LOG("mru-region", crab::outs() << res << "\n");
    return res;
  }

  mru_region_domain_t widening_thresholds(
      const mru_region_domain_t &o,
      const iterators::thresholds<number_t> &thresholds) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");

    // Trivial cases first: we don't cover cases where one operand is
    // top because is_top() calls the base domain which we don't know
    // whether it will perform some normalization or not.
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    }

    CRAB_LOG("mru-region", crab::outs() << "Widening with threshold " << *this
                                        << " and " << o << "=");

    auto base_dom_op = [&thresholds](const base_abstract_domain_t &v1,
                                     const base_abstract_domain_t &v2) {
      return v1.widening_thresholds(v2, thresholds);
    };
    auto base_uf_dom_op = [&thresholds](const base_uf_domain_t &v1,
                                        const base_uf_domain_t &v2) {
      return v1.widening_thresholds(v2, thresholds);
    };

    mru_region_domain_t res(std::move(do_join_or_widening(
        *this, o, false /*is widen*/, base_dom_op, base_uf_dom_op)));

    CRAB_LOG("mru-region", crab::outs() << res << "\n");
    return res;
  }

  mru_region_domain_t operator&&(const mru_region_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.narrowing");
    crab::ScopedCrabStats __st__(domain_name() + ".narrowing");

    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (o.is_bottom() || is_top()) {
      return o;
    }

    CRAB_LOG("mru-region", crab::outs()
                               << "Meet " << *this << " and " << o << "=");

    auto base_dom_op = [](const base_abstract_domain_t &v1,
                          const base_abstract_domain_t &v2) {
      return v1 && v2;
    };
    auto base_uf_dom_op = [](const base_uf_domain_t &v1,
                             const base_uf_domain_t &v2) { return v1 && v2; };

    mru_region_domain_t res(std::move(do_meet_or_narrowing(
        *this, o, false /*is widen*/, base_dom_op, base_uf_dom_op)));

    CRAB_LOG("mru-region", crab::outs() << res << "\n");
    return res;
  }

  // x := e operates on m_mem
  void assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    if (!is_bottom()) {
      auto b_e = rename_linear_expr(e);
      m_mem.assign(get_or_insert_gvars(x).get_var(), b_e);
    }
  }

  // x := e operates on m_cache_lines
  void cache_lines_assign(const variable_t &x, const linear_expression_t &e) {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");
    // TODO: must check x and e are region-based assignment
    if (!is_bottom()) {
      auto b_e = rename_linear_expr(e);
      const ghost_variables_t &new_gvars = get_or_insert_gvars(x);
      m_cache_var_map.insert({x, new_gvars});
      m_cache_lines.assign(new_gvars.get_var(), b_e);
    }
  }

  // x := e operates on m_regs
  void regs_assign(const variable_t &x, const linear_expression_t &e) {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");
    // precondition: x is a scalar variable and e is a region variable
    if (!is_bottom()) {
      auto b_e = rename_linear_expr(e);
      const ghost_variables_t &new_gvars = get_or_insert_gvars(x);
      m_cache_var_map.insert({x, new_gvars});
      m_regs.assign(new_gvars.get_var(), b_e);
    }
  }

  // x := e operates on m_addrs
  void addrs_assign(const variable_t &x, const variable_t &y) {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");
    // precondition: x is a scalar variable and e is a region variable
    if (!is_bottom()) {
      // auto b_e = rename_linear_expr(e);
      const ghost_variables_t &x_gvars = get_or_insert_gvars(x);
      const ghost_variables_t &y_gvars = get_or_insert_gvars(y);
      if (!x.name().str().compare("C_base")) {
        uninterpreted_symbol_t symb = base_uf_domain_t::make_uf(
            m_uninterpreted_symbol_index);
        m_uninterpreted_symbol_index++;
        m_base = std::make_pair(x_gvars, symb);
        get_or_insert_uninterpreted_symbol(y, symb);
      } else {
        if (!m_addr_var_map.contains(x)) {
          auto x_symb = get_or_insert_uninterpreted_symbol(x);
          get_or_insert_uninterpreted_symbol(y, x_symb);
        }
      }
      m_addrs.assign(x_gvars.get_var(), y_gvars.get_var());
    }
  }

  // add all constraints \in csts
  void operator+=(const linear_constraint_system_t &csts) override {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");

    if (!is_bottom()) {
      auto b_csts = rename_linear_cst_sys(csts);
      m_mem += b_csts;
      m_is_bottom = m_mem.is_bottom();
    }
  }

  void add_cons_into_cache_lines(const linear_constraint_system_t &csts) {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");

    if (!is_bottom()) {
      auto b_csts = rename_linear_cst_sys(csts);
      m_cache_lines += b_csts;
    }
  }

  // lhs := rhs
  void assign_bool_cst(const variable_t &lhs,
                       const linear_constraint_t &rhs) override {
    crab::CrabStats::count(domain_name() + ".count.assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_cst");

    if (!is_bottom()) {
      auto b_rhs = rename_linear_cst(rhs);
      m_mem.assign_bool_cst(get_or_insert_gvars(lhs).get_var(), b_rhs);

      if (crab_domain_params_man::get().region_tag_analysis()) {
        merge_tags(lhs, rhs.variables());
      }
    }
  }

  /***************** Regions and reference operations *****************/

  // Initialize a region
  void region_init(const variable_t &rgn) override {
    crab::CrabStats::count(domain_name() + ".count.region_init");
    crab::ScopedCrabStats __st__(domain_name() + ".region_init");

    ERROR_IF_NOT_REGION(rgn, __LINE__);

    if (is_bottom()) {
      return;
    }

    // check whether region has been initialized
    small_range count_num = m_rgn_counting_dom[rgn];
    if (count_num <= small_range::oneOrMore()) {
      CRAB_ERROR(domain_name(), "::init_region: ", rgn,
                 " cannot be initialized twice");
    }

    // Set to zero the number of references
    m_rgn_counting_dom.set(rgn, small_range::zero());

    // No stores in the region
    m_rgn_init_dom.set(rgn, boolean_value::get_false());

    if (is_tracked_unknown_region(rgn)) {
      m_rgn_type_dom.set(rgn, rgn.get_type());
    }

    if (!(rgn.get_type().is_unknown_region())) {
      // Assign ghost variables to rgn for modeling its content
      ghost_variables_t::create(m_alloc, rgn.get_type(), __LINE__);
      crab::CrabStats::count(domain_name() +
                             ".count.region_init.typed_regions");
    } else {
      // We cannot assign ghost variables if the region is unknown
      crab::CrabStats::count(domain_name() +
                             ".count.region_init.unknown_regions");
    }
    CRAB_LOG("mru-region", crab::outs()
                               << "After region_init(" << rgn << ":"
                               << rgn.get_type() << ")=" << *this << "\n";);
  }

  const uninterpreted_symbol_t &get_or_insert_uninterpreted_symbol(
      const variable_t &v, boost::optional<variable_t> equiv_v = boost::none) {
    return get_or_insert_uninterpreted_symbol(v, m_addr_var_map, equiv_v);
  }

  const uninterpreted_symbol_t &get_or_insert_uninterpreted_symbol(
      const variable_t &v, uf_var_map_t &addr_var_map,
      boost::optional<variable_t> equiv_v) {
    if (addr_var_map.contains(v)) {
      return addr_var_map.get(v);
    } else {
      if (equiv_v) {
        auto res = addr_var_map.join(v, *equiv_v);
        return res.first->second;
      } else {
        uninterpreted_symbol_t symb = base_uf_domain_t::make_uf(
            m_uninterpreted_symbol_index);
        m_uninterpreted_symbol_index++;
        auto res = addr_var_map.set(v, symb);
        return res.first->second;
      }
    }
  }

  // Create a new reference ref associated with as within region
  void ref_make(const variable_t &ref, const variable_t &rgn,
                /* size of the allocation in bytes */
                const variable_or_constant_t &size,
                /* identifier for the allocation site */
                const allocation_site &as) override {
    crab::CrabStats::count(domain_name() + ".count.ref_make");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_make");

    ERROR_IF_NOT_REGION(rgn, __LINE__);
    ERROR_IF_NOT_REF(ref, __LINE__);

    if (is_bottom()) {
      return;
    }

    // Update region counting
    auto num_refs = m_rgn_counting_dom[rgn];
    m_rgn_counting_dom.set(rgn, num_refs.increment());

    // Assign ghost variables to ref
    ghost_variables_t ref_gvars = get_or_insert_gvars(ref);

    // initialize ghost variables
    if (ref_gvars.has_offset_and_size()) {
      ref_gvars.get_offset_and_size().init(m_mem,
                                           rename_variable_or_constant(size));
    }

    // Assign uf symbol inside address domain
    uninterpreted_symbol_t symb = get_or_insert_uninterpreted_symbol(ref);
    m_addrs.set(ref_gvars.get_var(), symb);

    CRAB_LOG("mru-region", crab::outs()
                               << "After ref_make(" << ref << "," << rgn << ":"
                               << rgn.get_type() << "," << size << "," << as
                               << ")=" << *this << "\n";);
  }

  // Read the content of reference ref within rgn. The content is
  // stored in res.
  void ref_load(const variable_t &ref, const variable_t &rgn,
                const variable_t &res) override {
    crab::CrabStats::count(domain_name() + ".count.ref_load");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_load");

    ERROR_IF_NOT_REGION(rgn, __LINE__);
    ERROR_IF_ARRAY_REGION(rgn, __LINE__);
    ERROR_IF_NOT_REF(ref, __LINE__);

    if (is_bottom()) {
      return;
    }

    const ghost_variables_t &gvars_res = get_or_insert_gvars(res);

    if (is_null_ref(ref).is_true()) {
      CRAB_LOG("region-load", CRAB_WARN(domain_name(), "::ref_load: reference ",
                                        ref, " is null."););
      // set_to_bottom();
      crab::CrabStats::count(domain_name() +
                             ".count.ref_load.skipped.null_reference");
      gvars_res.forget(m_mem);
      m_regs.forget({gvars_res.get_var()});
      return;
    }

    // TODO: we now only assume the region is tracked. Either dynamic tracking
    // is not supported.

    // This is the core process during ref_store:
    //  Each load/store, we need to check whether the current reference
    //  refers to the same object stored in cache or not.
    // If the rgn belongs to some object that is stored in cache,
    // the following method does nothing.
    // Otherwise, we firstly commit the cache if the cache is dirty;
    // then update the cache by copying all regions' invariants into
    // m_cache_lines.
    invalidate_cache_if_miss(rgn, ref);

    if (auto region_gvars_opt = get_gvars(rgn)) {
      // if region exists
      // assign lhs := rgn
      m_regs.assign(gvars_res.get_var(), (*region_gvars_opt).get_var());
      m_cache_var_map.insert({res, gvars_res});
    } else {
      // if not, forget res
      gvars_res.forget(m_mem);
      m_regs.forget({gvars_res.get_var()});
    }

    CRAB_LOG("mru-region-load", crab::outs() << "After " << res << ":="
                                             << "ref_load(" << rgn << ":"
                                             << rgn.get_type() << "," << ref
                                             << ":" << ref.get_type()
                                             << ")=" << *this << "\n";);
  }

  // Write the content of val to the address pointed by ref in region.
  void ref_store(const variable_t &ref, const variable_t &rgn,
                 const variable_or_constant_t &val) override {
    crab::CrabStats::count(domain_name() + ".count.ref_store");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_store");

    ERROR_IF_NOT_REGION(rgn, __LINE__);
    ERROR_IF_ARRAY_REGION(rgn, __LINE__);
    ERROR_IF_NOT_REF(ref, __LINE__);
    // type checking: whether rgn and val are consistent
    if ((rgn.get_type().is_bool_region() && !val.get_type().is_bool()) ||
        (rgn.get_type().is_integer_region() && !val.get_type().is_integer()) ||
        (rgn.get_type().is_real_region() && !val.get_type().is_real()) ||
        (rgn.get_type().is_reference_region() &&
         !val.get_type().is_reference())) {
      CRAB_ERROR(domain_name(), "::ref_store: type of value ", val, " (",
                 val.get_type(), ") is not compatible with region ", rgn, " (",
                 rgn.get_type(), ")");
    }

    if (is_bottom()) {
      return;
    }

    // forget the ghost variables in the base domain
    auto forget_region_ghost_vars = [this](const variable_t &rgn) {
      if (boost::optional<ghost_variables_t> gvars = get_gvars(rgn)) {
        (*gvars).forget(m_mem);
        (*gvars).forget(m_cache_lines);
        m_regs.forget({(*gvars).get_var()});
      }
    };

    if (is_null_ref(ref).is_true()) {
      CRAB_LOG("mru-region-store",
               CRAB_WARN(domain_name(), "::ref_store: reference ", ref,
                         " is null. "););
      // set_to_bottom();
      crab::CrabStats::count(domain_name() +
                             ".count.ref_store.skipped.null_reference");
      forget_region_ghost_vars(rgn);
      return;
    }
    // TODO: we now only assume the region is tracked. Either dynamic tracking
    // is not supported.

    // This is the core process during ref_store:
    //  Each load/store, we need to check whether the current reference
    //  refers to the same object stored in cache or not.
    // If the rgn belongs to some object that is stored in cache,
    // the following method does nothing.
    // Otherwise, we firstly commit the cache if the cache is dirty;
    // then update the cache by copying all regions' invariants into
    // m_cache_lines.
    invalidate_cache_if_miss(rgn, ref);
    ghost_variables_t rgn_gvars = get_or_insert_gvars(rgn);

    // assign lhs:= rgn in regs domain
    // This step will always perform a strong update
    if (val.get_type().is_bool()) {
      if (val.is_constant()) {
        // if val is a constant, we directly update the cache line
        m_cache_lines.assign_bool_cst(
            rgn_gvars.get_var(),
            (val.is_bool_true() ? base_linear_constraint_t::get_true()
                                : base_linear_constraint_t::get_false()));
      } else { // val is a variable
        m_regs.assign_bool_var(
            rgn_gvars.get_var(),
            get_or_insert_gvars(val.get_variable()).get_var(), false);
        const ghost_variables_t &val_gvars =
            get_or_insert_gvars(val.get_variable());
        m_cache_var_map.insert({val.get_variable(), val_gvars});
      }
    } else if (val.get_type().is_integer() || val.get_type().is_real() ||
               val.get_type().is_reference()) {

      if (val.is_constant()) {
        // if val is a constant, we directly update the cache line
        m_cache_lines.assign(rgn_gvars.get_var(), val.get_constant());
        if (val.get_type().is_reference()) {
          // this one assume rgn is a reference region
          if (rgn_gvars.has_offset_and_size()) {
            // Storing NULL
            assert(val.get_constant() == number_t(0));
            rgn_gvars.get_offset_and_size().forget(m_mem);
          }
        }
      } else { // val is variable
        ghost_variables_t val_gvars = get_or_insert_gvars(val.get_variable());
        m_regs.assign(rgn_gvars.get_var(), val_gvars.get_var());

        if (val.get_type().is_reference()) {
          if (rgn_gvars.has_offset_and_size() &&
              val_gvars.has_offset_and_size()) {
            rgn_gvars.get_offset_and_size().assign(
                m_mem, val_gvars.get_offset_and_size());
          } else if (rgn_gvars.has_offset_and_size()) {
            rgn_gvars.get_offset_and_size().forget(m_mem);
          }
        }
        m_cache_var_map.insert({val.get_variable(), val_gvars});
      }
    } else {
      CRAB_ERROR(domain_name(), "::ref_store: unsupported type ",
                 val.get_type());
    }

    m_dirty = boolean_value::get_true(); // C_dirty = false

    CRAB_LOG("mru-region-store",
             crab::outs() << "After ref_store(" << rgn << ":" << rgn.get_type()
                          << "," << ref << ":" << ref.get_type() << "," << val
                          << ":" << val.get_type() << ")=" << *this << "\n";);
  }

  // Create a new reference ref2 to region rgn2.
  // The reference ref2 is created by adding offset to ref1.
  void ref_gep(const variable_t &ref1, const variable_t &rgn1,
               const variable_t &ref2, const variable_t &rgn2,
               const linear_expression_t &offset) override {
    crab::CrabStats::count(domain_name() + ".count.ref_gep");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_gep");

    ERROR_IF_NOT_REGION(rgn1, __LINE__);
    ERROR_IF_NOT_REGION(rgn2, __LINE__);
    ERROR_IF_NOT_REF(ref1, __LINE__);
    ERROR_IF_NOT_REF(ref2, __LINE__);

    if (is_bottom()) {
      return;
    }

    // ghost variable to ref1
    ghost_variables_t ref1_gvars = get_or_insert_gvars(ref1);
    uninterpreted_symbol_t ref1_symb = get_or_insert_uninterpreted_symbol(ref1);

    // ref2
    ghost_variables_t ref2_gvars = get_or_insert_gvars(ref2);
    get_or_insert_uninterpreted_symbol(ref2, ref1_symb);

    // assign pointer value for ref_2
    m_mem.assign(ref2_gvars.get_var(), ref1_gvars.get_var() + offset);

    // assign base address equality
    // ref2.base == ref1.base
    m_addrs.assign(ref2_gvars.get_var(), ref1_gvars.get_var());
    // m_addrs.set(ref2_gvars.get_var(), ref1_symb);

    CRAB_LOG("mru-region", crab::outs()
                               << "After (" << rgn2 << "," << ref2
                               << ") := ref_gep(" << rgn1 << "," << ref1
                               << " + " << offset << ")=" << *this << "\n";);
  }

  // FIXME: The followings are UNDEFINED METHODS
  bool get_allocation_sites(const variable_t &ref,
                            std::vector<allocation_site> &out) override {
    return false;
  }

  bool get_tags(const variable_t &rgn, const variable_t &ref,
                std::vector<uint64_t> &out) override {
    return false;
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {}
  // x := y op k
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {}

  // x := y op z
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {}
  // x := y op k
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {}
  // dst := src
  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {}

  // if(cond) lhs := e1 else lhs := e2
  void select(const variable_t &lhs, const linear_constraint_t &cond,
              const linear_expression_t &e1,
              const linear_expression_t &e2) override {}

  /********************** Boolean operations **********************/
  void assign_bool_ref_cst(const variable_t &lhs,
                           const reference_constraint_t &rhs) override {}
  // lhs := not(rhs) if is_not_rhs
  // lhs := rhs      otherwise
  void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                       bool is_not_rhs) override {}
  // x := y op z
  void apply_binary_bool(bool_operation_t op, const variable_t &x,
                         const variable_t &y, const variable_t &z) override {}
  // assume(not(v)) if is_negated
  // assume(v)      otherwise
  void assume_bool(const variable_t &v, bool is_negated) override {}

  // if(cond) lhs := b1 else lhs := b2
  // lhs, cond, b1, and b2 are boolean variables
  void select_bool(const variable_t &lhs, const variable_t &cond,
                   const variable_t &b1, const variable_t &b2) override {}

  /********************** Array operations **********************/
  // make a fresh array with contents a[j] initialized to val such that
  // j \in [lb_idx,ub_idx] and j % elem_size == val.
  // elem_size is in bytes.
  void array_init(const variable_t &a, const linear_expression_t &elem_size,
                  const linear_expression_t &lb_idx,
                  const linear_expression_t &ub_idx,
                  const linear_expression_t &val) override {}
  // lhs := a[i] where elem_size is in bytes
  void array_load(const variable_t &lhs, const variable_t &a,
                  const linear_expression_t &elem_size,
                  const linear_expression_t &i) override {}
  // a[i] := val where elem_size is in bytes
  void array_store(const variable_t &a, const linear_expression_t &elem_size,
                   const linear_expression_t &i, const linear_expression_t &val,
                   bool is_strong_update) override {}
  // forall i<=k<j and k % elem_size == 0 :: a[k] := val.
  // elem_size is in bytes
  void array_store_range(const variable_t &a,
                         const linear_expression_t &elem_size,
                         const linear_expression_t &i,
                         const linear_expression_t &j,
                         const linear_expression_t &val) override {}
  // forall i :: a[i] := b[i]
  void array_assign(const variable_t &a, const variable_t &b) override {}

  /********************* Regions and reference operations *********************/
  // Make a copy of a region
  void region_copy(const variable_t &lhs_reg,
                   const variable_t &rhs_reg) override {}
  // Cast between regions of different types
  void region_cast(const variable_t &src_reg,
                   const variable_t &dst_reg) override {}
  // Remove a reference ref within region reg
  void ref_free(const variable_t &reg, const variable_t &ref) override {}
  // Treat memory pointed by ref  as an array and perform an array load.
  void ref_load_from_array(const variable_t &lhs, const variable_t &ref,
                           const variable_t &region,
                           const linear_expression_t &index,
                           const linear_expression_t &elem_size) override {}
  // Treat region as an array and perform an array store.
  void ref_store_to_array(const variable_t &ref, const variable_t &region,
                          const linear_expression_t &index,
                          const linear_expression_t &elem_size,
                          const linear_expression_t &val) override {}
  // Add constraints between references
  void ref_assume(const reference_constraint_t &cst) override {}
  // Convert a reference to an integer variable
  void ref_to_int(const variable_t &reg, const variable_t &ref,
                  const variable_t &int_var) override {}
  // Convert an integer variable to a reference
  void int_to_ref(const variable_t &int_var, const variable_t &reg,
                  const variable_t &ref) override {}
  // if (cond) ref_gep(ref1, rgn1, lhs_ref, lhs_rgn, 0) else
  //           ref_gep(ref2, rgn2, lhs_ref, lhs_rgn, 0)
  // cond is a boolean variable
  void select_ref(const variable_t &lhs_ref, const variable_t &lhs_rgn,
                  const variable_t &cond, const variable_or_constant_t &ref1,
                  const boost::optional<variable_t> &rgn1,
                  const variable_or_constant_t &ref2,
                  const boost::optional<variable_t> &rgn2) override {}

  /********************** Backward numerical operations **********************/
  // x = y op z
  // Substitute x with y op z in the abstract value
  // The result is meet with invariant.
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const mru_region_domain_t &invariant) override {}
  // x = y op k
  // Substitute x with y op k in the abstract value
  // The result is meet with invariant.
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t k,
                      const mru_region_domain_t &invariant) override {}
  // x = e
  // Substitute x with e in the abstract value
  // The result is meet with invariant.
  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const mru_region_domain_t &invariant) override {}

  /********************** Backward boolean operations **********************/
  void backward_assign_bool_cst(const variable_t &lhs,
                                const linear_constraint_t &rhs,
                                const mru_region_domain_t &invariant) override {
  }
  void
  backward_assign_bool_ref_cst(const variable_t &lhs,
                               const reference_constraint_t &rhs,
                               const mru_region_domain_t &invariant) override {}
  void backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                                bool is_not_rhs,
                                const mru_region_domain_t &invariant) override {
  }
  void
  backward_apply_binary_bool(bool_operation_t op, const variable_t &x,
                             const variable_t &y, const variable_t &z,
                             const mru_region_domain_t &invariant) override {}

  /********************** Backward array operations **********************/
  void backward_array_init(const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &lb_idx,
                           const linear_expression_t &ub_idx,
                           const linear_expression_t &val,
                           const mru_region_domain_t &invariant) override {}
  void backward_array_load(const variable_t &lhs, const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const mru_region_domain_t &invariant) override {}
  void backward_array_store(const variable_t &a,
                            const linear_expression_t &elem_size,
                            const linear_expression_t &i,
                            const linear_expression_t &v, bool is_strong_update,
                            const mru_region_domain_t &invariant) override {}
  void backward_array_store_range(
      const variable_t &a, const linear_expression_t &elem_size,
      const linear_expression_t &i, const linear_expression_t &j,
      const linear_expression_t &v,
      const mru_region_domain_t &invariant) override {}
  void backward_array_assign(const variable_t &a, const variable_t &b,
                             const mru_region_domain_t &invariant) override {}

  /********************** Miscellaneous operations **********************/
  // Forget v
  void operator-=(const variable_t &v) override {}

  // Return an interval with the possible values of v if such notion
  // exists in the abstract domain.
  interval_t operator[](const variable_t &v) override {
    return interval_t::bottom();
  }

  // Normalize the abstract domain if such notion exists.
  void normalize() override {}

  // Reduce the size of the abstract domain representation.
  void minimize() override {}

  // Make a new copy of var without relating var with new_var
  void expand(const variable_t &var, const variable_t &new_var) override {}

  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const mru_region_domain_t &invariant) override {}

  // FIXME: The above methods are UNDEFINED METHODS

  void forget(const variable_vector_t &variables) override {
    if (is_bottom() || is_top()) {
      return;
    }

    std::vector<base_variable_t> base_vars;
    base_vars.reserve(variables.size());
    for (auto &v : variables) {
      if (v.get_type().is_region()) {
        m_rgn_counting_dom -= v;
        m_rgn_init_dom -= v;
        if (is_tracked_unknown_region(v)) {
          m_rgn_type_dom -= v;
        }
        if (crab_domain_params_man::get().region_deallocation()) {
          m_rgn_dealloc_dom.forget(v);
        }
        m_obj_rgn_map.forget(v);
      }
      if (crab_domain_params_man::get().region_tag_analysis()) {
        m_tag_env -= v;
      }
      auto it = m_var_map.find(v);
      if (it != m_var_map.end()) {
        it->second.add(base_vars);
        it->second.remove_rev_varmap(m_rev_var_map);
        m_var_map.erase(it);
      }
    }
    // forget variable into each subdomain
    m_mem.forget(base_vars);
    m_cache_lines.forget(base_vars);
    m_regs.forget(base_vars);
    m_addrs.forget(base_vars);
  }

  void project(const variable_vector_t &variables) override {
    crab::CrabStats::count(domain_name() + ".count.project");
    crab::ScopedCrabStats __st__(domain_name() + ".project");

    if (is_bottom() || is_top()) {
      return;
    }

    variable_vector_t sorted_variables(variables.begin(), variables.end());
    std::sort(sorted_variables.begin(), sorted_variables.end());

    // -- project in the base domain
    std::vector<base_variable_t> base_vars;
    base_vars.reserve(variables.size());
    for (auto const &v : variables) {
      if (v.get_type().is_region() && !is_tracked_region(v)) {
        continue;
      }
      get_or_insert_gvars(v).add(base_vars);
    }

    // project mem domain
    m_mem.project(base_vars);
    // project cache_lines domain
    m_cache_lines.project(base_vars);
    // project regs domain
    m_regs.project(base_vars);
    // project addrs domain
    m_addrs.project(base_vars);

    // -- update m_var_map and m_rev_var_map
    std::vector<variable_t> var_map_to_remove;
    var_map_to_remove.reserve(m_var_map.size());
    for (auto &kv : m_var_map) {
      if (!std::binary_search(sorted_variables.begin(), sorted_variables.end(),
                              kv.first)) {
        var_map_to_remove.push_back(kv.first);
        kv.second.remove_rev_varmap(m_rev_var_map);
      }
    }
    for (auto &v : var_map_to_remove) {
      m_var_map.erase(v);
      m_cache_var_map.erase(v);
    }

    m_rgn_counting_dom.project(sorted_variables);
    m_rgn_init_dom.project(sorted_variables);
    if (!crab_domain_params_man::get().region_skip_unknown_regions()) {
      m_rgn_type_dom.project(sorted_variables);
    }
    if (crab_domain_params_man::get().region_deallocation()) {
      m_rgn_dealloc_dom.project(sorted_variables);
    }
    if (crab_domain_params_man::get().region_tag_analysis()) {
      m_tag_env.project(sorted_variables);
    }
    m_obj_rgn_map.project(sorted_variables);
  }

  // NOTE: region might init without sorting but rename method required
  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    crab::CrabStats::count(domain_name() + ".count.rename");
    crab::ScopedCrabStats __st__(domain_name() + ".rename");

    if (is_bottom() || is_top()) {
      return;
    }

    if (from.size() != to.size()) {
      CRAB_ERROR(domain_name(), "::rename different lengths");
    }

    // update var_map and rev_var_map
    // extract base_variable_vector for from and to
    base_variable_vector_t old_base_vars, new_base_vars;
    for (unsigned i = 0; i < from.size(); ++i) {
      variable_t old_var = from[i];
      auto search = m_var_map.find(old_var);
      if (search == m_var_map.end()) { // NOT FOUND
        continue;
      }
      ghost_variables_t old_ghost_vars = search->second;
      // append base_var, [offset, size] into old_base_vars
      old_ghost_vars.add(old_base_vars);
      // remove from rev_var_map
      old_ghost_vars.remove_rev_varmap(m_rev_var_map);
      // remove from var_map
      m_var_map.erase(old_var);

      variable_t new_var = to[i];
      // create/get new ghost variable for new_var
      const ghost_variables_t &new_gvars = get_or_insert_gvars(new_var);
      // append base_var, [offset, size] into new_base_vars
      new_gvars.add(new_base_vars);
    }

    // rename mem domain
    m_mem.rename(old_base_vars, new_base_vars);
    // rename cache_lines domain
    m_cache_lines.rename(old_base_vars, new_base_vars);
    // rename regs domain
    m_regs.rename(old_base_vars, new_base_vars);
    // rename addrs domain
    m_addrs.rename(old_base_vars, new_base_vars);

    // rename the rest of environments
    m_rgn_counting_dom.rename(from, to);
    m_rgn_init_dom.rename(from, to);
    if (!crab_domain_params_man::get().region_skip_unknown_regions()) {
      m_rgn_type_dom.rename(from, to);
    }
    if (crab_domain_params_man::get().region_deallocation()) {
      m_rgn_dealloc_dom.rename(from, to);
    }
    if (crab_domain_params_man::get().region_tag_analysis()) {
      m_tag_env.rename(from, to);
    }
    m_obj_rgn_map.rename(from, to);
  }

  /* begin intrinsics operations */
  void intrinsic(std::string name, const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    //=================================================================//
    //       Special intrinsics supported by the region domain
    //=================================================================//
    // ---DSA region analysis---
    //      This analysis indicates which regions might belong to the
    //      same memory object. The intrinstic is added only if object
    //      has more than one field.
    // TODO: need to support other analysis
    auto error_if_not_rgn = [&name](const variable_or_constant_t &x) {
      if (!x.is_variable() || !x.get_type().is_region()) {
        // the input vector should only contains region variables
        CRAB_ERROR("Intrinsics ", name, " parameter ", x,
                   " should be a region");
      }
    };

    if (is_bottom()) {
      return;
    }

    if (name == "regions_from_memory_object") {
      // pass region varaibles into object region map
      // the intrinsics is only added if the object has more than one region.
      assert(inputs.size() >= 1);
      for (int i = 0, sz = inputs.size(); i < sz; ++i) {
        error_if_not_rgn(inputs[i]); // this should not happen
        if (inputs[i].get_type().is_unknown_region()) {
          // get_or_insert_gvars does not support unknown regions so we bail
          // out. The base domain shouldn't care about regions anyway.
          return;
        }
        // update object region map
        if (i == 0) { // also make it as merging representative
          m_obj_rgn_map.set(inputs[i].get_variable(),
                            boolean_value::get_false());
        } else {
          assert(m_obj_rgn_map.contains(inputs[0].get_variable()));
          m_obj_rgn_map.set(inputs[i].get_variable(),
                            boolean_value::get_false());
          assert(m_obj_rgn_map.contains(inputs[i].get_variable()));
          m_obj_rgn_map.join(inputs[0].get_variable(),
                             inputs[i].get_variable());
        }
      }
    }
  }

  // Return a 3-valued boolean. If true then ref is definitely
  // null. If false then ref is definitely non-null. Otherwise, we do
  // not know.
  boolean_value is_null_ref(const variable_t &ref) override {
    if (is_bottom()) {
      return boolean_value::bottom();
    }

    if (!ref.get_type().is_reference()) {
      return boolean_value::get_false();
    }

    if (auto gvars_opt = get_gvars(ref)) {
      interval_t ival = m_mem[(*gvars_opt).get_var()];
      number_t zero(0);

      if (!(interval_t(zero) <= ival)) {
        return boolean_value::get_false();
      }

      boost::optional<number_t> x = ival.lb().number();
      boost::optional<number_t> y = ival.ub().number();
      if (x && y && *x == zero && *y == zero) {
        return boolean_value::get_true();
      }
    }

    return boolean_value::top();
  }

  // Convert the abstract state into a conjunction of linear constraints.
  // Only integer or boolean variables in memory are exposed.
  //
  // Variables that shadow memory regions and reference variables are
  // ignored.
  linear_constraint_system_t to_linear_constraint_system() const override {
    if (is_bottom()) {
      return linear_constraint_t::get_false();
    } else if (is_top()) {
      return linear_constraint_t::get_true();
    } else {
      linear_constraint_system_t out_csts;
      const bool ignore_references = true;
      for (base_linear_constraint_t cst : m_mem.to_linear_constraint_system()) {
        if (boost::optional<linear_constraint_t> out_cst =
                rev_rename_linear_cst(cst, ignore_references)) {
          out_csts += *(out_cst);
        }
      }
      return out_csts;
    }
  }

  // Convert the abstract state into a disjunction of conjunction
  // of linear constraints.
  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    CRAB_ERROR(domain_name(), "::to_disjunctive_linear_constraint_system not "
                              "implemented");
  }

  std::string domain_name() const override {
    return "MRURegionDomain(" + m_cache_lines.domain_name() + ", " +
           m_regs.domain_name() + ")";
  }

  void write(crab_os &o) const override {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "{}";
    } else {
      CRAB_LOG("mru-region-print", o << "Mem= " << m_mem << ",\n"
                                     << "( CacheLine= " << m_cache_lines
                                     << ",\n"
                                     << "Regs= " << m_regs << ",\n"
                                     << "Addrs= " << m_addrs << ",\n"
                                     << "C_used= " << m_used << ",\n"
                                     << "C_dirty= " << m_dirty << " )\n";);

      std::unordered_map<std::string, std::string> renaming_map;
      auto get_type_fn = [](const variable_t &v) -> variable_type {
        return variable_type(variable_type_kind::UNK_TYPE);
      };
      // Assigns a string name to each ghost variable
      ghost_variables_t::mk_renaming_map(m_rev_var_map, get_type_fn,
                                         renaming_map);
      if (m_base) {
        base_variable_t c_base = (*m_base).first.get_var();
        renaming_map[c_base.name().str()] = "C_base";
      }
      m_alloc.add_renaming_map(renaming_map);
      o << "Mem= " << m_mem << ",\n"
        << "( CacheLine= " << m_cache_lines << ",\n"
        << "Regs= " << m_regs << ",\n"
        << "Addrs= " << m_addrs << ",\n"
        << "C_used= " << m_used << ",\n"
        << "C_dirty= " << m_dirty << " )\n";
      m_alloc.clear_renaming_map();
    }
  }
  /**------------- End domain API definitions -------------------**/

  /**------------- Begin MRU API definitions -------------------**/

  // The expand performs a copy of region vars from mem into cache_lines
  //  i.e. by given a most recently used object o, we copy all regions with o
  //  in mem into the cache:
  // For instance, o used V1, V2:
  //    mem: { 1 <= V1 <= 2; V1 < V2; ... }
  //    cache_lines: top // means nothing stored
  // The cache line will be:
  //    cache_lines: { 1 <= V1 <= 2; V1 < V2; }
  void expand(const variable_t &rgn) {
    variable_vector_t Vs_for_new_obj;

    // perform constructing a vector of region variables for new mru object
    if (!m_obj_rgn_map.contains(rgn)) {
      // this should not happen
      CRAB_ERROR(domain_name(), "::expand the mru object should be found");
    }
    m_obj_rgn_map.get_all_equiv_variables(rgn, Vs_for_new_obj);

    // expand does not need to perform base variable renaming
    mru_region_domain_t tmp = mru_region_domain_t(*this);
    tmp.project(Vs_for_new_obj);

    // perform expand
    base_abstract_domain_t res(m_cache_lines & tmp.m_mem);
    std::swap(m_cache_lines, res);

    // update cache_var map
    for (unsigned i = 0; i < Vs_for_new_obj.size(); ++i) {
      const ghost_variables_t &ghost_var =
          get_or_insert_gvars(Vs_for_new_obj[i]);
      m_cache_var_map.insert({Vs_for_new_obj[i], ghost_var});
    }

    CRAB_LOG("mru-region", crab::outs() << *this << "\n");
  }

  // The cache is used by given an allocated object.
  //  The fold operation performs smash all regions (fields) into the summary
  //  object in mem. That is, for each region used in cache line, join with
  //  region in mem.
  //  In addition, we also update regs used in cache via reduction.
  void fold() {
    // Put together cache lines and regs
    // Step 1.
    //  Perform regs: { r == Cj ... }
    //  constrain with cache_lines: { Ci <= Cj <= Ck ... }, i != j != k
    // the reduced cache_lines will be:
    //    { Ci <= Cj <= Ck; Ci <= r <= Ck; ...  }
    base_abstract_domain_t reduced = cache_lines_regs_reduce();

    // Step 2.
    //  Our design is joining cache line with memory.
    //  We need the projection before join because we only fold regions & regs
    //  We do not affect other properties in mem.
    variable_vector_t vars_cache;
    vars_in_cache(vars_cache); // vars for cache lines + registers
    mru_region_domain_t A = mru_region_domain_t(*this);
    A.project(vars_cache); // memory object represented in the cache
    A.m_mem |= reduced;    // join cache with mem,
                           // this step might lost other variables

    // Step 3.
    // So far we performed the fold operation, but the rest properties need to
    // put it back. This required we project all variables used in cache
    // Finally, we perform meet to recover.
    mru_region_domain_t B = mru_region_domain_t(*this);
    variable_vector_t mem_vars;
    vars(mem_vars); // All variables in mem
    variable_vector_t vars_cache_lines;
    region_vars_in_cache(vars_cache_lines); // All region variables in cache
    B.project(vec_minus(mem_vars,
                        vars_cache_lines)); // get all properties except those
                                            // relevant to the cache lines
    mru_region_domain_t res = A & B; // put back the rest of the properties

    CRAB_LOG("mru-region", crab::outs() << res << "\n");
    std::swap(*this, res);
  }

  /**------------- End MRU API definitions -------------------**/

}; // end class mru_region_domain

template <typename Params>
struct abstract_domain_traits<mru_region_domain<Params>> {
  using number_t = typename Params::number_t;
  using varname_t = typename Params::varname_t;
};

} // end namespace domains
} // end namespace crab