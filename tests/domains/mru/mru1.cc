#include "../../common.hpp"
#include "../../program_options.hpp"
#include "../../mru_crab_dom.hpp"


using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::mru_domain_impl;

namespace {
  using variable_vector_t = std::vector<z_var>;
  using variable_t = z_var;
}

void simulate_obj_map_update_per_object(z_mru_rgn_zones_t &mru_mem, 
                                        variable_vector_t rgn_vars) {
  variable_t representative = rgn_vars[0];
  for (int i = 0, sz = rgn_vars.size(); i < sz; ++i) {
    // update object region map
    if (i == 0) { // also make it as merging representative
      mru_mem.m_obj_rgn_map.set(rgn_vars[i], boolean_value::get_false());
    }
    else {
      mru_mem.m_obj_rgn_map.set(rgn_vars[i], boolean_value::get_false());
      mru_mem.m_obj_rgn_map.join(rgn_vars[0], rgn_vars[i]);
    }
  }
}

z_cfg_t *prog(z_mru_rgn_zones_t &init, variable_factory_t &vfac) {
  ////
  // Building the CFG
  ////

  // Definining program variables
  z_var A(vfac["objA"], crab::REF_TYPE, 32);
  z_var B(vfac["objB"], crab::REF_TYPE, 32);
  z_var a(vfac["a"], crab::INT_TYPE, 32);
  z_var b(vfac["b"], crab::INT_TYPE, 32);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var w(vfac["w"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var r1(vfac["r1"], crab::REF_TYPE, 32);
  z_var r2(vfac["r2"], crab::REF_TYPE, 32);
  z_var p1(vfac["p1"], crab::REF_TYPE, 32);
  z_var p2(vfac["p2"], crab::REF_TYPE, 32);
  z_var rgn1(vfac["V_a"], crab::REG_INT_TYPE, 32);
  z_var rgn2(vfac["V_b"], crab::REG_INT_TYPE, 32);
  z_var rgn3(vfac["V_c"], crab::REG_INT_TYPE, 32);

  // Make prestate for dom
  variable_vector_t obj1;
  init += (rgn1 >= z_number(0));
  init += (rgn1 <= z_number(1));
  obj1.push_back(rgn1);
  init += (rgn2 >= z_number(1));
  init += (rgn2 <= z_number(2));
  obj1.push_back(rgn2);
  init.assign(rgn3, rgn2);
  obj1.push_back(rgn3);
  simulate_obj_map_update_per_object(init, obj1);
  // mem: { 0 <= V_a <= 1 && V_b = V_c && 1 <= V_b <= 2 }, cache is empty

  // Defining size
  z_var_or_cst_t zero32(z_number(0), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t one32(z_number(1), crab::variable_type(crab::INT_TYPE, 32));  
  z_var_or_cst_t size4(z_number(4), crab::variable_type(crab::INT_TYPE, 32));
  // Create allocation sites
  crab::tag_manager as_man;
  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");

  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");

  // adding control flow
  entry.add_succ(bb1);
  bb1.add_succ(bb2);
  bb2.add_succ(ret);
  // adding statements

  // Intialization of memory regions
  // Assume regions have been initialized
//   entry.region_init(rgn1);
//   entry.region_init(rgn2);
//   entry.region_init(rgn3);
  // Create object
  entry.make_ref(A, rgn1, size4, as_man.mk_tag());
  entry.make_ref(B, rgn1, size4, as_man.mk_tag());
  bb1.gep_ref(r1, rgn2, A, rgn1, z_number(4)); // V_b, r1 = gep_ref(V_a, &A, 4)
  bb1.gep_ref(p1, rgn3, r1, rgn2, z_number(4)); // V_c, p1 = gep_ref(V_b, r1, 4)

  bb1.load_from_ref(x, r1, rgn2); // x = load_ref(V_b, r1)
  bb1.load_from_ref(y, p1, rgn3); // y = load_ref(V_c, p1)

  bb2.gep_ref(r2, rgn2, B, rgn1, z_number(4)); // V_b,r2 = gep_ref(V_a, &B, 4)
  bb2.gep_ref(p2, rgn3, r2, rgn2, z_number(4)); // V_c,p2 = gep_ref(V_b, r2, 4)
  bb2.load_from_ref(z, r2, rgn2); // z = load_ref(V_b, r2) // cache miss, update cache
  bb2.load_from_ref(w, p2, rgn3); // w = load_ref(V_c, p2)
  return cfg;
}

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  z_mru_rgn_zones_t init;
  variable_factory_t vfac;
  z_cfg_t *cfg = prog(init, vfac);
  cfg->simplify();
  crab::outs() << *cfg << "\n";

  run(cfg, cfg->entry(), init, false, 2, 2, 20, stats_enabled);

  return 0;
}