#include "../../program_options.hpp"
#include "../../common.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t* cfg1(variable_factory_t &vfac)  {

  /*
   int v1,v2,v3;
   int *i = &v1;
   int *x = &v2;
   int *y = &v3;

   *i := 0;
   *x := 1;
   *y := 0;

    while(*i <= 99) {
     *x = *x + *y;
     *y = *y + 1;
     *i = *i + 1;
    }

    assert(*x >= *y);
   */
  
  // === Define program variables
  z_var i(vfac["i"], crab::REF_TYPE);
  z_var x(vfac["x"], crab::REF_TYPE);
  z_var y(vfac["y"], crab::REF_TYPE);
  z_var deref_i(vfac["*i"], crab::INT_TYPE, 32);
  z_var deref_x(vfac["*x"], crab::INT_TYPE, 32);
  z_var deref_y(vfac["*y"], crab::INT_TYPE, 32);
  // === Define memory regions
  auto mem1 = crab::memory_region::make_int_memory_region(0, 32);
  auto mem2 = crab::memory_region::make_int_memory_region(1, 32);
  auto mem3 = crab::memory_region::make_int_memory_region(2, 32);  
  // === Create empty CFG
  z_cfg_t* cfg = new z_cfg_t("entry","ret",REF);
  // === Adding CFG blocks
  z_basic_block_t& entry = cfg->insert("entry");
  z_basic_block_t& bb1   = cfg->insert("bb1");
  z_basic_block_t& bb1_t = cfg->insert("bb1_t");
  z_basic_block_t& bb1_f = cfg->insert("bb1_f");
  z_basic_block_t& bb2   = cfg->insert("bb2");
  z_basic_block_t& bb3   = cfg->insert("bb3");  
  z_basic_block_t& ret   = cfg->insert("ret");
  // === Adding CFG edges
  entry.add_succ(bb1);
  bb1.add_succ(bb1_t);
  bb1.add_succ(bb1_f);
  bb1_t.add_succ(bb2);
  bb2.add_succ(bb1);
  bb1_f.add_succ(bb3);
  bb3.add_succ(ret);

  
  // === Adding statements

  // Intialization of memory regions
  entry.region_init(mem1);
  entry.region_init(mem2);
  entry.region_init(mem3);

  //// Create references 
  entry.make_ref(i, mem1);
  entry.make_ref(x, mem2);
  entry.make_ref(y, mem3);
  //// *i := 0;
  entry.store_to_ref(i, mem1, z_number(0));
  //// *x := 1;
  entry.store_to_ref(x, mem2, z_number(1));
  //// *y := 0;
  entry.store_to_ref(y, mem3, z_number(0));
  //// assume(*i <= 99);
  bb1_t.load_from_ref(deref_i, i, mem1);
  bb1_t.assume(deref_i <= 99);
  //// assume(*i >= 100);  
  bb1_f.load_from_ref(deref_i, i, mem1);
  bb1_f.assume(deref_i >= 100);
  //// *x = *x + *y
  bb2.load_from_ref(deref_x, x, mem2);
  bb2.load_from_ref(deref_y, y, mem3);
  bb2.add(deref_x, deref_x, deref_y);
  bb2.store_to_ref(x, mem2, deref_x);
  //// *y = *y + 1
  bb2.load_from_ref(deref_y, y, mem3);
  bb2.add(deref_y, deref_y, 1);
  bb2.store_to_ref(y, mem3, deref_y);  
  //// *i = *i + 1
  bb2.load_from_ref(deref_i, i, mem1);
  bb2.add(deref_i, deref_i, 1);
  bb2.store_to_ref(i, mem1, deref_i);  
  //// assert(*x >= *y)
  bb3.load_from_ref(deref_x, x, mem2);
  bb3.load_from_ref(deref_y, y, mem3);  
  ret.assertion(deref_x >= deref_y);
  return cfg;
}


int main(int argc, char** argv) {

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc,argv,stats_enabled)) {
      return 0;
  }

  variable_factory_t vfac;
  
  z_cfg_t *p1 = cfg1(vfac);
  crab::outs() << *p1 << "\n";
  run_and_check<z_ref_sdbm_t>(p1,p1->entry(),false,2,2,20,stats_enabled);
  delete p1;

  return 0;
}