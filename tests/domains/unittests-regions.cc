#include "../common.hpp"
#include "../program_options.hpp"

using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace ikos;

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  
  variable_factory_t vfac;
  crab::tag_manager as_man;
  
  { // join
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var ref(vfac["ref"], crab::REF_TYPE, 32);    
    z_var rgn1(vfac["region_0"], crab::REG_INT_TYPE, 32);
    z_var_or_cst_t n34_32(z_number(34), crab::variable_type(crab::INT_TYPE, 32));
    
    z_rgn_int_t inv1, inv2, inv3;
    inv1 += (x >= z_number(5));
    inv1 += (y >= z_number(10));
    inv1.region_init(rgn1);
    
    inv2 += (x >= z_number(5));
    inv2.region_init(rgn1);
    inv2.ref_make(ref, rgn1, as_man.mk_tag());
    inv2.ref_store(ref, rgn1, n34_32);


    crab::outs() << "Join of\n\t" << inv1  << "\nand\n\t" << inv2 << " =\n\t";
    inv3 = inv1 | inv2;
    crab::outs() << inv3 << "\n";
  }
  
  { // join
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var ref(vfac["ref"], crab::REF_TYPE, 32);    
    z_var rgn1(vfac["region_0"], crab::REG_INT_TYPE, 32);

    z_var_or_cst_t n34_32(z_number(34), crab::variable_type(crab::INT_TYPE, 32));
    z_rgn_int_t inv1, inv2;
    inv1 += (x >= z_number(5));
    inv1 += (y >= z_number(10));
    inv1.region_init(rgn1);
    
    inv2 += (x >= z_number(5));
    inv2.region_init(rgn1);
    inv2.ref_make(ref, rgn1, as_man.mk_tag());
    inv2.ref_store(ref, rgn1, n34_32);


    crab::outs() << "Join of\n\t" << inv1  << "\nand\n\t" << inv2 << " =\n\t";
    inv1 |= inv2;
    crab::outs() << inv1 << "\n";
  }


  { // select_ref
    
    { // both select values are possible
      z_cfg_t cfg("entry","exit");
      z_basic_block_t &entry = cfg.insert("entry");
      z_basic_block_t &exit = cfg.insert("exit");
      entry >> exit;
      z_var ref1(vfac["ref1"], crab::REF_TYPE, 32);
      z_var ref2(vfac["ref2"], crab::REF_TYPE, 32);
      z_var ref3(vfac["ref3"], crab::REF_TYPE, 32);
      z_var rgn1(vfac["rgn1"], crab::REG_INT_TYPE, 32);
      z_var rgn2(vfac["rgn2"], crab::REG_INT_TYPE, 32);
      z_var rgn3(vfac["rgn3"], crab::REG_INT_TYPE, 32);
      z_var cond(vfac["cond"], crab::BOOL_TYPE, 1);    
      z_rgn_int_t init;      
      entry.region_init(rgn1);
      entry.region_init(rgn2);
      entry.region_init(rgn3);      
      entry.assume_ref(z_ref_cst_t::mk_gt_null(ref2));
      entry.assume_ref(z_ref_cst_t::mk_gt_null(ref3));
      entry.havoc(cond);
      entry.select_ref(ref1, rgn1, cond, ref2, rgn2, ref3, rgn3);
      exit.assert_ref(z_ref_cst_t::mk_not_null(ref1));
      crab::outs() << "Unit test 1 for select_ref\n";
      crab::outs() << cfg;
      run_and_check(&cfg, cfg.entry(), init, false, 2, 2, 20, stats_enabled);
    }
    { // both select values are possible: the only expected warning
      z_cfg_t cfg("entry","exit");
      z_basic_block_t &entry = cfg.insert("entry");
      z_basic_block_t &exit = cfg.insert("exit");
      entry >> exit;
      z_var ref1(vfac["ref1"], crab::REF_TYPE, 32);
      z_var ref2(vfac["ref2"], crab::REF_TYPE, 32);
      z_var ref3(vfac["ref3"], crab::REF_TYPE, 32);
      z_var rgn1(vfac["rgn1"], crab::REG_INT_TYPE, 32);
      z_var rgn2(vfac["rgn2"], crab::REG_INT_TYPE, 32);
      z_var rgn3(vfac["rgn3"], crab::REG_INT_TYPE, 32);
      z_var cond(vfac["cond"], crab::BOOL_TYPE, 1);    
      z_rgn_bool_int_t init;
      entry.region_init(rgn1);
      entry.region_init(rgn2);
      entry.region_init(rgn3);
      entry.assume_ref(z_ref_cst_t::mk_gt_null(ref2));
      //entry.bool_assign(cond, z_lin_cst_t::get_true());
      entry.select_ref(ref1, rgn1, cond, ref2, rgn2, ref3, rgn3);
      exit.assert_ref(z_ref_cst_t::mk_not_null(ref1));
      crab::outs() << "Unit test 2 for select_ref\n";            
      crab::outs() << cfg;
      run_and_check(&cfg, cfg.entry(), init, false, 2, 2, 20, stats_enabled);
    }
    
    { // both select values are possible and one operand is null
      z_cfg_t cfg("entry","exit");
      z_basic_block_t &entry = cfg.insert("entry");
      z_basic_block_t &exit = cfg.insert("exit");
      entry >> exit;
      z_var ref1(vfac["ref1"], crab::REF_TYPE, 32);
      z_var ref2(vfac["ref2"], crab::REF_TYPE, 32);
      z_var ref3(vfac["ref3"], crab::REF_TYPE, 32);
      z_var rgn1(vfac["rgn1"], crab::REG_INT_TYPE, 32);
      z_var rgn2(vfac["rgn2"], crab::REG_INT_TYPE, 32);
      z_var rgn3(vfac["rgn3"], crab::REG_INT_TYPE, 32);
      z_var cond(vfac["cond"], crab::BOOL_TYPE, 1);    
      z_rgn_int_t init;      
      entry.region_init(rgn1);
      entry.region_init(rgn2);
      entry.region_init(rgn3);      
      entry.assume_ref(z_ref_cst_t::mk_null(ref2));
      entry.havoc(cond);
      entry.select_ref_null_true_value(ref1, rgn1, cond, ref2, rgn2);
      exit.assert_ref(z_ref_cst_t::mk_null(ref1));
      crab::outs() << "Unit test 3 for select_ref\n";      
      crab::outs() << cfg;
      run_and_check(&cfg, cfg.entry(), init, false, 2, 2, 20, stats_enabled);
    }
    { // only select true value
      z_cfg_t cfg("entry","exit");
      z_basic_block_t &entry = cfg.insert("entry");
      z_basic_block_t &exit = cfg.insert("exit");
      entry >> exit;
      z_var ref1(vfac["ref1"], crab::REF_TYPE, 32);
      z_var ref2(vfac["ref2"], crab::REF_TYPE, 32);
      z_var ref3(vfac["ref3"], crab::REF_TYPE, 32);
      z_var rgn1(vfac["rgn1"], crab::REG_INT_TYPE, 32);
      z_var rgn2(vfac["rgn2"], crab::REG_INT_TYPE, 32);
      z_var rgn3(vfac["rgn3"], crab::REG_INT_TYPE, 32);
      z_var cond(vfac["cond"], crab::BOOL_TYPE, 1);    
      z_rgn_bool_int_t init;
      entry.region_init(rgn1);
      entry.region_init(rgn2);
      entry.region_init(rgn3);
      entry.assume_ref(z_ref_cst_t::mk_gt_null(ref2));
      entry.bool_assign(cond, z_lin_cst_t::get_true());
      entry.select_ref(ref1, rgn1, cond, ref2, rgn2, ref3, rgn3);
      exit.assert_ref(z_ref_cst_t::mk_not_null(ref1));
      crab::outs() << "Unit test 4 for select_ref\n";
      crab::outs() << cfg;
      run_and_check(&cfg, cfg.entry(), init, false, 2, 2, 20, stats_enabled);
    }
    { // only select false value
      z_cfg_t cfg("entry","exit");
      z_basic_block_t &entry = cfg.insert("entry");
      z_basic_block_t &exit = cfg.insert("exit");
      entry >> exit;
      z_var ref1(vfac["ref1"], crab::REF_TYPE, 32);
      z_var ref2(vfac["ref2"], crab::REF_TYPE, 32);
      z_var ref3(vfac["ref3"], crab::REF_TYPE, 32);
      z_var rgn1(vfac["rgn1"], crab::REG_INT_TYPE, 32);
      z_var rgn2(vfac["rgn2"], crab::REG_INT_TYPE, 32);
      z_var rgn3(vfac["rgn3"], crab::REG_INT_TYPE, 32);
      z_var cond(vfac["cond"], crab::BOOL_TYPE, 1);    
      z_rgn_bool_int_t init;
      entry.assume_ref(z_ref_cst_t::mk_gt_null(ref3));
      entry.bool_assign(cond, z_lin_cst_t::get_false());
      entry.select_ref(ref1, rgn1, cond, ref2, rgn2, ref3, rgn3);
      exit.assert_ref(z_ref_cst_t::mk_not_null(ref1));
      crab::outs() << "Unit test 5 for select_ref\n";      
      crab::outs() << cfg;
      run_and_check(&cfg, cfg.entry(), init, false, 2, 2, 20, stats_enabled);
    }
  }
  return 0;
}
