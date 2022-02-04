#include "../../program_options.hpp"
#include "../../common.hpp"
#include <crab/domains/mru_region_domain.hpp>

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

namespace {
  using z_rgn_zones_params_t = TestRegionParams<z_soct_domain_t>;
  using z_mru_rgn_zones_t = mru_region_domain<z_rgn_zones_params_t>;
}

int main(int argc, char** argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc,argv,stats_enabled)) {
      return 0;
  }
  variable_factory_t vfac;
  
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var w(vfac["w"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var ref(vfac["ref"], crab::REF_TYPE, 32);
  z_var ref2(vfac["ref2"], crab::REF_TYPE, 32);
  z_var rgn1(vfac["V_1"], crab::REG_INT_TYPE, 32);
  z_var rgn2(vfac["V_2"], crab::REG_INT_TYPE, 32);
  z_var cache_addr(vfac["C_base"], crab::INT_TYPE, 32);

  crab::outs() << "=== Join === \n";

  // join with cache empty
  crab::outs() << "=== 1. join with empty caches === \n";
  {
    z_mru_rgn_zones_t left;
    z_mru_rgn_zones_t right;
    left += (x >= z_number(0));
    left += (x <= z_number(1));
    left += (rgn1 >= z_number(1));
    left += (rgn1 <= z_number(2));
    left += (rgn1 - rgn2 < z_number(0));
    // left: 
    //  mem: { 1 <= V_1 <= 2; V_1 < V_2; 0 <= x <= 1 }
    //  cache: empty

    right += (rgn1 >= z_number(0));
    right += (rgn2 <= y);
    right += (rgn1 - rgn2 <= z_number(0));
    right += (y <= z_number(5));
    // right: 
    //  mem: { 0 <= V_1; V_1 <= V_2; V_2 <= y; y <= 5}
    //  cache: empty

    z_mru_rgn_zones_t l_join_r = left | right;
    crab::outs() << left << " | \n" << right << " = \n" << l_join_r << "\n";
    // Mem= {-V_2 <= 0; -V_1 <= 0; V_1 <= 5; -V_2+V_1 <= 0},
    // ( CacheLine= {},
    // Regs= {},
    // Addrs= {},
    // C_used= false,
    // C_dirty= false )
  }
  // join with one cache is not empty
  crab::outs() << "=== 2. join with only one cache is dirty === \n";
  {
    z_mru_rgn_zones_t left;
    z_mru_rgn_zones_t right;

    left += (x >= z_number(0));
    left += (x <= z_number(1));
    left += (rgn1 >= z_number(1));
    left += (rgn1 <= z_number(2));
    left += (rgn1 - rgn2 < z_number(0));
    left += (rgn2 <= z_number(4));
    left.cache_lines_assign(rgn1, z_number(3));
    left.cache_lines_assign(rgn2, z_number(6));
    left.addrs_assign(cache_addr, ref);
    left.m_used = boolean_value::get_true();
    left.m_dirty = boolean_value::get_true();
    // left: 
    //  mem: { 1 <= V_1 <= 2; V_1 < V_2; V_2 <= 4; 0 <= x <= 1; }
    //  cache: ( 
    //   CacheLine: { V_1 = 3; V_2 = 6 },
    //   Regs: top,
    //   Addrs: { C_base == ref.base },
    //   C_dirty = true,
    //   C_used = true
    //  )

    right += (rgn1 >= z_number(0));
    right += (rgn2 <= y);
    right += (rgn1 - rgn2 <= z_number(0));
    right += (y <= z_number(5));
    // right: 
    //  mem: { 0 <= V_1; V_1 <= V_2; V_2 <= y; y <= 5}
    //  cache: empty

    z_mru_rgn_zones_t l_join_r = left | right;
    crab::outs() << left << " | " << right << " = " << l_join_r << "\n";
    // Mem= {-V_1 <= 0; V_1 <= 5; -V_2 <= 0; V_2 <= 6; -V_1+V_2 <= 5; V_1-V_2 <= 0; V_1+V_2 <= 10},
    // ( CacheLine= {},
    // Regs= {},
    // Addrs= {},
    // C_used= false,
    // C_dirty= false )
  }
  // join with caches referring the same object
  crab::outs() << "=== 3. join with two caches refer to same object === \n";
  {
    z_mru_rgn_zones_t left;
    z_mru_rgn_zones_t right;

    left += (x >= z_number(0));
    left += (x <= z_number(1));
    left += (rgn1 >= z_number(1));
    left += (rgn1 <= z_number(2));
    left += (rgn1 - rgn2 < z_number(0));
    left += (rgn2 <= z_number(4));
    left.cache_lines_assign(rgn1, z_number(3));
    left.cache_lines_assign(rgn2, z_number(6));
    left.addrs_assign(cache_addr, ref);
    left.m_used = boolean_value::get_true();
    left.m_dirty = boolean_value::get_true();
    // left: 
    //  mem: { 1 <= V_1 <= 2; V_1 < V_2; V_2 <= 4; 0 <= x <= 1; }
    //  cache: ( 
    //   CacheLine: { V_1 = 3; V_2 = 6 },
    //   Regs: top,
    //   Addrs: { C_base == ref.base },
    //   C_dirty = true,
    //   C_used = true
    //  )

    right += (rgn1 >= z_number(0));
    right += (rgn2 <= y);
    right += (rgn1 - rgn2 <= z_number(0));
    right += (y <= z_number(5));
    right.cache_lines_assign(rgn1, z_number(3));
    right.cache_lines_assign(rgn2, z_number(6));
    right.add_cons_into_cache_lines(rgn1 <= z_number(4));
    right.addrs_assign(cache_addr, ref);
    right.m_used = boolean_value::get_true();
    right.m_dirty = boolean_value::get_true();
    right.regs_assign(x, rgn1);
    // right: 
    //  mem: { 0 <= V_1; V_1 <= V_2; V_2 <= y; y <= 5}
    //  cache: ( 
    //   CacheLine: { V_1 = 3; V_2 = 6 },
    //   Regs: { C_1 == y },
    //   Addrs: { C_base == ref.base },
    //   C_dirty = true,
    //   C_used = true
    //  )

    z_mru_rgn_zones_t l_join_r = left | right;
    crab::outs() << left << " | \n" << right << " = \n" << l_join_r << "\n";
    // res:
    // Mem= {-V_2 <= 0; V_2 <= 5; -V_1 <= 0; V_1 <= 5; -V_2+V_1 <= 0},
    // ( CacheLine= {V_2 = 6; V_1 = 3},
    // Regs= {},
    // Addrs= {C_base -> $VAR_0},
    // C_used= true,
    // C_dirty= true )
  }
  // join with caches referring the different object
  crab::outs() << "=== 4. join with two caches refer to different object === \n";
  {
    z_mru_rgn_zones_t left;
    z_mru_rgn_zones_t right;

    left += (x >= z_number(0));
    left += (x <= z_number(1));
    left += (rgn1 >= z_number(1));
    left += (rgn1 <= z_number(2));
    left += (rgn1 - rgn2 < z_number(0));
    left += (rgn2 <= z_number(4));
    left.cache_lines_assign(rgn1, z_number(3));
    left.cache_lines_assign(rgn2, z_number(6));
    left.addrs_assign(cache_addr, ref);
    left.m_used = boolean_value::get_true();
    left.m_dirty = boolean_value::get_true();
    // left: 
    //  mem: { 1 <= V_1 <= 2; V_1 < V_2; V_2 <= 4; 0 <= x <= 1; }
    //  cache: ( 
    //   CacheLine: { V_1 = 3; V_2 = 6 },
    //   Regs: top,
    //   Addrs: { C_base == ref.base },
    //   C_dirty = true,
    //   C_used = true
    //  )

    right += (rgn1 >= z_number(0));
    right += (rgn2 <= y);
    right += (rgn1 - rgn2 <= z_number(0));
    right += (y <= z_number(5));
    right.cache_lines_assign(rgn1, z_number(3));
    right.cache_lines_assign(rgn2, z_number(6));
    right.add_cons_into_cache_lines(rgn1 <= z_number(4));
    right.addrs_assign(cache_addr, ref2);
    right.m_used = boolean_value::get_true();
    right.m_dirty = boolean_value::get_true();
    right.regs_assign(y, rgn1);
    // right: 
    //  mem: { 0 <= V_1; V_1 <= V_2; V_2 <= y; y <= 5}
    //  cache: ( 
    //   CacheLine: { V_1 = 3; V_2 = 6 },
    //   Regs: { C_1 == y },
    //   Addrs: { C_base == ref2.base },
    //   C_dirty = true,
    //   C_used = true
    //  )

    z_mru_rgn_zones_t l_join_r = left | right;
    crab::outs() << left << " | \n" << right << " = \n" << l_join_r << "\n";
  }

  crab::outs() << "=== Project === \n";
  {
    z_mru_rgn_zones_t dom;
    dom += (rgn1 >= z_number(0));
    dom += (rgn2 <= y);
    dom += (rgn1 - rgn2 <= z_number(0));
    dom += (y <= z_number(5));
    dom.cache_lines_assign(rgn1, z_number(3));
    dom.cache_lines_assign(rgn2, z_number(6));
    dom.add_cons_into_cache_lines(rgn1 <= z_number(4));
    dom.addrs_assign(cache_addr, ref2);
    dom.m_used = boolean_value::get_true();
    dom.m_dirty = boolean_value::get_true();
    dom.regs_assign(y, rgn1);
    // dom: 
    //  mem: { 0 <= V_1; V_1 <= V_2; V_2 <= y; y <= 5}
    //  cache: ( 
    //   CacheLine: { V_1 = 3; V_2 = 6 },
    //   Regs: { C_1 == y },
    //   Addrs: { C_base == ref2.base },
    //   C_dirty = true,
    //   C_used = true
    //  )
    std::vector<z_var> vs;
    vs.push_back(y);
    vs.push_back(rgn2);
    crab::outs() << "Projecting: ";
    for (unsigned i = 0; i < vs.size(); ++i)
      crab::outs() << vs[i] << " ";
    crab::outs() << "\n";
    crab::outs() << "Before project \n" << dom;
    dom.project(vs);
    crab::outs() << "After project \n" << dom;
  }

  crab::outs() << "==== Meet ====\n";
  {
    z_mru_rgn_zones_t left;
    z_mru_rgn_zones_t right;

    left += (x >= z_number(0));
    left += (x <= z_number(1));
    left += (rgn1 >= z_number(1));
    left += (rgn1 <= z_number(2));
    left += (rgn1 - rgn2 < z_number(0));
    left += (rgn2 <= z_number(4));
    left.cache_lines_assign(rgn1, z_number(3));
    left.cache_lines_assign(rgn2, z_number(6));
    left.addrs_assign(cache_addr, ref);
    left.m_used = boolean_value::get_true();
    left.m_dirty = boolean_value::get_true();
    // left: 
    //  mem: { 1 <= V_1 <= 2; V_1 < V_2; V_2 <= 4; 0 <= x <= 1; }
    //  cache: ( 
    //   CacheLine: { V_1 = 3; V_2 = 6 },
    //   Regs: top,
    //   Addrs: { C_base == ref.base },
    //   C_dirty = true,
    //   C_used = true
    //  )

    right += (rgn1 >= z_number(0));
    right += (rgn2 <= y);
    right += (rgn1 - rgn2 <= z_number(0));
    right += (y <= z_number(5));
    right.cache_lines_assign(rgn1, z_number(3));
    right.cache_lines_assign(rgn2, z_number(6));
    right.add_cons_into_cache_lines(rgn1 <= z_number(4));
    right.addrs_assign(cache_addr, ref2);
    right.m_used = boolean_value::get_true();
    right.m_dirty = boolean_value::get_true();
    right.regs_assign(y, rgn1);
    // right: 
    //  mem: { 0 <= V_1; V_1 <= V_2; V_2 <= y; y <= 5}
    //  cache: ( 
    //   CacheLine: { V_1 = 3; V_2 = 6 },
    //   Regs: { C_1 == y },
    //   Addrs: { C_base == ref2.base },
    //   C_dirty = true,
    //   C_used = true
    //  )

    z_mru_rgn_zones_t l_meet_r = left & right;
    crab::outs() << left << " & \n" << right << " = \n" << l_meet_r << "\n";
  }

  return 0;
}
