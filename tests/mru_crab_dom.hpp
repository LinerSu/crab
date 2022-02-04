#include "crab_dom.hpp"
#include <crab/domains/mru_region_domain.hpp>

namespace crab {
namespace mru_domain_impl {
  using namespace crab::domain_impl;
  using z_rgn_zones_params_t = TestRegionParams<z_soct_domain_t>;
  using z_mru_rgn_zones_t = mru_region_domain<z_rgn_zones_params_t>;
}
}