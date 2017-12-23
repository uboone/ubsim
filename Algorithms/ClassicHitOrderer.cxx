#ifndef CLASSICHITORDERER_CXX
#define CLASSICHITORDERER_CXX

#include "ClassicHitOrderer.h"


namespace cosmictag {

  static HitOrdererFactory __global_HitOrdererFactory__;

  ClassicHitOrderer::ClassicHitOrderer(const std::string name)
    : BaseFlashHypothesis(name)
  {}

  void ClassicHitOrderer::_Configure_(const Config_t &pset)
  {
    _global_qe = pset.get<double>("GlobalQE");
  }
  
  bool ClassicHitOrderer::OrderHits(const SimpleCluster& cluster) const
  {
 

    return true;
  }
}
#endif
