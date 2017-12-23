#ifndef CLASSICHITORDERER_CXX
#define CLASSICHITORDERER_CXX

#include "ClassicHitOrderer.h"


namespace cosmictag {

  static ClassicHitOrdererFactory __global_ClassicHitOrdererFactory__;

  ClassicHitOrderer::ClassicHitOrderer(const std::string name)
    : BaseHitOrdererAlgo(name)
  {}

  void ClassicHitOrderer::_Configure_(const Config_t &pset)
  {
    _global_qe = pset.get<double>("GlobalQE");
  }
  
  bool ClassicHitOrderer::OrderHits(const SimpleCluster& cluster) const
  {
 
    std::cout << "Test" << std::endl;
    return true;
  }
}
#endif
