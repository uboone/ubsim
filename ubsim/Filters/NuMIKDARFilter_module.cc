////////////////////////////////////////////////////////////////////////
// Class:       NuMIKDARFilter
// Plugin Type: filter (art v2_11_03)
// File:        NuMIKDARFilter_module.cc
//
// Generated at Thu Sep  6 14:51:50 2018 by Wesley Ketchum using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

namespace sim {
  class NuMIKDARFilter;
}


class sim::NuMIKDARFilter : public art::EDFilter {
public:
  explicit NuMIKDARFilter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuMIKDARFilter(NuMIKDARFilter const &) = delete;
  NuMIKDARFilter(NuMIKDARFilter &&) = delete;
  NuMIKDARFilter & operator = (NuMIKDARFilter const &) = delete;
  NuMIKDARFilter & operator = (NuMIKDARFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

private:

  // Declare member data here.

};


sim::NuMIKDARFilter::NuMIKDARFilter(fhicl::ParameterSet const & p)
: EDFilter(p)
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
}

bool sim::NuMIKDARFilter::filter(art::Event &)
{
  // Implementation of required member function here.

  return true;
}

DEFINE_ART_MODULE(sim::NuMIKDARFilter)
