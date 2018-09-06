////////////////////////////////////////////////////////////////////////
// Class:       NuMI_KDAR_Filter
// Plugin Type: filter (art v2_11_03)
// File:        NuMI_KDAR_Filter_module.cc
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
  class NuMI_KDAR_Filter;
}


class sim::NuMI_KDAR_Filter : public art::EDFilter {
public:
  explicit NuMI_KDAR_Filter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuMI_KDAR_Filter(NuMI_KDAR_Filter const &) = delete;
  NuMI_KDAR_Filter(NuMI_KDAR_Filter &&) = delete;
  NuMI_KDAR_Filter & operator = (NuMI_KDAR_Filter const &) = delete;
  NuMI_KDAR_Filter & operator = (NuMI_KDAR_Filter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

private:

  // Declare member data here.

};


sim::NuMI_KDAR_Filter::NuMI_KDAR_Filter(fhicl::ParameterSet const &)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
}

bool sim::NuMI_KDAR_Filter::filter(art::Event &)
{
  // Implementation of required member function here.

  return true;
}

DEFINE_ART_MODULE(sim::NuMI_KDAR_Filter)
