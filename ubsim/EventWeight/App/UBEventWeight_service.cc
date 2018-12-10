////////////////////////////////////////////////////////////////////////
// Class:       UBEventWeight
// Plugin Type: service (art v2_11_03)
// File:        UBEventWeight_service.cc
//
// Generated at Mon Nov 26 16:58:47 2018 by Andrew Mastbaum using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "fhiclcpp/ParameterSet.h"

class UBEventWeight {
public:
  explicit UBEventWeight(fhicl::ParameterSet const & p, art::ActivityRegistry & areg) {}
};

DECLARE_ART_SERVICE(UBEventWeight, LEGACY)
DEFINE_ART_SERVICE(UBEventWeight)

