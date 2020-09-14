////////////////////////////////////////////////////////////////////////
// Class:       EtaFilter
// Plugin Type: filter (art v3_05_01)
// File:        EtaFilter_module.cc
//
// Generated at Sun Sep  6 14:43:36 2020 by David Caratelli using cetskelgen
// from cetlib version v3_10_00.
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

#include "nusimdata/SimulationBase/MCTruth.h"

#include <memory>

class EtaFilter;


class EtaFilter : public art::EDFilter {
public:
  explicit EtaFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  EtaFilter(EtaFilter const&) = delete;
  EtaFilter(EtaFilter&&) = delete;
  EtaFilter& operator=(EtaFilter const&) = delete;
  EtaFilter& operator=(EtaFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  art::InputTag fMCTproducer;

  // Declare member data here.

};


EtaFilter::EtaFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fMCTproducer = p.get<art::InputTag>("MCTproducer");

}

bool EtaFilter::filter(art::Event& e)
{
  // Implementation of required member function here.

  // load MCTruth
  auto const &mct_h = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTproducer);

  auto mct = mct_h->at(0);

  int neta = 0; // initialize eta count

  size_t npart = mct.NParticles();

  for (size_t i = 0; i < npart; i++) {
    
    auto const &part = mct.GetParticle(i);
    
    // for eta does not have to be statuscode==1
    if (part.PdgCode() == 221)
      neta += 1; 

  }// for all MCParticles

  if (neta == 0) return false;

  return true;
}

void EtaFilter::beginJob()
{
  // Implementation of optional member function here.
}

void EtaFilter::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(EtaFilter)
