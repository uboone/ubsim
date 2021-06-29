////////////////////////////////////////////////////////////////////////
// Class:       NeutrinoTruthFilter
// Plugin Type: filter (art v3_01_01)
// File:        NeutrinoTruthFilter_module.cc
//
// Generated at Mon Feb 25 09:55:32 2019 by Wesley Ketchum using cetskelgen
// from cetlib version v3_05_01.
//
//
// 25 Feb 2019
// This is a super simple filter to just select neutrinos in an energy range.
// It can/could be extended to do a lot more!
// 
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

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCTruth.h"

namespace sim {
  class NeutrinoTruthFilter;
}


class sim::NeutrinoTruthFilter : public art::EDFilter {
public:
  explicit NeutrinoTruthFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NeutrinoTruthFilter(NeutrinoTruthFilter const&) = delete;
  NeutrinoTruthFilter(NeutrinoTruthFilter&&) = delete;
  NeutrinoTruthFilter& operator=(NeutrinoTruthFilter const&) = delete;
  NeutrinoTruthFilter& operator=(NeutrinoTruthFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:

  // Declare member data here.
  art::InputTag fMCTruthLabel;
  float fEMin;
  float fEMax;

  void reconfigure(fhicl::ParameterSet const& p);

};


sim::NeutrinoTruthFilter::NeutrinoTruthFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{
  reconfigure(p);
}

bool sim::NeutrinoTruthFilter::filter(art::Event& e)
{

  // Get MCTruths in the event ...
  art::Handle<std::vector<simb::MCTruth> > truthHandle;
  e.getByLabel(fMCTruthLabel, truthHandle);
  std::vector<simb::MCTruth> const& truthVec(*truthHandle);

  bool keep_event=false;

  // Loop over the MCTruth objects...
  for (auto const& mctruth : truthVec){

    if(mctruth.Origin()!=simb::Origin_t::kBeamNeutrino) continue;

    //get the neutrino MCParticle;
    auto const& nu = mctruth.GetNeutrino().Nu();

    //if we sit in our desired energy, we want to keep it!
    if(nu.E() > fEMin && nu.E() < fEMax) keep_event=true;

  }

  return keep_event;

}

void sim::NeutrinoTruthFilter::reconfigure(fhicl::ParameterSet const& p)
{
  fMCTruthLabel = p.get<art::InputTag>("MCTruthLabel");
  fEMin = p.get<float>("EMin");
  fEMax = p.get<float>("EMax");
}

DEFINE_ART_MODULE(sim::NeutrinoTruthFilter)
