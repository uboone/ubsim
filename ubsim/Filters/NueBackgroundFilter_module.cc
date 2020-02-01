////////////////////////////////////////////////////////////////////////
// Class:       NueBackgroundFilter
// Plugin Type: filter (art v3_01_02)
// File:        NueBackgroundFilter_module.cc
//
// Generated at Tue Dec 17 14:49:01 2019 by Giuseppe Cerati using cetskelgen
// from cetlib version v3_05_01.
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
#include "lardataobj/MCBase/MCShower.h"

class NueBackgroundFilter;


class NueBackgroundFilter : public art::EDFilter {
public:
  explicit NueBackgroundFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NueBackgroundFilter(NueBackgroundFilter const&) = delete;
  NueBackgroundFilter(NueBackgroundFilter&&) = delete;
  NueBackgroundFilter& operator=(NueBackgroundFilter const&) = delete;
  NueBackgroundFilter& operator=(NueBackgroundFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:

  // Declare member data here.
  art::InputTag fMCTruthLabel;
  art::InputTag fMCShowerLabel;
  float fMaxMuonEGeV;
  float fMinElecMCShwEMeV;
  float fMinProtonEGeV;
};


NueBackgroundFilter::NueBackgroundFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}
{
  fMCTruthLabel = p.get<art::InputTag>("MCTruthLabel");
  fMCShowerLabel = p.get<art::InputTag>("MCShowerLabel","");
  fMaxMuonEGeV = p.get<float>("MaxMuonEGeV",-1.);
  fMinElecMCShwEMeV = p.get<float>("MinElecMCShwEMeV",-1.);
  fMinProtonEGeV = p.get<float>("MinProtonEGeV",-1.);
}

bool NueBackgroundFilter::filter(art::Event& e)
{
  // Get MCTruths in the event
  art::Handle<std::vector<simb::MCTruth> > truthHandle;
  e.getByLabel(fMCTruthLabel, truthHandle);
  std::vector<simb::MCTruth> const& truthVec(*truthHandle);

  if (truthVec.size()!=1) { throw cet::exception("NueBackgroundFilter: expect MCTruth vector of size 1"); }

  auto mct = truthVec.at(0);
  auto neutrino = mct.GetNeutrino();
  auto lep_e = neutrino.Lepton().E();
  auto lep_pdg = neutrino.Lepton().PdgCode();

  if (fMaxMuonEGeV>=0.) {
    if (abs(lep_pdg)==13 && (lep_e-0.105)>0.02 && lep_e<fMaxMuonEGeV) {
      // if (1) std::cout << "returning true" << std::endl;
      return true;
    }
  }

  if (fMinProtonEGeV>=0.) {
    size_t npart = mct.NParticles();
    for (size_t i = 0; i < npart; i++) {
      auto const &part = mct.GetParticle(i);
      if (part.StatusCode() != 1) continue;
      if (part.PdgCode() == 2212 && part.StatusCode() == 1 && part.Momentum(0).E()>fMinProtonEGeV) return true;
    }
  }

  // Get MCShowers in the event
  if (fMCShowerLabel!="" && fMinElecMCShwEMeV>=0.) {
    art::Handle<std::vector<sim::MCShower> > mcShwowerHandle;
    e.getByLabel(fMCShowerLabel, mcShwowerHandle);
    std::vector<sim::MCShower> const&  mcShwowerVec(*mcShwowerHandle);
    for (auto& mcs : mcShwowerVec) {
      if (std::abs(mcs.PdgCode())==11) {
	// if (1) std::cout << "elec mcs etot=" << mcs.Start().E() << std::endl;
	if (mcs.Start().E()>fMinElecMCShwEMeV) {
	  // if (1) std::cout << "returning true" << std::endl;
	  return true;
	}
      }
    }
  }

  // if (1) std::cout << "returning false" << std::endl;
  return false;
}

DEFINE_ART_MODULE(NueBackgroundFilter)
