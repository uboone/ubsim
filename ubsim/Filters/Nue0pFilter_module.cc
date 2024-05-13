////////////////////////////////////////////////////////////////////////
// Class:       Nue0pFilter
// Plugin Type: filter (art v3_05_01)
// File:        Nue0pFilter_module.cc
//
// Generated at Sun Oct 1 00:00:00 2023 by David Caratelli using cetskelgen
// from cetlib version ????.
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

class Nue0pFilter;


class Nue0pFilter : public art::EDFilter {
public:
  explicit Nue0pFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Nue0pFilter(Nue0pFilter const&) = delete;
  Nue0pFilter(Nue0pFilter&&) = delete;
  Nue0pFilter& operator=(Nue0pFilter const&) = delete;
  Nue0pFilter& operator=(Nue0pFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  art::InputTag fMCTproducer;
  float fProtonKEThreshold;
  bool fProtonBool; // 0 -> 0p, 1 -> Np

  // Declare member data here.

};


Nue0pFilter::Nue0pFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fMCTproducer = p.get<art::InputTag>("MCTproducer");
  // threshold of 0 should allow any event with protons regardless of KE
  // threshold of 0 should return true if no true protons are present
  fProtonKEThreshold = p.get<float>("ProtonKEThreshold",0.);
  fProtonBool = p.get<bool>("ProtonBool",0);

}

bool Nue0pFilter::filter(art::Event& e)
{
  // Implementation of required member function here.

  // load MCTruth
  auto const &mct_h = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTproducer);

  auto mct = mct_h->at(0);

  // check that neutrino is nue
  auto neutrino = mct.GetNeutrino();
  //auto lep_e = neutrino.Lepton().E();
  //auto lep_pdg = neutrino.Lepton().PdgCode();
  auto nu_pdg = neutrino.Nu().PdgCode();

  if (abs(nu_pdg) != 12) return false;

  // check hadronic activity looking for highest energy proton
  float LeadingProtonKE = -1; // initialize leading proton's KE [GeV]

  size_t npart = mct.NParticles();

  for (size_t i = 0; i < npart; i++) {
    
    auto const &part = mct.GetParticle(i);
    
    if (part.PdgCode() == 2212){
      
      // GeV
      float protonKE = part.Momentum(0).E() - 0.938;
      if (protonKE > LeadingProtonKE) {LeadingProtonKE = protonKE;}

    }//if proton
    
  }// for all MCParticles

  if ( (!fProtonBool)  && (LeadingProtonKE > fProtonKEThreshold) ) return false;
  if ( (fProtonBool) && (LeadingProtonKE < fProtonKEThreshold) ) return false;

  std::cout << "[Nue0pFilter] pass with leading proton energy of " << LeadingProtonKE << std::endl;

  return true;
}

void Nue0pFilter::beginJob()
{
  // Implementation of optional member function here.
}

void Nue0pFilter::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(Nue0pFilter)
