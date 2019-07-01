////////////////////////////////////////////////////////////////////////
// Class:       LYSimPhotonScaling
// Plugin Type: producer (art v3_01_02)
// File:        LYSimPhotonScaling_module.cc
//
// Generated at Mon Jun 24 11:07:27 2019 by David Caratelli using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/Simulation/SimPhotons.h"

#include "nutools/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandFlat.h"

#include "larcore/Geometry/Geometry.h"

#include "ubevt/Database/LightYieldService.h"
#include "ubevt/Database/LightYieldProvider.h"
#include "ubevt/Database/UbooneLightYieldProvider.h"

#include <memory>

class LYSimPhotonScaling;


class LYSimPhotonScaling : public art::EDProducer {
public:
  explicit LYSimPhotonScaling(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  LYSimPhotonScaling(LYSimPhotonScaling const&) = delete;
  LYSimPhotonScaling(LYSimPhotonScaling&&) = delete;
  LYSimPhotonScaling& operator=(LYSimPhotonScaling const&) = delete;
  LYSimPhotonScaling& operator=(LYSimPhotonScaling&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  CLHEP::HepRandomEngine& fEngine;

  art::InputTag fSimPhotonProducer;


};


LYSimPhotonScaling::LYSimPhotonScaling(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fEngine(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, "HepJamesRandom", "LYscaling", p, "Seed"))
  // More initializers here.
{

  produces<std::vector<sim::SimPhotons> >();
  fSimPhotonProducer = p.get< art::InputTag >("SimPhotonProducer");

}

void LYSimPhotonScaling::produce(art::Event& e)
{

  // load LY provider
  const lariov::LightYieldProvider& ly_provider = art::ServiceHandle<lariov::LightYieldService>()->GetProvider();

  // produce output SimPhotons
  std::unique_ptr< std::vector<sim::SimPhotons> > SimPhoton_v(new std::vector<sim::SimPhotons> );  

  // load input SimPhotons
  auto const& simphoton_h = e.getValidHandle<std::vector<sim::SimPhotons> >(fSimPhotonProducer);



  // loop through simphotons
  for (size_t s=0; s < simphoton_h->size(); s++) {

    auto const& simphoton = simphoton_h->at(s);

    // create new simphoton which is an identical copy of this one
    // then clear it of all the individual photons
    sim::SimPhotons newsimphoton(simphoton);
    newsimphoton.clear();

    // what's the OpChannel?
    auto opchannel = simphoton.OpChannel();

    // what's the LY drop for this OpChannel?
    auto LYscaling = ly_provider.LYScaling(opchannel);

    std::cout << "For OpChannel " << opchannel << " LY scaling is " << LYscaling << std::endl;

    int ntot = 0;
    int nfin = 0;

    // each SimPhoton is a vector of photons. Loop through each
    for (size_t p=0; p < simphoton.size(); p++) {

      ntot += 1;

      // draw random number
      float prob = CLHEP::RandFlat::shoot(&fEngine, 0., 1.);
      
      if (prob < LYscaling) {

	newsimphoton.emplace_back( simphoton[p] );

	nfin += 1;

      }// if random number above LY for this PMT

    }// for all photosn on this PMT

    std::cout << "\t " << ntot << " -> " << nfin << " simphotons simulated" << std::endl;

    SimPhoton_v->emplace_back( newsimphoton );

  }// for all SimPhotons

    
  e.put(std::move(SimPhoton_v));  

}

void LYSimPhotonScaling::beginJob()
{
  // Implementation of optional member function here.
}

void LYSimPhotonScaling::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(LYSimPhotonScaling)
