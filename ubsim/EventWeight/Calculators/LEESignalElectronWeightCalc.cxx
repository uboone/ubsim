/**
 * \class evwgh::LEESignalElectronWeightCalc
 * \brief eLEE Signal reweighting weightcalc
 * \author W. Ketchum <wketchum@fnal.gov> (24 Feb 2019)
 * 
 * Shamelessly copied from ReinteractionWeightCalc by A. Mastbaum.
 * Shamelessly incorporating the work of the MiniBooNE LEE reweight 
 * team.
 *
 * Reweight nu_e events to reproduce an LEE signal like model
 */

#include <map>
#include <string>

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Provenance/ModuleContext.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/EventWeight/Base/WeightCalcCreator.h"
#include "larsim/EventWeight/Base/WeightCalc.h"

namespace evwgh {

class LEESignalElectronWeightCalc : public WeightCalc {
public:
  LEESignalElectronWeightCalc() {}

  void Configure(fhicl::ParameterSet const& p, CLHEP::HepRandomEngine& );

  std::vector<std::vector<double> > GetWeight(art::Event& e);

private:
  art::InputTag fMCTruthProducer; //!< Label for MCTruth producer
  std::vector<float>  fEnergyBinEdges; //!< list of the bin edges to use in energy
  std::vector<double> fWeights; //!< List of the weights to apply per bin.

  DECLARE_WEIGHTCALC(LEESignalElectronWeightCalc)
};


void LEESignalElectronWeightCalc::Configure(fhicl::ParameterSet const& p,
                                        CLHEP::HepRandomEngine&)
{
  // Get configuration for this function
  fhicl::ParameterSet const& pset = p.get<fhicl::ParameterSet>(GetName());

  fMCTruthProducer = pset.get<art::InputTag>("MCTruthProducer");
  fEnergyBinEdges = pset.get< std::vector<float> >("EnergyBinEdges");
  fWeights = pset.get< std::vector<double> >("Weights");

  if(fEnergyBinEdges.size()+1!=fWeights.size())
    throw cet::exception("LEESignalElectronWeightCalc")
      << "BinEdges size is " << fEnergyBinEdges.size()
      << " so weights size should be " << fEnergyBinEdges.size()+1
      << " but it's " << fWeights.size() << ".";

  for(size_t i_ebin=0; i_ebin<fEnergyBinEdges.size()-1; ++i_ebin)
    if(fEnergyBinEdges[i_ebin]>fEnergyBinEdges[i_ebin+1])
      throw cet::exception("LEESignalElectronWeightCalc")
	<< "Energy bin edges are not monotonically increasing!";
}
  
std::vector<std::vector<double> >
LEESignalElectronWeightCalc::GetWeight(art::Event& e) {

  // Get MCTruths in the event ...
  art::Handle<std::vector<simb::MCTruth> > truthHandle;
  e.getByLabel(fMCTruthProducer, truthHandle);
  std::vector<simb::MCTruth> const& truthVec(*truthHandle);

  // Initialize the vector of event weights
  std::vector<std::vector<double> > weight(truthVec.size());

  // Loop over the MCTruth objects...
  for (size_t itruth=0; itruth<truthVec.size(); ++itruth){

    // Initialize weight vector to zero for this MCTruth
    weight[itruth].resize(1, 0.0);

    auto const& mctruth = truthVec[itruth];
    if(mctruth.Origin()!=simb::Origin_t::kBeamNeutrino) continue;

    //get the neutrino MCParticle;
    auto const& nu = mctruth.GetNeutrino().Nu();
    
    //if not electron or anti-electron neutrino, get out
    if(std::abs(nu.PdgCode())!=12) continue;
    
    double nu_energy = nu.E();
    size_t nu_ebin=0;
    for(size_t i_ebin=0; i_ebin<fEnergyBinEdges.size(); ++i_ebin){
      if(nu_energy<fEnergyBinEdges[i_ebin]) break;
      ++nu_ebin;
    }
    weight[itruth][0] = fWeights[nu_ebin];
    
  }

  return weight;
}

REGISTER_WEIGHTCALC(LEESignalElectronWeightCalc)

}  // namespace evwgh
