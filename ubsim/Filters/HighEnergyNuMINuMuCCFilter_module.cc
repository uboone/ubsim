////////////////////////////////////////////////////////////////////////
// Class:       HighEnergyNuMINuMuCCFilter
// Plugin Type: filter (art v3_01_02)
// File:        HighEnergyNuMINuMuCCFilter_module.cc
//
// Generated at Wed Dec 16 09:12:59 2020 by Christopher Barnes using cetskelgen
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

// nusimdata                                                                                                                                                                                              
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// MCTracks.
#include "lardataobj/MCBase/MCTrack.h"

#include <memory>

namespace lar {
  class HighEnergyNuMINuMuCCFilter;
}


class lar::HighEnergyNuMINuMuCCFilter : public art::EDFilter {
public:
  explicit HighEnergyNuMINuMuCCFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HighEnergyNuMINuMuCCFilter(HighEnergyNuMINuMuCCFilter const&) = delete;
  HighEnergyNuMINuMuCCFilter(HighEnergyNuMINuMuCCFilter&&) = delete;
  HighEnergyNuMINuMuCCFilter& operator=(HighEnergyNuMINuMuCCFilter const&) = delete;
  HighEnergyNuMINuMuCCFilter& operator=(HighEnergyNuMINuMuCCFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:

  // Declare member data here.
 
  // Parent vertex z position.                                                                                                                                                                            
  double parent_fvz;

  // Neutrino PDG.
  int    neutrino_PDG;

  // Interaction Channel
  int    neutrino_interaction_channel;

  // Neutrino Energy
  double neutrino_energy;

  // Neutrino Truth Coordinates.
  double nu_vtx_x_truth;
  double nu_vtx_y_truth;
  double nu_vtx_z_truth;
  
};


lar::HighEnergyNuMINuMuCCFilter::HighEnergyNuMINuMuCCFilter(fhicl::ParameterSet const& p)
  : EDFilter{p} 
{
}

bool lar::HighEnergyNuMINuMuCCFilter::filter(art::Event& e)
{

  art::Handle<std::vector<simb::MCFlux> > parent_h;
  e.getByLabel("generator", parent_h);

  // make sure parent info looks good                                                                                                                                                                     
  if(!parent_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Parent!" << std::endl;
    throw std::exception();
  }

  // Declare a vector of pointers with the 'MCFlux' object.                                                                                                                                             
  std::vector<art::Ptr<simb::MCFlux> > ParentVec;
  art::fill_ptr_vector(ParentVec, parent_h);

  for (auto& parent : ParentVec){

    parent_fvz   = parent->fvz;

  }

  if ( parent_fvz > 15600. ) return false;

  art::Handle<std::vector<simb::MCTruth> > neutrino_h;
  e.getByLabel("generator", neutrino_h);

  // make sure MCTruth info looks good                                                                                                                                                                     
  if(!neutrino_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Neutrino!"<<std::endl;
    throw std::exception();
  }

  std::vector<art::Ptr<simb::MCTruth> >NeutrinoVec;
  art::fill_ptr_vector(NeutrinoVec, neutrino_h);

  // loop through neutrinos..                                                                                                                                                                  
  for (auto& neutrino : NeutrinoVec ) {

    // Unpack the neutrino object to find an MCParticle.                                                                                                                                                  
    const simb::MCNeutrino& truth_neutrino = neutrino->GetNeutrino();
    const simb::MCParticle& truth_particle = truth_neutrino.Nu();

    // Unpack the neutrino PDG code.
    neutrino_PDG = truth_particle.PdgCode();

    // Unpack the interaction channel of the neutrino.
    neutrino_interaction_channel = truth_neutrino.CCNC();

    // Unpack the energy of the neutrino (in MeV).
    neutrino_energy = truth_particle.E( 0 ) * 1000.;

    // Unpack the coordinates for the vertex.
    nu_vtx_x_truth                         = truth_particle.Vx(0);
    nu_vtx_y_truth                         = truth_particle.Vy(0);
    nu_vtx_z_truth                         = truth_particle.Vz(0);

  }

  if ( neutrino_interaction_channel != 0 )
    return false;

  if ( neutrino_PDG != -14 && neutrino_PDG != 14 ) 
    return false;

  if ( nu_vtx_x_truth < 20. || nu_vtx_x_truth > 236.35 || nu_vtx_y_truth < -96.5 || nu_vtx_y_truth > 96.5 || nu_vtx_z_truth < 20. || nu_vtx_z_truth > 1016.8 )
    return false;

  if ( neutrino_energy < 1000. )
    return false;

  return true;

}
DEFINE_ART_MODULE(lar::HighEnergyNuMINuMuCCFilter)
