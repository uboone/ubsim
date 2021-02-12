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

  // Truth Muon Track Containment Requirement.
  bool is_truth_muon_contained;

  // Number of Points on Muon MCTrack.
  int num_mctrack_points;

  // Truth Muon Endpoint Coordinates.
  double truth_muon_starting_x_coordinate;
  double truth_muon_starting_y_coordinate;
  double truth_muon_starting_z_coordinate;
  double truth_muon_ending_x_coordinate;
  double truth_muon_ending_y_coordinate;
  double truth_muon_ending_z_coordinate;
  
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

  if ( parent_fvz > 15600. ) 
    return false;

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

  // Load in the MCTrack information.                                                                                                                                                                      
  art::Handle<std::vector<sim::MCTrack> > mctrack_h;
  e.getByLabel("mcreco", mctrack_h);

  // make sure MCTrack info looks good                                                                                                                                                                     
  if(!mctrack_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate MCTrack!"<<std::endl;
    throw std::exception();
  }

  std::cout << "An event passed my filter!" << std::endl;
  std::cout << "Number of mctracks = " << mctrack_h->size() << "." << std::endl;

  is_truth_muon_contained = false;

  for ( size_t mctrack_iter = 0; mctrack_iter < mctrack_h->size(); mctrack_iter++ ) {

    num_mctrack_points = mctrack_h->at( mctrack_iter ).size();

    // At the top of the loop, find out which type of track this is.                                                                                    
    if ( mctrack_h->at( mctrack_iter ).PdgCode() == 13 || mctrack_h->at( mctrack_iter ).PdgCode() == -13 ) {

      if ( num_mctrack_points > 0 ) {

        if ( fabs( nu_vtx_x_truth - mctrack_h->at( mctrack_iter ).at( 0 ).X() ) < 0.001 && fabs( nu_vtx_y_truth - mctrack_h->at( mctrack_iter ).at( 0 ).Y() ) < 0.001 && fabs( nu_vtx_z_truth - mctrack_h->at( mctrack_iter ).at( 0 ).Z() ) < 0.001 ) {

	  truth_muon_starting_x_coordinate = mctrack_h->at( mctrack_iter ).at( 0 ).X();
          truth_muon_starting_y_coordinate = mctrack_h->at( mctrack_iter ).at( 0 ).Y();
          truth_muon_starting_z_coordinate = mctrack_h->at( mctrack_iter ).at( 0 ).Z();
          truth_muon_ending_x_coordinate   = mctrack_h->at( mctrack_iter ).at( mctrack_h->at( mctrack_iter ).size() - 1 ).X();
          truth_muon_ending_y_coordinate   = mctrack_h->at( mctrack_iter ).at( mctrack_h->at( mctrack_iter ).size() - 1 ).Y();
          truth_muon_ending_z_coordinate   = mctrack_h->at( mctrack_iter ).at( mctrack_h->at( mctrack_iter ).size() - 1 ).Z();

          if ( truth_muon_starting_x_coordinate > 3.0 && truth_muon_starting_x_coordinate < 253.35 && truth_muon_starting_y_coordinate > -113.5 && truth_muon_starting_y_coordinate < 113.5 && truth_muon_starting_z_coordinate > 3.0 && truth_muon_starting_z_coordinate < 1033.8 && truth_muon_ending_x_coordinate > 3.0 && truth_muon_ending_x_coordinate < 253.35 && truth_muon_ending_y_coordinate > -113.5 && truth_muon_ending_y_coordinate < 113.5 && truth_muon_ending_z_coordinate > 3.0 && truth_muon_ending_z_coordinate < 1033.8 ) {
	   
	    is_truth_muon_contained = true;
	    break;
	  
	  }

	}

      }  

    }

  }

  // If the track is not contained, return false.
  if ( is_truth_muon_contained == false )
    return false;

  return true;

}
DEFINE_ART_MODULE(lar::HighEnergyNuMINuMuCCFilter)
