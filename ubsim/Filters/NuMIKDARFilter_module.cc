////////////////////////////////////////////////////////////////////////
// Class:       NuMIKDARFilter
// Plugin Type: filter (art v2_11_03)
// File:        NuMIKDARFilter_module.cc
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
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// services etc...                                                                                                                                                                                     
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// data-products                                                                                                                                                                                            
// lardataobj                                                                                                                                                                                           
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/Utilities/AssociationUtil.h"

// nusimdata                                                                                                                                                                                           
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// ROOT                                                                                                                                                                                                  
#include "TVector3.h"

// C++                                                                                                                                                                                                  
#include <memory>
#include <iostream>
#include <utility>

namespace sim {
  class NuMIKDARFilter;
}


class sim::NuMIKDARFilter : public art::EDFilter {
public:
  explicit NuMIKDARFilter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuMIKDARFilter(NuMIKDARFilter const &) = delete;
  NuMIKDARFilter(NuMIKDARFilter &&) = delete;
  NuMIKDARFilter & operator = (NuMIKDARFilter const &) = delete;
  NuMIKDARFilter & operator = (NuMIKDARFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

private:

  // Declare member data here.
  // producer of 3D reconstructed track to be used                                                                                                                                              
  int         event_ctr; // Count the events you are looping over.                                                                                                                                   

  bool        _debug;

  bool        kaon_parent_from_dump;
  bool        kdar_energy;
  bool        CCNC_interaction;
  bool        neutrino_contained;

  // Declare variables for the neutrino coordinates.                                                                                                                                                      
  double nu_vtx_x_truth;
  double nu_vtx_y_truth;
  double nu_vtx_z_truth;

};


sim::NuMIKDARFilter::NuMIKDARFilter(fhicl::ParameterSet const & p)
: EDFilter(p)
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  event_ctr = 0;
  _debug    = false;
}

bool sim::NuMIKDARFilter::filter(art::Event &e)
{
  // Set the four booleans to false;                                                                                                                                                                  
  kaon_parent_from_dump = false;
  kdar_energy           = false;
  CCNC_interaction      = false;
  neutrino_contained    = false;

  std::cout << "Event #" << event_ctr << "." << std::endl;
  event_ctr++;

  // load neutrino information.                                                                                                                                                                         
  if (_debug) { std::cout << "loading neutrino from producer generator." << std::endl; }
  art::Handle<std::vector<simb::MCTruth> > neutrino_h;
  e.getByLabel("generator", neutrino_h);

  // make sure tracks look good                                                                                                                                                                        
  if(!neutrino_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Neutrino!"<<std::endl;
    throw std::exception();
  }

  // Declare a vector of pointers for the neutrino object.                                                                                                                                             
  std::vector<art::Ptr<simb::MCTruth> >NeutrinoVec;
  art::fill_ptr_vector(NeutrinoVec, neutrino_h);

  // load the truth MCFlux information.                                                                                                                                                                
  if (_debug) { std::cout << "loading MCFlux object from producer generator." << std::endl; }
  art::Handle<std::vector<simb::MCFlux> > parent_h;
  e.getByLabel("generator", parent_h);

  // make sure tracks look good                                                                                                                                                                        
  if(!parent_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Parent!" << std::endl;
    throw std::exception();
  }

  // Declare a vector of pointers with the 'MCFlux' object.                                                                                                                                             
  std::vector<art::Ptr<simb::MCFlux> > ParentVec;
  art::fill_ptr_vector(ParentVec, parent_h);

  // Print out the lengths of 'neutrino_h' and 'parent_h'.                                                                                                                                               
  std::cout << "Length of 'MCFlux' object in this event = " << parent_h->size() << "." << std::endl;
  std::cout << "Length of 'MCTruth' object in this event = " << neutrino_h->size() << "." << std::endl;

  // Loop through the neutrino parents.
  kaon_parent_from_dump = false;

  for (auto& parent : ParentVec){

    // Continue if the parent is not a kaon and/or it does not decay at the dump.                                                                                                                       
    if ( parent->fptype == 321 && parent->fvz > 72300 && parent->fvz < 72800 ) {
      kaon_parent_from_dump = true;
      break;
    }

  }

  if ( kaon_parent_from_dump == false ) 
    return false;

  // loop through neutrinos themselves.                                                                                                                                                                
  for (auto& neutrino : NeutrinoVec ) {

    // Reset 'kdar_energy', 'CCNC_interaction', and 'neutrino_contained' to true.                                                                                  
    kdar_energy        = true;
    CCNC_interaction   = true;
    neutrino_contained = true;

    // Unpack the neutrino object to find an MCParticle.                                                                                                                                              
    const simb::MCNeutrino& truth_neutrino = neutrino->GetNeutrino();
    const simb::MCParticle& truth_particle = truth_neutrino.Nu();

    // Unpack the coordinates for the vertex as well.                                                                                                                                                 
    nu_vtx_x_truth               = truth_particle.Vx(0);
    nu_vtx_y_truth               = truth_particle.Vy(0);
    nu_vtx_z_truth               = truth_particle.Vz(0);

    // Only look at those events which decay via the charged current channel.                                                                                                                         
    if ( truth_neutrino.CCNC() != 0 )
      CCNC_interaction = false;
    
    // Continue if the particle's energy is outside of the correct peak.                                                                                                                               
    if ( truth_particle.E(0) < 0.2355 || truth_particle.E(0) > 0.2356 ) {
      kdar_energy = false;
    }

    // Check to see if the neutrino vertex is contained.                                                                                                                                              
    if ( nu_vtx_x_truth < 0.0 || nu_vtx_x_truth > 256.4 || nu_vtx_y_truth < -116.5 || nu_vtx_y_truth > 116.5 || nu_vtx_z_truth < 0.0 || nu_vtx_z_truth > 1036.8 ) {
      neutrino_contained = false;
    }

    // Return 'true' if all of the conditions are met.
    if ( kdar_energy && CCNC_interaction && neutrino_contained ) {
      std::cout << "Hurray! We found an event with a KDAR neutrino." << std::endl;
      return true;
    }

  } // End of the loop over the neutrino candidates in the event.

  // Return false otherwise - the conditions in the neutrino loop were not met.
  return false;

}

DEFINE_ART_MODULE(sim::NuMIKDARFilter)
