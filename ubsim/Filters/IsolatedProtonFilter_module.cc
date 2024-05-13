////////////////////////////////////////////////////////////////////////
// Class:       IsolatedProtonFilter
// Plugin Type: filter (art v3_01_02)
// File:        IsolatedProtonFilter_module.cc
//
// This is the filter module for the Isolated Proton Candidate 
// pre-selections. 
//
// Author: Jiaoyang Li/李 娇瑒 (Jiaoyang.li@ed.ac.uk)
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

#include "nusimdata/SimulationBase/MCParticle.h"
#include <memory>

namespace sim {
  class IsolatedProtonFilter;
}


class sim::IsolatedProtonFilter : public art::EDFilter {
public:
  explicit IsolatedProtonFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  IsolatedProtonFilter(IsolatedProtonFilter const&) = delete;
  IsolatedProtonFilter(IsolatedProtonFilter&&) = delete;
  IsolatedProtonFilter& operator=(IsolatedProtonFilter const&) = delete;
  IsolatedProtonFilter& operator=(IsolatedProtonFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:

  // Declare member data here.
  art::InputTag fMCParticleLabel;
  int fPDG_code;
  float fXmin, fXmax, fYmin, fYmax, fZmin, fZmax;
  float ftrack_len_cut, ftrack_len_x_cut;

};


sim::IsolatedProtonFilter::IsolatedProtonFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{
  fMCParticleLabel = p.get<art::InputTag>("MCParticleLabel");
  fPDG_code        = p.get<int>("PDGCode");
  fXmin            = p.get<float>("Xmin");
  fXmax            = p.get<float>("Xmax");
  fYmin            = p.get<float>("Ymin");
  fYmax            = p.get<float>("Ymax");
  fZmin            = p.get<float>("Zmin");
  fZmax            = p.get<float>("Zmax");
  ftrack_len_cut   = p.get<float>("track_len_cut");
  ftrack_len_x_cut = p.get<float>("track_len_x_cut");
}

bool sim::IsolatedProtonFilter::filter(art::Event& e)
{
  // Get MCParticles in the event ...
  art::Handle<std::vector<simb::MCParticle>> mcparHandle;
  e.getByLabel(fMCParticleLabel, mcparHandle);
  std::vector<simb::MCParticle> const& mcparVec(*mcparHandle);
  bool isTargetedProton = false;
  
  for (auto const& mcpar : mcparVec){
    float track_x_length_abs=std::abs(mcpar.Vx()-mcpar.EndX());
    float track_length=std::hypot(mcpar.Vx()-mcpar.EndX(), mcpar.Vy()-mcpar.EndY(), mcpar.Vz()-mcpar.EndZ());
    
    if (mcpar.PdgCode()==fPDG_code && mcpar.Mother()==0){ // check if it is the primary proton
      // If particle is within the fiducial volume cut
      if ((mcpar.Vx() > fXmin && mcpar.EndX() > fXmin && mcpar.Vx() < fXmax && mcpar.EndX() < fXmax) &&
          (mcpar.Vy() > fYmin && mcpar.EndY() > fYmin && mcpar.Vy() < fYmax && mcpar.EndY() < fYmax) &&
          (mcpar.Vz() > fZmin && mcpar.EndZ() > fZmin && mcpar.Vz() < fZmax && mcpar.EndZ() < fZmax)){ 
          // If the particle passes the track length selection and the selection of track length in the x-direction.  
          if (track_length < ftrack_len_cut && track_x_length_abs < ftrack_len_x_cut) { 
            isTargetedProton=true;
          }
      }
    }
  }
  return isTargetedProton;
}

DEFINE_ART_MODULE(sim::IsolatedProtonFilter)
