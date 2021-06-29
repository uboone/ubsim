////////////////////////////////////////////////////////////////////////
// Class:       StopMuFilter
// Plugin Type: filter (art v3_01_02)
// File:        StopMuFilter_module.cc
//
// Generated at Fri Feb  5 09:29:37 2021 by David Caratelli using cetskelgen
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

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larcore/Geometry/Geometry.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include <memory>

class StopMuFilter;


class StopMuFilter : public art::EDFilter {
public:
  explicit StopMuFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  
  // Plugins should not be copied or assigned.
  StopMuFilter(StopMuFilter const&) = delete;
  StopMuFilter(StopMuFilter&&) = delete;
  StopMuFilter& operator=(StopMuFilter const&) = delete;
  StopMuFilter& operator=(StopMuFilter&&) = delete;
  
  // Required functions.
  bool filter(art::Event& e) override;
  
  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  
private:
  
  // Declare member data here.
  art::InputTag fMCPproducer;
  
  bool Contained(const float& x, const float& y, const float& z);
  
};


StopMuFilter::StopMuFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
// More initializers here.
{
  fMCPproducer = p.get< art::InputTag > ("MCPproducer");
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

bool StopMuFilter::filter(art::Event& e)
{
  // Implementation of required member function here.
  
  
  // load  MCParticle from largeant
  auto const& mcp_h = e.getValidHandle<std::vector<simb::MCParticle> >(fMCPproducer);
  
  // loop through MCParticles and cut on whether end position is in TPC boundary
  // and on in-detector track-length
  bool stopmu = false;
  
  for (size_t p=0; p < mcp_h->size(); p++) {
    
    auto const& mcp = mcp_h->at(p);
    
    if ( fabs(mcp.PdgCode()) != 13) continue;
    
    //std::cout << "Found a muon!" << std::endl;
    
    // does the track cross the cathode?
    //for each trajectory step, check if cathode is crossed
    
    // get particle time for drift velocity offset
    //auto time = mcp.T(); // in ns
    //auto xoffset = 1.1098 * (1e-4) * time; // drift velocity in cm / ns
    
    // grab last trajectory point of muon
    int npt = mcp.NumberTrajectoryPoints();
    
    auto xstop = mcp.Vx(npt-1);
    auto ystop = mcp.Vy(npt-1);
    auto zstop = mcp.Vz(npt-1);
    
    //auto tstop = mcp.T(npt-1);
    
    // is this last point of the muon within the TPC volume?
    if ( Contained(xstop,ystop,zstop) == false) continue;
    
    // made it this far? stopping muon
    stopmu = true;
    
  }// for all MCParticles
  
  if (stopmu == true) return true;

  return false;
  
}

bool StopMuFilter::Contained(const float& x, const float& y, const float& z) {
  
  art::ServiceHandle<geo::Geometry> geo;
  geo::TPCGeo const &thisTPC = geo->TPC();
  geo::BoxBoundedGeo theTpcGeo = thisTPC.ActiveBoundingBox();
  std::vector<double> bnd = {theTpcGeo.MinX(), theTpcGeo.MaxX(), theTpcGeo.MinY(), theTpcGeo.MaxY(), theTpcGeo.MinZ(), theTpcGeo.MaxZ()};
  bool is_x = (x > bnd[0] && x < bnd[1]);
  bool is_y = (y > bnd[2] && y < bnd[3]);
  bool is_z = (z > bnd[4] && z < bnd[5]);
  
  return is_x && is_y && is_z;
}


void StopMuFilter::beginJob()
{
  // Implementation of optional member function here.
}

void StopMuFilter::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(StopMuFilter)
