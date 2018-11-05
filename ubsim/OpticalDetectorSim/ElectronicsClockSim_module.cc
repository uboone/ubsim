/**
 * \file ElectronicsClockSim_module.cc
 *
 * \ingroup OpticalDetectorSim
 * 
 * \brief Class def header for a class ElectronicsClockSim_module.cc
 *
 * @author kazuhiro
 */

/** \addtogroup OpticalDetectorSim

@{*/

#ifndef ELECTRONICSCLOCKSIM_MODULE_CC
#define ELECTRONICSCLOCKSIM_MODULE_CC

// LArSoft includes
#include "larcoreobj/SummaryData/RunData.h"
#include "lardataobj/Simulation/ElectronicsClockInfo.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "larcore/Geometry/Geometry.h" // larcore

// ART includes.
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/Exception.h"

// C++ language includes
#include <iostream>
#include <sstream>
#include <cstring>
#include <vector>
#include <map>
#include <memory>
#include <cmath>

// Package include
#include "CLHEP/Random/RandFlat.h"
#include "TRandom.h"
#include "nutools/RandomUtils/NuRandomService.h"

namespace opdet {

  class ElectronicsClockSim : public art::EDProducer{
  public:
    
    ElectronicsClockSim(const fhicl::ParameterSet&);
    virtual ~ElectronicsClockSim() {delete _rand_gen;};
    
    // This method reads in any parameters from the .fcl files. 
    void reconfigure(fhicl::ParameterSet const& pset);
    
    void produce(art::Event&);

    void beginRun(art::Run&);
      
  private:  

    /// Photo 2D array [opdet id][G4 ns]
    double _start_g4time;
    double _random_shift;
    int _make_rundata;
    CLHEP::RandFlat* _rand_gen;
  };
} // namespace opdet


// Required for any LArSoft module.
namespace opdet {
  DEFINE_ART_MODULE(ElectronicsClockSim)
} 


// Implementation
namespace opdet {
  
  /// ------------------------------------------------------------------------------------
  /// Constructor
  ElectronicsClockSim::ElectronicsClockSim(fhicl::ParameterSet const& parameterSet)
  {
    _rand_gen = nullptr;
    this->reconfigure(parameterSet);
    produces< sim::ElectronicsClockInfo >();
    if(_make_rundata)
      produces< sumdata::RunData, art::InRun >();
  }

  /// ------------------------------------------------------------------------------------

  void ElectronicsClockSim::reconfigure(fhicl::ParameterSet const& p)
  {
    _start_g4time = p.get<double>("G4RefTime");
    _random_shift = p.get<double>("RandomShift",-1.);
    _make_rundata = p.get<int>("MakeRunData",0);
    art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, p, "Seed");
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    delete _rand_gen;
    _rand_gen = new CLHEP::RandFlat(engine);

    return;
  }

  //-------------------------------------------------

  void ElectronicsClockSim::beginRun(art::Run& run)
  {
    // grab the geometry object to see what geometry we are using
    if(_make_rundata) {
      ::art::ServiceHandle<geo::Geometry> geo;
      
      std::unique_ptr<sumdata::RunData> runData(new sumdata::RunData(geo->DetectorName()));
      
      run.put(std::move(runData));
    }
    return;
  }

  //-------------------------------------------------

  void ElectronicsClockSim::produce(art::Event& event)
  {
    double tstart = _start_g4time;

    if(_random_shift>0)
      tstart += _rand_gen->fire(0.,_random_shift);

    std::unique_ptr< sim::ElectronicsClockInfo > eclock(new sim::ElectronicsClockInfo(tstart));

    event.put(std::move(eclock));
  }
}

#endif
/** @} */ // end of doxygen group 
