/**
 * \file FakePhotons_module.cc
 *
 * \ingroup OpticalDetectorSim
 * 
 * \brief Class def header for a class FakePhotons_module.cc
 *
 * @author kazuhiro
 */

/** \addtogroup OpticalDetectorSim

@{*/

#ifndef FAKEPHOTONS_MODULE_CC
#define FAKEPHOTONS_MODULE_CC

// LArSoft includes
#include "lardataobj/Simulation/SimPhotons.h"
#include "larcoreobj/SummaryData/RunData.h"
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

namespace opdet {

  class FakePhotons : public art::EDProducer{
  public:
    
    FakePhotons(const fhicl::ParameterSet&);
    virtual ~FakePhotons() {};
    
    // This method reads in any parameters from the .fcl files. 
    void reconfigure(fhicl::ParameterSet const& pset);
    
    void produce(art::Event&);

    void beginRun(art::Run&);
      
  private:  
    void addPE(sim::SimPhotons& simph,double g4time);
    /// Photo 2D array [opdet id][G4 ns]
    std::vector<std::vector<double> > _photons;
    double _start, _end, _period;
    int _num_pe;
    std::vector<int> _ch_v;
    int _make_rundata;
  };
} // namespace opdet


// Required for any LArSoft module.
namespace opdet {
  DEFINE_ART_MODULE(FakePhotons)
} 


// Implementation
namespace opdet {
  
  /// ------------------------------------------------------------------------------------
  /// Constructor
  FakePhotons::FakePhotons(fhicl::ParameterSet const& parameterSet)
  {
    this->reconfigure(parameterSet);
    produces<std::vector<sim::SimPhotons> >();
    if(_make_rundata)
      produces< sumdata::RunData, art::InRun >();
  }

  /// ------------------------------------------------------------------------------------

  void FakePhotons::reconfigure(fhicl::ParameterSet const& p)
  {
    _num_pe = p.get<int>("PELevel",10);
    _make_rundata = p.get<int>("MakeRunData",0);
    ::art::ServiceHandle<geo::Geometry> geo;
    //auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();
    _photons.clear();
    _photons.resize((int)geo->NOpDets());
    std::vector<double> empty;
    for(size_t i=0; i<_photons.size(); ++i) {
      std::string key = "PMT" + std::to_string(i);
      _photons[i] = p.get< std::vector<double> >(key,empty);
    }

    _start = p.get<double>("Start",-1.e6);
    _end   = p.get<double>("End",_start - 1.);
    _period = p.get<double>("Period",0.);
    _ch_v  = p.get<std::vector<int> >("OpDetList",_ch_v);
    return;
  }

  //-------------------------------------------------

  void FakePhotons::beginRun(art::Run& run)
  {
    // grab the geometry object to see what geometry we are using

    if(_make_rundata) {
      ::art::ServiceHandle<geo::Geometry> geo;
      
      std::unique_ptr<sumdata::RunData> runData(new sumdata::RunData(geo->DetectorName()));
      
      run.put(std::move(runData));
    }
    return;
  }

  void FakePhotons::addPE(sim::SimPhotons& simph,double g4time) {
    sim::OnePhoton ph;
    ph.Time = g4time;    
    for(int i=0; i<_num_pe; ++i)
      simph.push_back(ph);
  }

  //-------------------------------------------------

  void FakePhotons::produce(art::Event& event)
  {
    std::unique_ptr< std::vector<sim::SimPhotons> > simph_v(new std::vector<sim::SimPhotons>);

    for(size_t opch=0; opch<_photons.size(); ++opch) {
      sim::SimPhotons simph;
      simph.SetChannel(opch);
      simph.reserve(_photons[opch].size() * _num_pe);
      for(auto const& t : _photons[opch])
	addPE(simph,t);
      simph_v->emplace_back(simph);
    }
    
    ::art::ServiceHandle<geo::Geometry> geo;
    if(_start < _end && _period > 0) {
      //std::cout<<"Firing periodic pes..."<<std::endl;
      double t = _start;
      while(t < _end) {
	//std::cout<<"Firing at " << t << std::endl;
	for(auto const& ch : _ch_v) {
	  if(ch<0 || ch >= (int)(geo->NOpDets())) continue;
	  //std::cout<<"Channel "<<ch<<std::endl;
	  auto& simph = (*simph_v)[ch];
	  addPE(simph, t);
	}
	t += _period;
      }
    }
    event.put(std::move(simph_v));
  }
}

#endif
/** @} */ // end of doxygen group 
