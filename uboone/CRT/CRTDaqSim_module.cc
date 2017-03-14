#include "uboone/CRT/CRTData.hh"
#include "uboone/CRT/CRTDaqSim.hh"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "nutools/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "uboone/CRT/CRTData.hh"
#include "artdaq-core/Data/Fragments.hh"
#include "bernfebdaq-core/Overlays/BernZMQFragment.hh"


#include <vector>
#include <memory>
#include <sstream>
#include <string>
#include <math.h>

namespace crt{

  CRTDaqSim::CRTDaqSim(const fhicl::ParameterSet& pSet)
  {
    produces< std::vector<artdaq::Fragment> >();

  }

  CRTDaqSim::~CRTDaqSim()
  {

  }

  void CRTDaqSim::produce(art::Event& evt)
  {
    std::unique_ptr<std::vector<artdaq::Fragment> > frags(
        new std::vector<artdaq::Fragment>);

    //make the metadata 
    /*
      BernZMQFragmentMetadata(uint32_t ts_s, uint32_t ts_ns, 
        uint32_t te_s, uint32_t te_ns,
        int t_c, uint64_t t_o,
        uint32_t r, uint32_t seq,
        uint64_t fid, uint32_t rid,
        uint32_t nch, uint32_t nadc)


    bernfebdaq::BernZMQEvent zmqEvent;
      uint16_t mac5;
      uint16_t flags;
      uint16_t lostcpu;
      uint16_t lostfpga;
      uint32_t ts0;
      uint32_t ts1;
      uint16_t adc[32];

  frags.emplace_back( artdaq::Fragment::FragmentBytes(initial_payload_size,
                                                      ev_counter(), fragment_id(),
                                                      fragment_type_, metadata) );

    // Performs the hardware coincidence
  */
    evt.put(std::move(frags));
  }


}

DEFINE_ART_MODULE(crt::CRTDaqSim)