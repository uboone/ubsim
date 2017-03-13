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

#include <vector>
#include <memory>
#include <sstream>
#include <string>
#include <math.h>

namespace crt{

  CRTDaqSim::CRTDaqSim(const fhicl::ParameterSet& pSet)
  {
    produces< std::vector<CRTData> >();

  }

  CRTDaqSim::~CRTDaqSim()
  {

  }

  void CRTDaqSim::reconfigure(fhicl::ParameterSet const & pSet) {

    fProducerName = pSet.get<std::string>("ProducerName", "crtdetsim");
  }

  void CRTDaqSim::produce(art::Event& evt)
  {
    std::unique_ptr<std::vector<artdaq::Fragment> > crtHits(
        new std::vector<artdaq::Fragment>);

    evt.put(std::move(crtHits));
  }


}

DEFINE_ART_MODULE(crt::CRTDaqSim)