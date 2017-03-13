/**
  Eats up CRTData and spits out BernZMQFragments

**/

#ifndef CRTDaqSim_HH_
#define CRTDaqSim_HH_

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"

#include <string>

namespace crt{
  class CRTDaqSim :  public art:: EDProducer{
    std::string fProducerName;

  public:

    /// Default ctor
    CRTDaqSim(const fhicl::ParameterSet&);

    /// Default dtor
    ~CRTDaqSim();

    /// art::EDProducer::produce implementation
    virtual void produce (art::Event&);

  };
}


#endif  //CRTDaqSim_HH_
