////////////////////////////////////////////////////////////////////////
// Class:       SimWireMicroBooNEAna
// Module Type: analyzer
// File:        SimWireMicroBooNEAna_module.cc
//
// Generated at Wed May 21 14:57:20 2014 by Matthew Toups using artmod
// from cetpkgsupport v1_05_04.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include <iostream>

namespace detsim {
  class SimWireMicroBooNEAna;
}

class detsim::SimWireMicroBooNEAna : public art::EDAnalyzer {
public:
  explicit SimWireMicroBooNEAna(fhicl::ParameterSet const & p);
  virtual ~SimWireMicroBooNEAna();

  void analyze(art::Event const & evt) override;


private:

  std::string fDigitModuleLabel;
  // Declare member data here.

};


detsim::SimWireMicroBooNEAna::SimWireMicroBooNEAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  fDigitModuleLabel= p.get< std::string >("DigitModuleLabel");
}

detsim::SimWireMicroBooNEAna::~SimWireMicroBooNEAna() {
  // Clean up dynamic memory and other resources here.
}

void detsim::SimWireMicroBooNEAna::analyze(art::Event const & evt) {

   art::Handle< std::vector<raw::RawDigit> > digitVecHandle;

   evt.getByLabel(fDigitModuleLabel, digitVecHandle);
   if(!digitVecHandle.isValid()) throw cet::exception("") << "NO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";

   if (!digitVecHandle->size())  std::cout << "BLAH\n";

   mf::LogInfo("CalWireMicroBooNE") << "CalWireMicroBooNE:: digitVecHandle size is " << digitVecHandle->size();

   std::cout << "CalWireMicroBooNE:: digitVecHandle size is " << digitVecHandle->size() << std::endl;

   //std::cout << "DigSize: " << digitVecHandle->size() << std::endl;
   // How many raw digits are there? 8192?
   for(raw::RawDigit const& digit : *digitVecHandle) {

     unsigned short mySize = digit.Samples();
     //std::cout << "Size: " << mySize << std::endl;
     std::vector<float> holder(mySize);
     std::vector<short> rawadc(mySize);

     //uint32_t channel = digit.Channel(); // unused

     for(size_t i = 0; i<digit.NADC(); i++) {
       //std::cout << "i: " << i << "\tfADC[" << i << "]: " << digit.ADC(i) << std::endl;
     }
     //std::cout << "fADC size: " << digit.NADC() << "\tEntry 0: " << digit.ADC(1236) << std::endl;
     // uncompress the data
     raw::Uncompress(digit.ADCs(), rawadc, digit.Compression());

     // loop over all adc values and subtract the pedestal
     // float pdstl = digit.GetPedestal(); // unused

     std::string adcinfo="";
     for(size_t bin = 0; bin < mySize; ++bin) {
       //if(rawadc[bin]>0) std::cout << "Bin: " << bin << "\trawadc[" << bin << "]: " << rawadc[bin] << std::endl;
       adcinfo += rawadc[bin];
       if(bin+1<mySize)
	 adcinfo += ",";
     }
     // channel += pdstl; // unused
     //std::cout << printf("Channel: %d\tPedestal: %f\tADCs:\n%s\n",channel,pdstl,adcinfo.c_str());
   }

   // Implementation of required member function here.
}

DEFINE_ART_MODULE(detsim::SimWireMicroBooNEAna)
