#ifndef CRT_MERGER_HH
#define CRT_MERGER_HH

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "artdaq-core/Data/Fragments.hh"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/InputTag.h"
//#include "gallery/Event.h"
//#include "gallery/ValidHandle.h"
#include <string>
#include <istream>
//#include "ifdh.h"
#include "IFDH_service.h"

using namespace art;
using namespace std;

namespace crt
{
	class CRTMerger : public art::EDProducer
	{
		std::vector<std::string> fFileNames;
		std::vector<std::string> fTestFiles;
		 
		// Producer tag of the CRT events
		art::InputTag fTag;
		std::string Merged_Object;
		
		// Time window
		unsigned long fTimeWindow;
		unsigned int fMaxCount;
		bool _debug;
		
		public:
		CRTMerger(const fhicl::ParameterSet&);
		~CRTMerger();
		//explicit CRTMerger(fhicl::ParameterSet const &p);
		
		void produce( art::Event &evt ) override;
		
		void reconfigure(fhicl::ParameterSet const & p) override;
		void bernfebdaq_ReadEvents(std::string filename, std::string output_file="output.root", int max_events=-1);
		private:
		
		std::vector< std::vector< artdaq::Fragment > > w;
		ifdh_ns::ifdh* fIFDH=0;
		ifdh_ns::ifdh* tIFDH=0;
		
		std::vector< std::string > fCRTFile;
		
		std::string fInputFileName;
		uint32_t fRun;
		uint32_t fsubRun;
		uint32_t fDataStartTime;
		uint32_t fDataEndTime;
		
  std::string          fSamDefname;
  std::string          fSamProject;
  std::string          fSamStation;
  std::string          fSamAppFamily;
  std::string          fSamAppName;
  std::string          fSamAppVersion;
  std::string          fSamUser;
  std::string          fSamDescription;
  int                  fSamFileLimit;
  std::string          fSamSchema;

  std::string          fSamProjectURI;
  std::string          fSamProcessID;
  std::string          fSamCurrentFileURI;
  std::string          fSamCurrentFileName;
	};
}
#endif // CRT_MERGER_HH
