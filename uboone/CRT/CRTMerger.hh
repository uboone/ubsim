#ifndef CRT_MERGER_HH
#define CRT_MERGER_HH

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "artdaq-core/Data/Fragments.hh"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/InputTag.h"
//#include "gallery/Event.h"
//#include "gallery/ValidHandle.h"
#include <string>
#include <istream>

using namespace art;

namespace crt
{
	class CRTMerger : public art::EDProducer
	{
		std::vector<std::string> fFileNames;
		
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
		
		void produce(art::Event & evt) override;
		
		void beginJob() override;
		
		void endJob() override;
		
		void reconfigure(fhicl::ParameterSet const & p) override;
		
		private:
		
		std::vector< std::vector< artdaq::Fragment > > w;
	};
}
#endif // CRT_MERGER_HH
