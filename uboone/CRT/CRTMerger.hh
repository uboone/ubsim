#ifndef CRT_MERGER_HH
#define CRT_MERGER_HH

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include <string>
#include <istream>

namespace crt
{
	class CRTMerger : public art::EDProducer
	{
		std::vector<std::string> fFileNames;
		
		// Producer tag of the CRT events
		art::InputTag fTag;
		std::string Merged_Object;
		
		// Time window
		unsigned fTimeWindow;
		
		bool _debug;
		
		public:
		CRTMerger(const fhicl::ParameterSet&);
		
		~CRTMerger();
		
		void produce(art::Event & evt) override;
		
		void beginJob() override;
		
		void endJob() override;
		
		void reconfigure(fhicl::ParameterSet const & p) override;
	};
}
#endif // CRT_MERGER_HH
