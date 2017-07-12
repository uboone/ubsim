#include "uboone/CRT/CRTMerger.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "uboone/CRT/CRTProducts/MSetCRTFrag.hh"
#include "artdaq-core/Data/Fragments.hh"
#include "bernfebdaq-core/Overlays/BernZMQFragment.hh"

crt::CRTMerger::CRTMerger(const fhicl::ParameterSet& pset): fFileNames(pset.get<std::vector<std::string> >("InputFilenames"))
{
	this->reconfigure(pset);
	produces< std::vector<crt::MSetCRTFrag> >();	
	fTag		= pset.get<art::InputTag> ("InputTagName");
	//fMerged_Object= pset.get<std::string> ("ObjectProducer");
	_debug		= pset.get<bool>		 ("debug");
}

crt::CRTMerger::~CRTMerger()
{
}

void crt::CRTMerger::beginJob()
{
}

void crt::CRTMerger::endJob()
{
}

void crt::CRTMerger::produce(art::Event& event)
{
	if (_debug)
	{
		std::cout << "NEW EVENT" << std::endl;
	}
	
	std::unique_ptr<std::vector<crt::MSetCRTFrag> > MergedCRTFragSet(new std::vector<crt::MSetCRTFrag>);
	
	//For this event
	art::Timestamp evtTime = event.time();
	auto evt_time_sec = evtTime.timeHigh();	
	auto evt_time_nsec = evtTime.timeLow();
	auto event_t =  event.time().timeHigh();
	//std::cout<<"TPC event time"<<event_t<<"   event time (s) "<<evt_time_sec<<"   event time (ns) "<<evt_time_nsec<<std::endl;
	
	mf::LogInfo("CRTMerger")<<"Finding CRT Events for time : "<<event_t<<std::endl;;
	
	for(gallery::Event fCRTEvent(fFileNames); !fCRTEvent.atEnd(); fCRTEvent.next())
	{
		int xxx = 0;
		auto const& crtFragSet = *(fCRTEvent.getValidHandle< std::vector<artdaq::Fragment> >(fTag));
		for(auto const& i_crtFrag : crtFragSet)
		{
			bernfebdaq::BernZMQFragment bernfrag(i_crtFrag);
			auto i_Frag_metadata = bernfrag.metadata();
			auto crt_bf_time_s = i_Frag_metadata->time_start_seconds();
			auto crt_bf_time_start_ns = i_Frag_metadata->time_start_nanosec();
			auto crt_bf_time_end_ns = i_Frag_metadata->time_end_nanosec();
			
			std::cout<<"TPC (s) "<<evt_time_sec<<"     CRT(s) "<<crt_bf_time_s<<"     TPC (ns) "<<evt_time_nsec<<"     CRT_start(ns) "<<crt_bf_time_start_ns<<"     CRT_end(ns) "<<crt_bf_time_end_ns<<std::endl;
			
			if ((evt_time_sec - crt_bf_time_s)!=0)
			{
				xxx=1;
				std::cout<<"escape 1"<<std::endl;
				break;
			}
			if ((evt_time_sec - crt_bf_time_s)==0)
			{
				std::cout<<"---> "<<evt_time_sec<<"     "<<crt_bf_time_s<<"    "<<evt_time_nsec<<"     ("<<crt_bf_time_start_ns<<", "<<crt_bf_time_end_ns<<")     "<<xxx<<std::endl;
				
				if(abs(crt_bf_time_start_ns - evt_time_nsec)<fTimeWindow)
				{
					std::cout<<"merge:"<<std::endl;
					//crt::MSetCRTFrag FragCRT2Merge(bernfrag, evt_time_sec, evt_time_nsec, crt_bf_time_start_ns, crt_bf_time_end_ns);
					crt::MSetCRTFrag FragCRT2Merge(i_crtFrag);
					MergedCRTFragSet->emplace_back(FragCRT2Merge);
					xxx=1;
					break;
				}		
			}
		}
		if (xxx==1)
		{
			std::cout<<"escape 2"<<std::endl;
			break;
		}
	}
	if (MergedCRTFragSet->size()>0)
	std::cout<<"----> "<<MergedCRTFragSet->size()<<std::endl;
	event.put(std::move(MergedCRTFragSet));
}

void crt::CRTMerger::reconfigure(fhicl::ParameterSet const & pset)
{
	fTag = {pset.get<std::string>("InputTagName","crtdaq")};
	fTimeWindow = pset.get<unsigned>("TimeWindow",5000000);
}

DEFINE_ART_MODULE(crt::CRTMerger)
