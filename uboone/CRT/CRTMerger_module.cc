#include "uboone/CRT/CRTMerger.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "uboone/CRT/CRTProducts/CRTSimData.hh"
#include "uboone/CRT/CRTProducts/MSetCRTFrag.hh"
#include "artdaq-core/Data/Fragments.hh"
#include "bernfebdaq-core/Overlays/BernZMQFragment.hh"

namespace crt{

	CRTMerger::CRTMerger(const fhicl::ParameterSet& pset): 
		fFileNames(pset.get<std::vector<std::string> >("InputFilenames"))
	{
		this->reconfigure(pset);
		produces< std::vector<CRTSimData> >();
		
		fTag		= pset.get<art::InputTag> ("InputTagName");
		fMerged_Object	= pset.get<std::string> ("ObjectProducer");
		_debug		= pset.get<bool>		 ("debug");
	}

	CRTMerger::~CRTMerger()
	{

	}


	void beginRun(art::Run &)
	{

	}

	void endRun(art::Run &)
	{

	}

	void CRTMerger::produce(art::Event& event)
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
		std::cout<<"TPC event time"<<event_t<<"   event time (s) "<<evt_time_sec<<"   event time (ns) "<<evt_time_nsec<<std::endl;
		
		mf::LogInfo("CRTMerger")<<"Finding CRT Events for time : "<<event_t<<std::endl;;
		
		//while(!fCRTEvent.atEnd())
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
				
				//std::cout<<"TPC (s) "<<evt_time_sec<<"     CRT(s) "<<crt_bf_time_s<<"     TPC (ns) "<<evt_time_nsec<<"     CRT_start(ns) "<<crt_bf_time_start_ns<<"     CRT_end(ns) "<<crt_bf_time_end_ns<<std::endl;
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
						crt::MSetCRTFrag FragCRT2Merge(bernfrag, evt_time_sec, evt_time_nsec, crt_bf_time_start_ns, crt_bf_time_end_ns);
						MergedCRTFragSet->push_back(FragCRT2Merge);
					}
				}
				/*
				if (evt_time_sec < crt_bf_time_s)
				{
					std::cout<<"Go to next TPC event"<<std::endl;
					return;
				}
				else if ((evt_time_sec - crt_bf_time_s)==0)
				{
					if (crt_bf_time_start_ns < evt_time_nsec)
					{
						std::cout<<"crt_ns less than tpc_ns "<<crt_bf_time_start_ns<<std::endl;
						if (abs(crt_bf_time_start_ns - evt_time_nsec)>fTimeWindow)
						{
							//fCRTEvent.next();
							std::cout<<"After calling next, start"<<crt_bf_time_start_ns<<"     stop "<<crt_bf_time_end_ns<<std::endl;
						}
						else if (abs(crt_bf_time_start_ns - evt_time_nsec)<fTimeWindow)
						{
							std::cout<<"perform merging"<<std::endl;
							break;
						}
						else{}
					}
					else
					break;	
				}
				else
				{
					std::cout<<"Go to next TPC event"<<std::endl;
					break;
				}
				*/
			}
			if (xxx==1)
			{
				std::cout<<"escape 2"<<std::endl;
				break;
			}
		}
		
		
		/*
		while(!fCRTEvent.atEnd()){
			for(auto const& crtFrag : crthits){
				bernfebdaq::BernZMQFragment bernfrag(crtFrag);
				auto bfmetadata = bernfrag.metadata();
				//auto feb_id = mtdata->feb_id();
				//auto pc_id = bernfrag.metadata()->reader_id();
				
				//std::cout<<"TPC event time:(s) "<<evt_time_sec<<",   (ns) "<<evt_time_nsec<<"     "<<evt_time_sec+evt_time_nsec<<std::endl;
				//std::cout<<"CRT event time:(s) "<<metadata->time_start_seconds()<<",   (s) "<<metadata->time_end_seconds()<<",   (ns) "<<metadata->time_start_nanosec()<<"     "<<metadata->time_end_nanosec()<<"     "<<metadata->time_start_seconds()+metadata->time_start_nanosec()<<std::endl;
				
				for(size_t i_e=0; i_e<bfmetadata->n_events(); ++i_e){
					bernfebdaq::BernZMQEvent const* this_event = bernfrag.eventdata(i_e); //single hit/event 
					auto time_crt_ns = this_event->Time_TS0();
					
					if (abs(evt_time_sec-bfmetadata->time_start_seconds())<5)
					std::cout<<"TPC (s) "<<evt_time_sec<<"     CRT(s) "<<bfmetadata->time_start_seconds()<<"     TPC (ns) "<<evt_time_nsec<<"     CRT(ns) "<<time_crt_ns<<std::endl;
				}
				
				
				//std::cout<<feb_id<<"     "<<pc_id<<"     "<<bernfrag.metadata()->n_events()<<std::endl;
				for(size_t i = 0; i<bernfrag.metadata()->n_events() ; i++){
					bernfebdaq::BernZMQEvent const* this_event = bernfrag.eventdata(i);
					
					bool sptevent0 = this_event->IsReference_TS0();
					bool sptevent1 = this_event->IsReference_TS1();
					std::cout<<"..................................................."<<std::endl;
					std::cout<<"sptevent0 "<<sptevent0<<"     sptevent1 "<<sptevent1<<std::endl;
					
					if ((!sptevent0) && (!sptevent1))
					{
						auto time_ts0 = this_event->Time_TS0();
						//double corrected_time = GetCorrectedTime(time_ts0,*mtdata);
						//double time_tevt_s = mtdata->time_start_seconds();
						//double time_tevt_ns = corrected_time;
						std::cout<<"CRT event timing: "<<mtdata->time_start_seconds()<<"     "<<mtdata->time_end_seconds()<<"     "<<this_event->Time_TS0()<<"     "<<GetCorrectedTime(time_ts0,*mtdata)<<"   --> "<<mtdata->time_start_seconds()+this_event->Time_TS0()<<std::endl;
						std::cout<<"TPC event timing: "<<evt_time_sec<<"     "<<evt_time_nsec<<"   --> "<<evt_time_sec+evt_time_nsec<<std::endl;
					}
					auto crt_time = bernfrag.eventdata(i)->Time_TS0();
					
					//auto strttime = mtdata->time_start_seconds();
					//auto end_time = mtdata->time_end_seconds();
					//double tt = bernfebdaq::GetCorrectedTime(crt_time,*mtdata);
					
					if(!_debug)
					std::cout<<evt_time_sec<<"     "<<evt_time_nsec<<"     "<<crt_time<<"     "<<bernfrag.eventdata(i)->Time_TS1()<<std::endl;
					if (crt_time>= evt_time_nsec/10.0-fTimeWindow/2 && crt_time< evt_time_nsec/10.0+fTimeWindow/2 ){
						for(auto adc_channel =0; adc_channel<32; adc_channel++){

							unsigned channel = pc_id*1000+feb_id*100+adc_channel;
							crt::CRTSimData myNewHit(channel, bernfrag.eventdata(i)->Time_TS0(), bernfrag.eventdata(i)->Time_TS1(), bernfrag.eventdata(i)->ADC(adc_channel));
							crtHits->push_back(myNewHit);std::cout<<"new hit"<<std::endl;
						}
					}
					else if (crt_time>=evt_time_nsec/10.0+fTimeWindow/2){
						event.put(std::move(crtHits));
						return;
					}
				}
			}
			fCRTEvent.next();
		}*/
		//std::cout<<"size of crtHits "<<crtHits.size()<<std::endl;
		event.put(std::move(MergedCRTFragSet));
	}

	void CRTMerger::reconfigure(fhicl::ParameterSet const & pset)
	{
		//fCRTEvent.toBegin();
		fTag = {pset.get<std::string>("InputTagName","crtdaq")};
		fTimeWindow = pset.get<unsigned>("TimeWindow",5000000);
		crt_start_ns= pset.get<unsigned>("start_ns", 0);
		crt_end_ns  = pset.get<unsigned>("end_ns", 5000000);
	}


	/*
	   void CRTMerger::beginRun(art::Run const&)
	   {

	   }
	   */
}

DEFINE_ART_MODULE(crt::CRTMerger)
