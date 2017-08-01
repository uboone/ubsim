#include "uboone/CRT/CRTMerger.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "uboone/CRT/CRTProducts/MSetCRTFrag.hh"
#include "artdaq-core/Data/Fragments.hh"
#include "gallery/Event.h"
#include "bernfebdaq-core/Overlays/BernZMQFragment.hh"
crt::CRTMerger::CRTMerger(const fhicl::ParameterSet& pset): fFileNames(pset.get<std::vector<std::string> >("InputFilenames"))
{
	this->reconfigure(pset);
	//produces< std::vector<crt::MSetCRTFrag> >();
	produces< std::vector<artdaq::Fragment> >();
	//produces< std::vector< std::vector<artdaq::Fragment> > >();
	fTag		= pset.get<art::InputTag> ("InputTagName");
	_debug		= pset.get<bool>		 ("debug");
}

crt::CRTMerger::~CRTMerger()
{
	//fIFDH->cleanup();
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
	//std::unique_ptr<std::vector<artdaq::Fragment> > ThisFragment(new std::vector<artdaq::Fragment>);
	//std::unique_ptr< std::vector <std::vector<artdaq::Fragment> > > VecofVecFrags(new std::vector<std::vector<artdaq::Fragment>>);
	//std::unique_ptr<std::vector<artdaq::Fragment> > LastFragment(new std::vector<artdaq::Fragment>);
	//std::unique_ptr<std::vector<artdaq::Fragment> > NextFragment(new std::vector<artdaq::Fragment>);
	
	std::unique_ptr<std::vector<artdaq::Fragment> > ThisFragSet(new std::vector<artdaq::Fragment>);
	std::unique_ptr<std::vector<artdaq::Fragment> > LastFragSet(new std::vector<artdaq::Fragment>);
	std::unique_ptr<std::vector<artdaq::Fragment> > NextFragSet(new std::vector<artdaq::Fragment>);
	
	//For this event
	//std::cout<<fMaxCount<<std::endl;
	std::vector<artdaq::Fragment>  ThisFragment;
	art::Timestamp evtTime = event.time();
	unsigned long evt_time_sec = evtTime.timeHigh();	
	unsigned long evt_time_nsec = evtTime.timeLow();
	//std::cout<<"TPC event time"<<event_t<<"   event time (s) "<<evt_time_sec<<"   event time (ns) "<<evt_time_nsec<<std::endl;
	
	unsigned int count = 0;
	unsigned int n = 0;
	//gallery::Event fCRTEvent(fFileNames);
	gallery::Event fCRTEvent(fTestFiles);
	for(fCRTEvent.toBegin(); !fCRTEvent.atEnd(); ++fCRTEvent)
	{
		int xxx = 0;
		
		auto const& crtFragSet = *(fCRTEvent.getValidHandle< std::vector<artdaq::Fragment> >(fTag));
		n=0;
		//for(auto const& i_crtFrag : crtFragSet)
		for(auto iv = begin(crtFragSet); iv != end(crtFragSet); ++iv)
		{
			auto const& i_crtFrag = *iv;
			//std::cout<<"count "<<count<<", n "<<n<<std::endl;
			bernfebdaq::BernZMQFragment bernfrag(i_crtFrag);
			
			auto i_Frag_metadata = bernfrag.metadata();
			unsigned long crt_bf_time_s = i_Frag_metadata->time_start_seconds();
			unsigned long crt_bf_time_e = i_Frag_metadata->time_end_seconds();
			unsigned long crt_bf_time_start_ns = i_Frag_metadata->time_start_nanosec();
			unsigned long crt_bf_time_end_ns = i_Frag_metadata->time_end_nanosec();
			unsigned long total_TPC_time = evt_time_sec*1000000000+evt_time_nsec;
			unsigned long time_CRT_start = crt_bf_time_s*1000000000+crt_bf_time_start_ns;
			unsigned long time_CRT_end   = crt_bf_time_e*1000000000+crt_bf_time_end_ns;
			
			//std::cout<<"TPC (s) "<<evt_time_sec<<"     CRT(s) "<<crt_bf_time_s<<"     TPC (ns) "<<evt_time_nsec<<"     CRT_start(ns) "<<crt_bf_time_start_ns<<"     CRT_end(ns) "<<crt_bf_time_end_ns<<std::endl;
			////////////////////////////////////////////////////////////////////////////////////////////////////
			if ((total_TPC_time>=time_CRT_start) && (total_TPC_time<=time_CRT_end))
                        {
				if (evt_time_sec != crt_bf_time_s)
				{
					xxx = 1;
					if (_debug)
					std::cout<<"escape 1"<<std::endl;
					break;
				}
				if (_debug)
				std::cout<<"merge:"<<evt_time_sec<<"     "<<crt_bf_time_s<<"    "<<evt_time_nsec<<"     ("<<crt_bf_time_start_ns<<", "<<crt_bf_time_end_ns<<")"<<std::endl;
				//ThisFragSet->emplace_back(i_crtFrag);
				xxx=1;
				
				artdaq::Fragment last_artdaqFrag = (count > 1) ? w[count-1][n]:i_crtFrag;
				bernfebdaq::BernZMQFragment last_bernfrag(last_artdaqFrag);
				auto last_bernfrag_metadata = last_bernfrag.metadata();
				if (_debug)
				std::cout<<"last Frag: start "<<last_bernfrag_metadata->time_start_nanosec()<<", stop: "<<last_bernfrag_metadata->time_end_nanosec()<<std::endl;
				LastFragSet->emplace_back(last_artdaqFrag);
				
				artdaq::Fragment next_artdaqFrag = (count < fMaxCount) ? w[count+1][n]:i_crtFrag;
				bernfebdaq::BernZMQFragment next_bernfrag(next_artdaqFrag);
				auto next_bernfrag_metadata = next_bernfrag.metadata();
				if (_debug)
				std::cout<<"next Frag: start: "<<next_bernfrag_metadata->time_start_nanosec()<<", stop: "<<next_bernfrag_metadata->time_end_nanosec()<<std::endl;
				NextFragSet->emplace_back(next_artdaqFrag);	
				ThisFragSet->emplace_back(last_artdaqFrag);
				ThisFragSet->emplace_back(i_crtFrag);
				ThisFragSet->emplace_back(next_artdaqFrag);
				//ThisFragSet->emplace_back(last_artdaqFrag);
				//ThisFragSet->emplace_back(i_crtFrag);
				//ThisFragSet->emplace_back(next_artdaqFrag);
				crt::MSetCRTFrag testSet(last_artdaqFrag,i_crtFrag,next_artdaqFrag);
				MergedCRTFragSet->emplace_back(testSet);
			}
			n++;
		}
		
		if (xxx==1)
		{
			if (_debug)
			std::cout<<"escape 2"<<std::endl;
			break;
		}
		count++;
	} // end CRT evt loop
	if (_debug)
	//std::cout<<"Set of Fragments: "<<ThisFragSet->size()<<std::endl;
	std::cout<<"Set of Fragments: ("<<LastFragSet->size()<<", "<<ThisFragSet->size()<<", "<<NextFragSet->size()<<"), MergedCRTFragSet: "<<MergedCRTFragSet->size()<<std::endl;
	//ec//std::cout << "Size of VecofVecFrags " << VecofVecFrags->size() << std::endl;
	//event.put(std::move(VecofVecFrags));
	event.put(std::move(ThisFragSet));
	//event.put(std::move(MergedCRTFragSet));
	if (_debug)
	std::cout<<"---X---"<<std::endl;
}

void crt::CRTMerger::reconfigure(fhicl::ParameterSet const & pset)
{
	fTag = {pset.get<std::string>("InputTagName","crtdaq")};
	fTimeWindow = pset.get<unsigned>("TimeWindow",5000000);
	fCRTFile = pset.get< std::vector < std::string > >("InputFilenames");
	if ( ! fIFDH ) fIFDH = new ifdh_ns::ifdh;
	std::string fetchedfile(fIFDH->fetchInput(fCRTFile[0]));
	fTestFiles.push_back(fetchedfile);
	gallery::Event testGallery4ifdh(fTestFiles);
	unsigned int N = 0;
	for(testGallery4ifdh.toBegin(); !testGallery4ifdh.atEnd(); ++testGallery4ifdh)
	{
		auto const& TryCRT_frags = *(testGallery4ifdh.getValidHandle< std::vector<artdaq::Fragment> >(fTag));
		N++;
		std::vector< artdaq::Fragment > f;
		
		for(auto iv=begin(TryCRT_frags); iv!=end(TryCRT_frags); ++iv)
		{
			auto const& fg = *iv;
			f.push_back(fg);
		}
		w.push_back(f);
	}
	fMaxCount=N;
	
	if ( ! fIFDH )
	delete fIFDH;
	//std::cout<<"fMaxCount "<<fMaxCount<<std::endl;
	
	/*
	gallery::Event fCRTEvent(fFileNames);
	unsigned int row = 0;
	unsigned int col = 0;
	for(fCRTEvent.toBegin(); !fCRTEvent.atEnd(); ++fCRTEvent)
	{
		//u.emplace_back(fCRTEvent.getValidHandle< std::vector<artdaq::Fragment> >(fTag));	
		//std::cout<<""<<typeid(fCRTEvent.getValidHandle< std::vector<artdaq::Fragment> >(fTag)).name()<<std::endl;
		//auto const& AllCRT_frags = *(fCRTEvent.getValidHandle< std::vector<artdaq::Fragment> >(fTag));
		//v.emplace_back(AllCRT_frags);
		//std::cout<<""<<typeid(AllCRT_frags).name()<<std::endl;
		row++;
		
		auto const& AllCRT_frags = *(fCRTEvent.getValidHandle< std::vector<artdaq::Fragment> >(fTag));
		//auto this_frag = begin(AllCRT_frags);
		//auto const& test = *this_frag;
		//v.push_back(test);
		
		std::vector< artdaq::Fragment > f;
		
		for(auto iv=begin(AllCRT_frags); iv!=end(AllCRT_frags); ++iv)
		{
			col++;
			auto const& fg = *iv;
			//u.push_back(fg);
			f.push_back(fg);
			
		}
		w.push_back(f);
	}
	fMaxCount=row;
	std::cout<<"fMaxCount "<<fMaxCount<<std::endl;
	*/
}

DEFINE_ART_MODULE(crt::CRTMerger)
