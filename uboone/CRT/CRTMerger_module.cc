#include "uboone/CRT/CRTMerger.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "uboone/CRT/CRTProducts/MSetCRTFrag.hh"
#include "artdaq-core/Data/Fragments.hh"
#include "gallery/Event.h"
#include "bernfebdaq-core/Overlays/BernZMQFragment.hh"
#include "art/Framework/Services/System/FileCatalogMetadata.h"
#include "art/Framework/Principal/Handle.h"
#include "IFDH_service.h"

#include <memory>
#include <string>
#include <vector>
#include <exception>
#include <sstream>
#include <unistd.h>
#include <ctime>

//some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include <boost/date_time/c_local_time_adjustor.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

crt::CRTMerger::CRTMerger(const fhicl::ParameterSet& pset): fFileNames(pset.get<std::vector<std::string> >("InputFilenames")),
  fSamDefname(pset.get<std::string>("SamDefname", "")),
  fSamProject(pset.get<std::string>("SamProject", "")),
  fSamStation(pset.get<std::string>("SamStation", "")),
  fSamAppFamily(pset.get<std::string>("SamAppFamily", "art")),
  fSamAppName(pset.get<std::string>("SamAppName", "mix")),
  fSamAppVersion(pset.get<std::string>("SamAppVersion", "1")),
  fSamUser(pset.get<std::string>("SamUser", "")),
  fSamDescription(pset.get<std::string>("SamDescription", "sam-description")),
  fSamFileLimit(pset.get<int>("SamFileLimit", 100)),
  fSamSchema(pset.get<std::string>("SamSchema", "root"))
{
	std::cout<<"1 crt::CRTMerger::CRTMerger"<<std::endl;
	setenv("TZ", "CST+6CDT", 1);
	tzset();
	
	setenv("IFDH_DATA_DIR","/uboone/data/users/kolahalb/MicroBooNE/",1);
	//std::cout<<"ifdh_data_dir "<<getenv("IFDH_DATA_DIR")<<std::endl;
	
	this->reconfigure(pset);
	produces< std::vector<artdaq::Fragment> >();
	fTag		= pset.get<art::InputTag> ("InputTagName");
	_debug		= pset.get<bool>		 ("debug");
	if ( ! tIFDH ) tIFDH = new ifdh_ns::ifdh;
}

crt::CRTMerger::~CRTMerger()
{
	//fIFDH->cleanup();
	std::cout<<"8 crt::CRTMerger::~CRTmerger"<<std::endl;
	if(!fSamProcessID.empty())
	{
		if(!fSamCurrentFileName.empty())
		{
			tIFDH->updateFileStatus(fSamProjectURI,
			fSamProcessID,
			fSamCurrentFileName,
			"consumed");
		}
		tIFDH->endProcess(fSamProjectURI, fSamProcessID);
	}
	
	if ( ! tIFDH )
	delete tIFDH;
}

void crt::CRTMerger::produce(art::Event& event)
{
	if (_debug)
	{
		std::cout <<"5 crt::CRTMerger::produce, NEW EVENT" << std::endl;
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
	unsigned long time_tpc = evt_time_sec*1000000000 + evt_time_nsec;
	unsigned long time_tpc1= time_tpc/1000;
	std::cout<<"time_tpc1 "<<time_tpc1<<std::endl;
	
	const char* tz = getenv("TZ");
	std::string tzs(tz); std::cout<<"time-zone "<<tzs<<std::endl;
	/*
	if(tz != 0 && *tz != 0)
	{
		std::string tzs(tz);
		if(tzs != std::string("CST+6CDT"))
		{
			// Timezone is wrong, throw exception.
			throw cet::exception("CRTMerger") << "Wrong timezone: " << tzs;
		}
	}
	else
	{
		// Timezone is not set.  Throw exception.
		throw cet::exception("CRTMerger") << "Timezone not set.";
	}
	*/
	boost::posix_time::ptime const time_epoch(boost::gregorian::date(1970, 1, 1));
	boost::posix_time::ptime this_event_time = time_epoch + boost::posix_time::microseconds(time_tpc1);
	std::cout << boost::posix_time::to_iso_extended_string(this_event_time)<<std::endl;
	std::string stringTime = boost::posix_time::to_iso_extended_string(this_event_time);
	stringTime = "'"+stringTime+"'";
	struct tm tm;
	memset(&tm, 0, sizeof(tm));
	tm.tm_isdst = -1;
	
	strptime(stringTime.c_str(), "%Y-%m-%dT%H:%M:%S%Z", &tm);
	int d = tm.tm_mday,
            m = tm.tm_mon + 1,
            y = tm.tm_year + 1900,
            h = tm.tm_hour,
            M = tm.tm_min,
            s = tm.tm_sec;
        
	std::cout<<"y "<<y<<"   m "<<m<<"   d "<<d<<"   h "<<h<<"   M "<<M<<"   s "<<s<<std::endl;
	//time_t tstart = mktime(&tm);
	
	//std::cout<<"TPC event time"<<event_t<<"   event time (s) "<<evt_time_sec<<"   event time (ns) "<<evt_time_nsec<<std::endl;
	
	//At this stage, find out the CRT files that may be time co-incident
	
	//string filename = "hello.txt";
	//string extension = boost::filesystem::extension(filename);
	std::string teststring = "2017-06-29T16:09:57+00:00";
	std::ostringstream dim;
	dim<<"file_format "<<"crt-binaryraw"<<" and file_type "<<"data"<<" and start_time < "<<stringTime<<" and end_time > "<<stringTime;
	std::cout<<"dim = "<<dim.str()<<std::endl;
	std::vector< std::string > crtfiles = tIFDH->translateConstraints(dim.str());
	
	std::string testfn = "ProdRun20170604_001007-crt04.1.crtdaq";
	std::ostringstream dim1;
	dim1<<"file_format "<<"artroot"<<" and ischildof: (file_name "<<testfn<<")"<<std::endl;
	std::cout<<"dim1 = "<<dim1.str()<<std::endl;
	std::vector< std::string > crtrootfiles = tIFDH->translateConstraints(dim1.str());
	std::cout<<crtrootfiles[0]<<std::endl;
	
	std::cout<<"# of legitimate CRT files "<<crtrootfiles.size()<<std::endl;
	for(size_t i=0; i<crtfiles.size(); i++)
	{
		std::cout<<crtfiles[i]<<std::endl;
		std::string schema = "gsiftp";
		std::vector< std::string > crtdaqFile_url = tIFDH->locateFile(crtfiles[i], schema);
		std::cout<<"gsiftp: "<<crtdaqFile_url[0]<<std::endl;
		std::cout<<"ifdh_data_dir "<<getenv("IFDH_DATA_DIR")<<std::endl;
		//TFile* f = TFile::Open(crtdaqFile_url[0],"");
		//std::cout<<"url 1 "<<crtdaqFile_url[0]<<"\n"<<"url 2 "<<crtdaqFile_url[1]<<std::endl;
		//std::string url_1 = crtdaqFile_url[0];
		//std::string url_2 = crtdaqFile_url[1];
		
		//std::string fetchedcrtdaqfile = (tIFDH->fetchInput(crtdaqFile_url[0]));
		std::cout<<"i "<<i<<std::endl;
		//std::cout<<"---> "<<fetchedcrtdaqfile<<std::endl;
		
		//url_1 = std::regex_replace(url_1, std::regex(R"(\([^()]*\))"), "/");
		//url_1+=	crtfiles[i];
		//url_2 = url_2 + "/" + crtfiles[i];
		//std::cout<<"Modified URLs:\n "<<url_1<<"\n"<<url_2<<std::endl;
		
		//std::string url_1 = crtdaqFile_url[0];
		//auto begin = url_1.find_first_of("(");
		//auto end   = url_1.find_last_of(")");
		//if (std::string::npos!=begin && std::string::npos!=end && begin <= end)
		//url_1.erase(begin, end-begin);
		
		//url_1.erase(url_1.begin()+url_1.find("("), url_1.begin()+url_1.find(")"));
		//std::string str = "At(Robot,Room3)";
		//str = std::regex_replace(str, std::regex("[^(]*)\\([^)]*\\)(.*)"),"/");
		//std::cout<<"modilfied url_1 "<<str<<std::endl;
		//url_1=url_1+"/"+crtfiles[i];
		//std::string fileAddress = url_1;
		//std::string fileAddress = "ftp://"+crtdaqFile_url[1]+"/"+crtfiles[i];
		//std::cout<<"full address: "<<fileAddress<<std::endl;
		//bernfebdaq_ReadEvents(fetchedcrtdaqfile,"output.root",100);
	}
	std::cout<<"next"<<std::endl;
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
	std::cout<<"2 crt::CRTMerger::reconfigure"<<std::endl;
	fTag = {pset.get<std::string>("InputTagName","crtdaq")};
	fTimeWindow = pset.get<unsigned>("TimeWindow",5000000);
	fCRTFile = pset.get< std::vector < std::string > >("InputFilenames");
}
DEFINE_ART_MODULE(crt::CRTMerger)
