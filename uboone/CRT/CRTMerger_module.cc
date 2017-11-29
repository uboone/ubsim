#include "uboone/CRT/CRTMerger.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "uboone/CRT/CRTProducts/MSetCRTFrag.hh"
#include "artdaq-core/Data/Fragments.hh"
#include "gallery/Event.h"
#include "CRTBernFEBDAQCore/Overlays/BernZMQFragment.hh"
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

crt::CRTMerger::CRTMerger(const fhicl::ParameterSet& pset)//: fFileNames(pset.get<std::vector<std::string> >("InputFilenames"))
{
	std::cout<<"1 crt::CRTMerger::CRTMerger"<<std::endl;
	
	setenv("TZ", "CST+6CDT", 1);
	tzset();
	
	//setenv("IFDH_DATA_DIR","/uboone/data/users/kolahalb/MicroBooNE/",1);
	//std::cout<<"ifdh_data_dir "<<getenv("IFDH_DATA_DIR")<<std::endl;
	
	this->reconfigure(pset);
	produces< std::vector<artdaq::Fragment> >();
	fTag		= pset.get<art::InputTag> ("InputTagName");
	_debug		= pset.get<bool>		 ("debug");
	fTimeOffSet	= pset.get<std::vector< unsigned long > > ("test_t_offset");
	previouscrtrootfile = "";
	if ( ! tIFDH ) tIFDH = new ifdh_ns::ifdh;
}

crt::CRTMerger::~CRTMerger()
{
	//fIFDH->cleanup();
	std::cout<<"4 crt::CRTMerger::~CRTmerger"<<std::endl;
	
	if ( ! tIFDH )
	delete tIFDH;
}

void crt::CRTMerger::produce(art::Event& event)
{
	if (_debug)
	{
		std::cout <<"3 crt::CRTMerger::produce, NEW EVENT" << std::endl;
	}
	
	std::unique_ptr<std::vector<crt::MSetCRTFrag> > MergedCRTFragSet(new std::vector<crt::MSetCRTFrag>);
	
	std::unique_ptr<std::vector<artdaq::Fragment> > ThisFragSet(new std::vector<artdaq::Fragment>);
	std::unique_ptr<std::vector<artdaq::Fragment> > LastFragSet(new std::vector<artdaq::Fragment>);
	std::unique_ptr<std::vector<artdaq::Fragment> > NextFragSet(new std::vector<artdaq::Fragment>);
	
	std::unique_ptr< std::vector < std::vector <artdaq::Fragment> > > FragMatrix(new std::vector < std::vector <artdaq::Fragment> >);
	
	//For this event
	if(_debug)
	std::cout<<fMaxCount<<std::endl;
	
	std::vector<artdaq::Fragment>  ThisFragment;
	
	// First find the art event time stamp
	art::Timestamp evtTime = event.time();
	
	unsigned long evt_time_sec = evtTime.timeHigh();	
	unsigned long evt_time_nsec = evtTime.timeLow();
	
	// Use the information for configuring SAM query about the coincident crt-binrary-raw files
	unsigned long time_tpc = evt_time_sec*1000000000 + evt_time_nsec;
	unsigned long time_tpc1= time_tpc/1000;
	
	const char* tz = getenv("TZ");
	std::string tzs(tz);
	if (_debug)
	std::cout<<"time-zone "<<tzs<<std::endl;
	
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
	
	boost::posix_time::ptime const time_epoch(boost::gregorian::date(1970, 1, 1));
	boost::posix_time::ptime this_event_time = time_epoch + boost::posix_time::microseconds(time_tpc1);
	
	if (_debug)
	std::cout << boost::posix_time::to_iso_extended_string(this_event_time)<<std::endl;
	
	std::string stringTime = boost::posix_time::to_iso_extended_string(this_event_time);
	stringTime = "'"+stringTime+"'";
	
	struct tm tm;
	memset(&tm, 0, sizeof(tm));
	tm.tm_isdst = -1;
	
	/*
	strptime(stringTime.c_str(), "%Y-%m-%dT%H:%M:%S%Z", &tm);
	int d = tm.tm_mday,
            m = tm.tm_mon + 1,
            y = tm.tm_year + 1900,
            h = tm.tm_hour,
            M = tm.tm_min,
            s = tm.tm_sec;
        
	//std::cout<<"y "<<y<<"   m "<<m<<"   d "<<d<<"   h "<<h<<"   M "<<M<<"   s "<<s<<std::endl;
	*/
	//time_t tstart = mktime(&tm);
	
	if (_debug)
	std::cout<<"TPC event time (s) "<<evt_time_sec<<"   event time (ns) "<<evt_time_nsec<<std::endl;
	
	// At this stage, find out the crtdaq files that may be time co-incident
	std::ostringstream dim;
	dim<<"file_format "<<"crt-binaryraw"<<" and file_type "<<"data"<<" and start_time < "<<stringTime<<" and end_time > "<<stringTime;
	
	if (_debug)
	std::cout<<"dim = "<<dim.str()<<std::endl;
	
	// List those crtdaq files:
	std::vector< std::string > crtfiles = tIFDH->translateConstraints(dim.str());
	
	// Find the corresponding child artroot file
	std::vector< std::string > crtrootfile;
	for(unsigned k =0; k<crtfiles.size(); k++)
	{
		std::ostringstream dim1;
		dim1<<"file_format "<<"artroot"<<" and ischildof: (file_name "<<crtfiles[k]<<")"<<std::endl;
		
		if (_debug)
		std::cout<<"dim1 = "<<dim1.str()<<std::endl;
		
		crtrootfile = tIFDH->translateConstraints(dim1.str());
		
		if (crtrootfile.size()>0)
		{
			if (_debug)
			std::cout<<"found the child of the parent crtdaq file:"<<std::endl;
			break;
		}
	}
	if (_debug)
	std::cout<<"The child artroot file is "<<crtrootfile[0]<<std::endl;
	
	// Read off the root file by streaming over internet. Use xrootd URL
	// This is alternative to use gsiftp URL. In that approach we use the URL to ifdh::fetchInput the crt file to local directory and keep it open as long as we make use of it
	//xrootd URL
	std::string schema = "root";
	std::vector< std::string > crtrootFile_xrootd_url = tIFDH->locateFile(crtrootfile[0], schema);
	std::cout<<"xrootd URL: "<<crtrootFile_xrootd_url[0]<<std::endl;
	
	// gallery, when fed the list of the xrootd URL, internally calls TFile::Open() to open the file-list
	// In interactive mode, you have to get your proxy authenticated by issuing:
	// voms-proxy-init -noregen -rfc -voms fermilab:/fermilab/uboone/Role=Analysis
	// when you would like to launch a 'lar -c ... ... ...'
	// In batch mode, this step is automatically done
	
	gallery::Event fCRTEvent(crtrootFile_xrootd_url);
	std::cout<<"Opened the CRT root file from xrootd URL"<<std::endl;
	
	unsigned nCRT = 0;
        for(fCRTEvent.toBegin(); !fCRTEvent.atEnd(); ++fCRTEvent)
        {
		auto const& TryCRT_frags = *(fCRTEvent.getValidHandle< std::vector<artdaq::Fragment> >(fTag));
		std::vector< artdaq::Fragment > f;
		f.clear();
                for(auto iv=begin(TryCRT_frags); iv!=end(TryCRT_frags); ++iv)
                {
                        auto const& fg = *iv;
                        f.emplace_back(fg);
                }
                FragMatrix->emplace_back(f);
		nCRT++;
        }
	if (_debug)
	std::cout<<"Filled FragMatrix, nCRT "<<nCRT<<std::endl;
	
	/*
	///////////////////////////////////////////////////////////////////////////////////
	////gsiftp URL
	//std::string schema = "gsiftp";
	//std::vector< std::string > crtrootFile_url = tIFDH->locateFile(crtrootfile[0], schema);
	//std::cout<<"gsiftp URL: "<<crtrootFile_url[0]<<std::endl;
	//
	////std::cout<<"ifdh_data_dir "<<getenv("IFDH_DATA_DIR")<<std::endl;
	
	////std::string fetchedcrtrootfile = "";
	////bool fetch = false;
	
	////std::cout<<"### "<<fTestFiles.size()<<std::endl;
	
	////if (previouscrtrootfile != "")
	////{
	////	std::cout<<"previouscrtrootfile: "<<previouscrtrootfile<<std::endl;
	////	if (previouscrtrootfile==crtrootfile[0])
	////	{
	////		fetchedcrtrootfile = crtrootfile[0];
	////		std::cout<<"fetchedcrtrootfile "<<fetchedcrtrootfile<<std::endl;
	////		fetch = false;
	////	}
	////	else
	////	fetch = true;
	////}
	////else
	////std::cout<<"Need to fill fTestFiles through fetchInput"<<std::endl;
	
	////if ((fetchedcrtrootfile == "") | (fetch==true))
	////fetchedcrtrootfile = (tIFDH->fetchInput(crtrootFile_url[0]));
	
	////fTestFiles.push_back(fetchedcrtrootfile);
	//////////////////////////////////////////////////////////////////////////////////
	////std::cout<<"# of legitimate CRT files "<<crtrootfile.size()<<std::endl;
	*/
	
	unsigned int count = 0;
	unsigned int n = 0;
	//gallery::Event fCRTEvent(fFileNames);
	
	if (_debug)
	std::cout<<"Start merging attempts"<<std::endl;
	
	for(fCRTEvent.toBegin(); !fCRTEvent.atEnd(); ++fCRTEvent)
	{
		int xxx = 0;
		
		auto const& crtFragSet = *(fCRTEvent.getValidHandle< std::vector<artdaq::Fragment> >(fTag));
		n=0;
		
		for(auto iv = begin(crtFragSet); iv != end(crtFragSet); ++iv)
		{
			auto const& i_crtFrag = *iv;
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
			//std::cout<<"   TPC "<<total_TPC_time<<"   CRT_start "<<time_CRT_start<<"   CRT_end "<<time_CRT_end<<std::endl;
			
			//if (total_TPC_time < time_CRT_start)
			//break;
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
				
				artdaq::Fragment last_artdaqFrag = (count > 1) ? FragMatrix->at(count-1)[n]:i_crtFrag;
				bernfebdaq::BernZMQFragment last_bernfrag(last_artdaqFrag);
				auto last_bernfrag_metadata = last_bernfrag.metadata();
				if (_debug)
				std::cout<<"last Frag: start "<<last_bernfrag_metadata->time_start_nanosec()<<", stop: "<<last_bernfrag_metadata->time_end_nanosec()<<std::endl;
				LastFragSet->emplace_back(last_artdaqFrag);
				
				artdaq::Fragment next_artdaqFrag = (count < nCRT) ? FragMatrix->at(count+1)[n]:i_crtFrag;
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
	}// end CRT evt loop
	
	if (_debug)
	std::cout<<"Set of Fragments: ("<<LastFragSet->size()<<", "<<ThisFragSet->size()<<", "<<NextFragSet->size()<<"), MergedCRTFragSet: "<<MergedCRTFragSet->size()<<std::endl;
	
	event.put(std::move(ThisFragSet));
	
	if (_debug)
	std::cout<<"---X---"<<std::endl;
}

void crt::CRTMerger::reconfigure(fhicl::ParameterSet const & pset)
{
	std::cout<<"2 crt::CRTMerger::reconfigure"<<std::endl;
	
	fTag = {pset.get<std::string>("InputTagName","crtdaq")};
	fTimeWindow = pset.get<unsigned>("TimeWindow",5000000);
	
	//fCRTFile = pset.get< std::vector < std::string > >("InputFilenames");
	
	//if ( ! fIFDH ) fIFDH = new ifdh_ns::ifdh;
	//std::string fetchedfile(fIFDH->fetchInput(fFileNames[0]));
	//fTestFiles.push_back(fetchedfile);
	/*
	gallery::Event testGallery4ifdh(fFileNames);
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
	std::cout<<"fMaxCount "<<fMaxCount<<std::endl;
	//if ( ! fIFDH )
        //delete fIFDH;
	*/
}
DEFINE_ART_MODULE(crt::CRTMerger)
