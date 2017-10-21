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
#include "boost/date_time/local_time_adjustor.hpp"
#include <boost/date_time/c_local_time_adjustor.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
using namespace boost::posix_time;
crt::CRTMerger::CRTMerger(const fhicl::ParameterSet& pset)//: fFileNames(pset.get<std::vector<std::string> >("InputFilenames"))
{
	std::cout<<"crt::CRTMerger::CRTMerger"<<std::endl;
	
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
	std::cout<<"crt::CRTMerger::~CRTmerger"<<std::endl;
	
	if ( ! tIFDH )
	delete tIFDH;
}

void crt::CRTMerger::produce(art::Event& event)
{
	if (_debug)
	{
		std::cout <<"crt::CRTMerger::produce, NEW EVENT" << std::endl;
	}
	
	std::unique_ptr< std::vector <artdaq::Fragment> > MergedSet (new std::vector <artdaq::Fragment>);
	
	//For this event
	if(_debug)
	std::cout<<fMaxCount<<std::endl;
	
	std::vector<artdaq::Fragment>  ThisFragment;
	
	// First find the art event time stamp
	art::Timestamp evtTime = event.time();
	
	unsigned long evt_time_sec = evtTime.timeHigh()+fTimeOffSet[2];	
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
	typedef boost::date_time::c_local_adjustor<ptime> local_adj;
	ptime this_event_localtime = local_adj::utc_to_local(this_event_time);
	std::cout<<"Local time of event: "<< boost::posix_time::to_iso_extended_string(this_event_localtime)<<std::endl;
	
	if (_debug)
	std::cout << boost::posix_time::to_iso_extended_string(this_event_time)<<std::endl;
	
	std::string stringTime = boost::posix_time::to_iso_extended_string(this_event_localtime);
	stringTime = "'"+stringTime+"'";
	
	struct tm tm;
	memset(&tm, 0, sizeof(tm));
	tm.tm_isdst = -1;
	
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
	std::vector< std::string > tmprootfile;
	
	for(unsigned k =0; k<crtfiles.size(); k++)
	{
		std::ostringstream dim1;
		dim1<<"file_format "<<"artroot"<<" and ischildof: (file_name "<<crtfiles[k]<<")"<<std::endl;
		
		if (_debug)
		std::cout<<"dim1 = "<<dim1.str()<<std::endl;
		
		tmprootfile = tIFDH->translateConstraints(dim1.str());
		
		if (tmprootfile.size()>0)
		crtrootfile.push_back(tmprootfile[0]);
		/*
		crtrootfile = tIFDH->translateConstraints(dim1.str());
		
		if (crtrootfile.size()>0)
		{
			if (_debug)
			std::cout<<"found the child of the parent crtdaq file:"<<std::endl;
			break;
		}
		*/
	}
	std::cout<<"total: "<<crtrootfile.size()<<std::endl;
	
	for(unsigned crf_index = 0; crf_index < crtrootfile.size(); crf_index++)
	{
	if (_debug)
	std::cout<<"The child artroot file is "<<crtrootfile[crf_index]<<std::endl;
	
	// Read off the root file by streaming over internet. Use xrootd URL
	// This is alternative to use gsiftp URL. In that approach we use the URL to ifdh::fetchInput the crt file to local directory and keep it open as long as we make use of it
	//xrootd URL
	std::string schema = "root";
	std::vector< std::string > crtrootFile_xrootd_url = tIFDH->locateFile(crtrootfile[crf_index], schema);
	std::cout<<"xrootd URL: "<<crtrootFile_xrootd_url[0]<<std::endl;
	
	// gallery, when fed the list of the xrootd URL, internally calls TFile::Open() to open the file-list
	// In interactive mode, you have to get your proxy authenticated by issuing:
	// voms-proxy-init -noregen -rfc -voms fermilab:/fermilab/uboone/Role=Analysis
	// when you would like to launch a 'lar -c ... ... ...'
	// In batch mode, this step is automatically done
	
	gallery::Event fCRTEvent(crtrootFile_xrootd_url);
	std::cout<<"Opened the CRT root file from xrootd URL"<<std::endl;
	
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
	//unsigned int n = 0;
	//gallery::Event fCRTEvent(fFileNames);
	
	if (_debug)
	std::cout<<"Start merging attempts"<<std::endl;
	
	unsigned int merging = 0;
	
	for(fCRTEvent.toBegin(); !fCRTEvent.atEnd(); ++fCRTEvent)
	{
		auto const& artdaqFragSet = *(fCRTEvent.getValidHandle< std::vector<artdaq::Fragment> >(fTag));
		auto iv = begin(artdaqFragSet);
		
		auto const& i_crtFrag= *iv;
		bernfebdaq::BernZMQFragment bernfrag(i_crtFrag);
		auto i_Frag_metadata = bernfrag.metadata();
		
                unsigned long crt_bf_time_s             = i_Frag_metadata->time_start_seconds();
                unsigned long crt_bf_time_e             = i_Frag_metadata->time_end_seconds();
                unsigned long crt_bf_time_start_ns      = i_Frag_metadata->time_start_nanosec();
                unsigned long crt_bf_time_end_ns        = i_Frag_metadata->time_end_nanosec();
                unsigned long total_TPC_time            = evt_time_sec*1000000000+evt_time_nsec;
                unsigned long time_CRT_start            = crt_bf_time_s*1000000000+crt_bf_time_start_ns;
                unsigned long time_CRT_end              = crt_bf_time_e*1000000000+crt_bf_time_end_ns;
		
                if ((total_TPC_time>=time_CRT_start) && (total_TPC_time<=time_CRT_end) && (merging != 2))
                {
                        std::cout<<"merging happens"<<std::endl;
                        merging = 1;
                }
		
		if (merging != 0)
		std::cout<<count<<"     "<<merging<<"     "<<MergedSet->size()<<std::endl;
		
		if (merging == 0)
		{
			if (MergedSet->size()>0)
			{
				MergedSet->clear();
				//LastFragSet->clear();
			}
			
			for(auto i_f = begin(artdaqFragSet); i_f != end(artdaqFragSet); ++i_f)
			{
				auto const& ifrag = *i_f;
				MergedSet->push_back(ifrag);
				//LastFragSet->push_back(ifrag);
			}
		}
		else if (merging == 1)
		{
			std::cout<<"merging=1"<<std::endl;
                        for(auto i_f = begin(artdaqFragSet); i_f != end(artdaqFragSet); ++i_f)
                        {
                                auto const& ifrag = *i_f;
                                MergedSet->push_back(ifrag);
				//ThisFragSet->push_back(ifrag);
                        }
			merging = 2;
		}
		else if (merging == 2)
		{
			std::cout<<"merging=2"<<std::endl;
                        for(auto i_f = begin(artdaqFragSet); i_f != end(artdaqFragSet); ++i_f)
                        {
                                auto const& ifrag = *i_f;
                                MergedSet->push_back(ifrag);
				//NextFragSet->push_back(ifrag);
                        }
			break;
		}
		else
		{
		}
		++count;
	}
	}
	event.put(std::move(MergedSet));
	
	if (_debug)
	std::cout<<"---X---"<<std::endl;
}

void crt::CRTMerger::reconfigure(fhicl::ParameterSet const & pset)
{
	std::cout<<"crt::CRTMerger::reconfigure"<<std::endl;
	
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
