////////////////////////////////////////////////////////////////////////
// Class:       CRTCalibration
// Module Type: analyzer
// File:        CRTCalibration_module.cc
// Description: Module to obtain pedestals and gains from SiPMs (and do tests so far)
// Generated at Fri Nov 25 16:03:09 2016 by dalorga using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
//#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "bernfebdaq-core/Overlays/BernZMQFragment.hh"
#include "artdaq-core/Data/Fragments.hh"

#include "art/Framework/Services/Optional/TFileService.h"

#include "uboone/CRT/CRTAuxFunctions.hh"

#include "TTree.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TDatime.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <fstream>


namespace bernfebdaq {
  class CRTCalibration;
}

class bernfebdaq::CRTCalibration : public art::EDAnalyzer {
public:
  explicit CRTCalibration(fhicl::ParameterSet const & pset);
  
  // Required functions.
  void analyze(art::Event const & evt) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  //void reconfigure(fhicl::ParameterSet const & p) override;

  //int getFEBN(uint64_t febid);
  //void FillPos(std::string filePos,   std::map <int, std::vector<double> >& sensor_pos); //key = FEB*100+ch

private:

  // Declare member data here.

  //this is the magic stuff that sets up the tree.
  art::ServiceHandle<art::TFileService> tfs;
  
  std::string  raw_data_label_;
  //unsigned int max_time_diff_;
  std::string  SiPMpositions_;
  
  TTree*       my_tree_;
  
  Double_t Hist_min;
  Double_t Hist_max;
  Int_t Hist_bin;

  std::map <int, TH1F*> histos; //key = FEB*100+ch
  //  std::map <int, TH1F*>::iterator it;
  TH1F* h0;
  TH2F* ht;
  //TH2F* htdis;
  TH1F* hpedpos;
  TH1F* hpedpos_err;
  int key;

};


bernfebdaq::CRTCalibration::CRTCalibration(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset),
    raw_data_label_(pset.get<std::string>("raw_data_label")),
    SiPMpositions_(pset.get<std::string>("CRTpositions_file"))//,


{
  //std::cout<<"read parameters"<<std::endl;
  //std::cout<<"CRTpositionsfile: "<<SiPMpositions_<<std::endl;
  // getchar();
}

void bernfebdaq::CRTCalibration::analyze(art::Event const & evt)
{

  //time_t evtTime = evt.time().value();
  //   std::cout<<"time: "<<evtTime<<std::endl;

  // look for raw BernFEB data  
  art::Handle< std::vector<artdaq::Fragment> > rawHandle;
  evt.getByLabel(raw_data_label_, "BernZMQ", rawHandle);
  
  //check to make sure the data we asked for is valid
  if(!rawHandle.isValid()){
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
              << ", event " << evt.event() << " has zero"
              << " BernFEB fragments " << " in module " << raw_data_label_ << std::endl;
    std::cout << std::endl;
    return;
  }
  
  //get better access to the data
  std::vector<artdaq::Fragment> const& rawFragments(*rawHandle);

  //loop over the raw data fragments
  //There should be one fragment per FEB in each event.
  for(auto const& frag : rawFragments){
    
    //overlay this so it's in the "BernFragment" format. Same data!
    BernZMQFragment bfrag(frag);
    
    //Grab the metadata.
    //See bernfebdaq-core/bernfebdaq-core/Overlays/BernFEBFragment.hh
    auto bfrag_metadata = bfrag.metadata();
 
    size_t   nevents    = bfrag_metadata->n_events();   //number of BernFEBEvents in this packet

    auto FEB_ID =  bfrag_metadata->feb_id();     //mac addresss of this packet
    auto IDN = crt::auxfunctions::getFEBN(FEB_ID); //FEB ID   

    for(size_t i_e=0; i_e<nevents; ++i_e){
      BernZMQEvent const* this_event = bfrag.eventdata(i_e); //get the single hit/event

      //      auto time_ts0 = this_event->Time_TS0();         //grab the event time                                                                       
      //      double corrected_time = GetCorrectedTime(time_ts0,*bfrag_metadata);


      // double tdis = time_ts0 - corrected_time;
      //htdis->Fill(tdis);      

      for(size_t i_chan=0; i_chan<32; ++i_chan){
  
  if((this_event->adc[i_chan])>32 && (this_event->adc[i_chan])<4080){//remove fake data
    
    key = IDN*100 + i_chan;  
    
    auto it = histos.find(key);
    TH1F* h = (*it).second;
    h->Fill(this_event->adc[i_chan]);
    
    h0->Fill(this_event->adc[i_chan]);
    
    //ht->Fill(evtTime, this_event->adc[i_chan]);
    
  }
      }
    }
  }
}//end analyzer

void bernfebdaq::CRTCalibration::beginJob()
{
  // Implementation of optional member function here.
  
  my_tree_ = tfs->make<TTree>("my_tree","CRT Analysis Tree");

  Hist_min = 1;
  Hist_max = 4096;
  Hist_bin = 4096;

  h0 = tfs->make<TH1F>("h0","h0", Hist_bin, Hist_min, Hist_max);

  /*
  //time
  TDatime tmin("2016-12-02 02:00:00");
  TDatime tmax("2016-12-02 14:00:00"); 
  
  Double_t   Time_Min =  tmin.Convert();
  Double_t   Time_Max =  tmax.Convert(); 
  Int_t   Time_bin =  tmax.Convert()-tmin.Convert();

  ht = tfs->make<TH2F>("ht","ht",Time_bin, Time_Min, Time_Max, Hist_bin, Hist_min, Hist_max);
  ht->GetXaxis()->SetTitle("Date");
  ht->GetYaxis()->SetTitle("ADC");
  */

  //  htdis = tfs->make<TH2F>("ht","ht",Time_bin, Time_Min, Time_Max, 1000, -500, 500);
  //htdis->GetXaxis()->SetTitle("Time event - corrected (ns)");
  //htdis->GetYaxis()->SetTitle("ADC");


  //  ht = tfs->make<TH2F>("ht","ht", Hist_bin, Hist_min, Hist_max);

  char *namehist = new char[128];

  for (int FEB = 11; FEB<62; FEB++){

    if(FEB==13){FEB++;}
    if(FEB==25){FEB++;}

    for(int ch=0; ch<32; ch++){

      key = FEB*100 + ch;

      sprintf (namehist, "FEB_%d_SiPM_%d", FEB, ch);

      TH1F* h = tfs->make<TH1F>(namehist,namehist, Hist_bin, Hist_min, Hist_max);
      h->GetXaxis()->SetTitle("ADC Channel");
      h->GetYaxis()->SetTitle("Entries/bin");
      histos.emplace(key,h);
    }
  }
  
}

void bernfebdaq::CRTCalibration::endJob(){
  
  int bin;
  hpedpos = tfs->make<TH1F>("hpedpos","CRT_pedestal_positions", 5200, 1000, 6200);
  hpedpos->GetXaxis()->SetTitle("IDN (FEB*100+SiPMch)");
  hpedpos->GetYaxis()->SetTitle("ADC Counts");

  hpedpos_err = tfs->make<TH1F>("hpedpos_err","CRT_pedestal_position_error", 5200, 1000, 6200);
  hpedpos_err->GetXaxis()->SetTitle("IDN (FEB*100+SiPM)");
  hpedpos_err->GetYaxis()->SetTitle("ADC Counts");

  std::map<int,std::pair<double,double> > values;  
  
  for (int FEB = 11; FEB<62; FEB++){
    if(FEB==13){FEB++;}
    if(FEB==25){FEB++;}
    for(int ch=0; ch<32; ch++){

      key = FEB*100 + ch;
      auto it = histos.find(key);
      TH1F* h = (*it).second;

      //auto entries = h->GetEntries();

      double maxat= h->GetMaximumBin();
      //double maxval = h->GetBinContent(h->GetMaximumBin());

      TF1 *fitped = new TF1("fitped","gaus", maxat-20, maxat+20);
      // TF1 *fitped = new TF1("fitped","gaus", 400, 2000);
      h->Fit(fitped, "R+");
      double ped = fitped->GetParameter(1);
      double ped_err= fitped->GetParError(1);

      //      std::cout <<"Entries for FEB "<<FEB<<" SiPM "<<ch<<" : "<<entries<<std::endl;
      //      std::cout <<"Max value "<<maxval<<" at "<<maxat<<std::endl;
      //      std::cout <<"Pedestal value "<<ped<<" with error "<<ped_err<<std::endl;
      //      std::cout <<key<<"   " <<ped<<" "<<ped_err<<std::endl;
      std::pair<double,double> chped;
      chped = std::make_pair(ped, ped_err);
      values[key]=chped;

      bin = key - 1000;
      hpedpos->SetBinContent(bin,ped);
      hpedpos->SetBinError(bin,ped_err);
      hpedpos_err->SetBinContent(bin,ped_err);
      //hpedpos->Fill(key, ped);
      // hpedpos->SetBinError(key, ped_err);




      //calcular gain
      //      TH1F *hnew = (TH1F*)h->Clone("hnew");
      //TF1 *fitmean = new TF1("fitmean","gaus", maxat-20, maxat+20);
      //h->Fit(fitmean, "R+");
      //double mean = fitmean->GetParameter(1);
      //double mean_err= fitmean->GetParError(1);



    }
  }

  hpedpos->SetMarkerStyle(2);
  hpedpos_err->SetMarkerStyle(2);


  for(auto itA = begin(values); itA != end(values); ++itA){
    int key = (*itA).first;
    std::pair<double,double> val = (*itA).second;

    std::cout<<key<<"   "<<val.first<<"   "<<val.second<<std::endl;
   
  
  }

  
  
}


DEFINE_ART_MODULE(bernfebdaq::CRTCalibration)
