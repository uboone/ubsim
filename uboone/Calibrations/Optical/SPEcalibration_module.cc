////////////////////////////////////////////////////////////////////////
// Class:       analyzer
// Module Type: SPEcalibration
// File:        SPEcalibration_module.cc
//
// Generated at Tue Oct 13 13:14:22 2015 by Jarrett Moon using artmod
// from cetpkgsupport v1_08_06.
////////////////////////////////////////////////////////////////////////

// ART libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "art/Framework/Services/Optional/TFileService.h"

// Supporting libraries
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "LogicPulseFinder.h"

// C++ libraries
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

// ROOT
#include "TTree.h"
#include "TH1D.h"

// LArSoft 
#include "SimpleTypesAndConstants/RawTypes.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "Utilities/TimeService.h"

//Optical Channel Maps
#include "uboone/Geometry/UBOpChannelTypes.h"
#include "uboone/Geometry/UBOpReadoutMap.h"

//RawDigits
#include "RawData/raw.h"
#include "RawData/RawDigit.h"
#include "RawData/OpDetWaveform.h"

// Pulse finding
#include "uboone/OpticalDetectorAna/OpticalSubEvents/cfdiscriminator_algo/cfdiscriminator.hh"



class SPEcalibration;

class SPEcalibration : public art::EDAnalyzer {
public:
  explicit SPEcalibration(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SPEcalibration(SPEcalibration const &) = delete;
  SPEcalibration(SPEcalibration &&) = delete;
  SPEcalibration & operator = (SPEcalibration const &) = delete;
  SPEcalibration & operator = (SPEcalibration &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;


private:

  int irun;
  int isubrun;
  int ievent;
  int opchannel;
  int readout_ch;
  double charge;
  double maxamp;
  double ave_baseline;
  double rms_baseline;
  double ave_baseline2;
  double rms_baseline2;
  int slot;
  int ch;

  std::vector<unsigned short> flasher_adcs;
  std::vector<unsigned short> adcs;
  std::vector<int>* event_nchfires;
  std::vector<double>* event_lastcosmicwin;
  std::vector<double>* event_lastsatcosmicwin;
  int nsamples = 0;
  double event_maxamp;

  unsigned short LOGIC_CH;
  unsigned short LOGIC_SLOT;
  unsigned short LOGIC_THRESHOLD;
  unsigned short SPE_SLOT;
  unsigned int WINDOW_MINSIZE;
  int COSMIC_THRESHOLD;
  unsigned short baseline_start;
  unsigned short baseline_end;
  unsigned short pulse_start;
  unsigned short pulse_end;
  unsigned short baseline2_start;
  unsigned short baseline2_end;
  double cfd_threshold;
  unsigned short cfd_delay;
  unsigned short cfd_deadtime;
  unsigned short cfd_width;
  bool SAVE_AVE_SPE;
  double AVESPE_RMSLIMIT;
  std::string OpDataModule;

  TTree* fTouttree;
  TTree* fTpulsetree;
  TTree* fTevent;
  TH1D** hSPE_ave;
  //TH1D* hSPE_rms;
  int hspe_nfills[32];

  // LED Pulser Mode
  void anayzeLEDPulserMode( const art::Event& evt );

  // Pulse Finding Mode
  void analyzePulseFindingMode( const art::Event& evt );

};


SPEcalibration::SPEcalibration(fhicl::ParameterSet const& p)
  :EDAnalyzer(p)
{

  // Setup Analyzer File
  art::ServiceHandle<art::TFileService> out_file;

  // Pulse Tree
  fTouttree = out_file->make<TTree>("outtree","Analyzed Flasher Run Data");
  fTouttree->Branch("run",        &irun,        "run/I");
  fTouttree->Branch("subrun",        &isubrun,        "subrun/I");
  fTouttree->Branch("event",        &ievent,        "event/I");
  fTouttree->Branch("opchannel",    &opchannel,     "opchannel/I");
  fTouttree->Branch("charge",       &charge ,       "charge/D");
  fTouttree->Branch("maxamp",       &maxamp,        "maxamp/D");
  fTouttree->Branch("baseline",     &ave_baseline,  "baseline/D");
  fTouttree->Branch("baselinerms",  &rms_baseline,  "baselinerms/D");
  fTouttree->Branch("baseline2",    &ave_baseline2, "baseline2/D");
  fTouttree->Branch("baselinerms2", &rms_baseline2, "baselinerms2/D");

  fTpulsetree = out_file->make<TTree>("pulsetree","Analyzed Pulses");
  fTpulsetree->Branch("run",        &irun,        "run/I");
  fTpulsetree->Branch("subrun",        &isubrun,        "subrun/I");
  fTpulsetree->Branch("event",        &ievent,        "event/I");
  fTpulsetree->Branch("opchannel",    &opchannel,     "opchannel/I");
  fTpulsetree->Branch("charge",       &charge ,       "charge/D");
  fTpulsetree->Branch("maxamp",       &maxamp,        "maxamp/D");
  fTpulsetree->Branch("baseline",     &ave_baseline,  "baseline/D");
  fTpulsetree->Branch("baselinerms",  &rms_baseline,  "baselinerms/D");
  fTpulsetree->Branch("baseline2",    &ave_baseline2, "baseline2/D");
  fTpulsetree->Branch("baselinerms2", &rms_baseline2, "baselinerms2/D");
  fTpulsetree->Branch("chmaxamp",     &event_maxamp,  "chmaxamp/D" );

  // Data for each Beam Window
  fTevent = out_file->make<TTree>("eventtree", "Data about the event" );
  event_nchfires = new std::vector<int>;
  event_lastcosmicwin = new std::vector<double>;
  event_lastsatcosmicwin = new std::vector<double>;
  fTevent->Branch("run",        &irun,        "run/I");
  fTevent->Branch("subrun",        &isubrun,        "subrun/I");
  fTevent->Branch("event",    &ievent,       "event/I");
  fTevent->Branch("nsamples", &nsamples,     "nsamples/I" );
  fTevent->Branch("chmaxamp", &event_maxamp, "chmaxamp/D");
  fTevent->Branch("nchfires", "vector<int>", &event_nchfires );
  fTevent->Branch("dt_lastcosmicwin_usec",    "vector<double>", &event_lastcosmicwin );
  fTevent->Branch("dt_lastsatcosmicwin_usec", "vector<double>", &event_lastsatcosmicwin );

  // Tree for PMT slam events (planned)
  // fTslam = out_file->make<TTree>("slamtree", "Data about fully staurated PMT events" );
  // slam_chtstamp = new std::vector<double>;
  // fTslam->Branch("event",     &ievent,          "event/I");
  // fTslam->Branch( "chtstamp", "vector<double>", &slam_chtstamp );

  // Setup Analyzer
  LOGIC_CH         = p.get<unsigned short>( "LogicChannel", 39 );
  LOGIC_SLOT       = p.get<unsigned short>( "LogicSlot", 6 );
  LOGIC_THRESHOLD  = p.get<unsigned short>( "LogicThreshold", 200 );
  COSMIC_THRESHOLD  = p.get<int>( "LogicThreshold", 2348 );  // ~15 pe (with pe=20 ADC counts)
  SPE_SLOT         = p.get<unsigned short>( "WaveformSlot", 5 );
  WINDOW_MINSIZE   = p.get<unsigned int>( "WindowMinSize", 500 );
  baseline_start   = p.get<unsigned short>( "LeadingBaselineStart", 0 );   // first tick after logic pulse start that defines baseline window
  baseline_end     = p.get<unsigned short>( "LeadingBaselineEnd", 180 );   // last tick after lofic pulses start that defines baseline window
  baseline2_start  = p.get<unsigned short>( "TrailingBaselineStart", 210 );   // first tick after logic pulse start that defines baseline window
  baseline2_end    = p.get<unsigned short>( "TrailingBaselineEnd", 260  );   // last tick after lofic pulses start that defines baseline window
  pulse_start      = p.get<unsigned short>( "PulseStart", 180 );    // first tick after logic pulse start that defines pulse window
  pulse_end        = p.get<unsigned short>( "PulseEnd", 210 );      // first tick after logic pulse start that defines pulse window
  OpDataModule     = p.get<std::string>( "OpDataModule", "pmtreadout" );
  SAVE_AVE_SPE     = p.get<bool>( "SaveAverageSPE", true );
  //USE_PULSE_FINDER = p.get<bool>( "UsePulseFinder", false );
  cfd_threshold    = p.get<double>( "CFDThreshold", 10.0 );
  cfd_delay        = p.get<unsigned short>( "CFDDelay", 4 );
  cfd_deadtime     = p.get<unsigned short>( "CFDDeadtime", 20 );
  cfd_width        = p.get<unsigned short>( "CFDWidth", 20 );
  AVESPE_RMSLIMIT  = p.get<double>( "AveSPEBaselineRMScut", 2.0 );

  // Ave SPE waveform
  hSPE_ave = new TH1D*[32];
  for (int i=0; i<32; i++) {
    hspe_nfills[i] = 0;
    char hname[32];
    sprintf( hname, "hSPE_ave_femch%02d", i );
    hSPE_ave[i] = out_file->make<TH1D>(hname, "Average SPE waveform;ticks from disc. fire-delay-deadtime", (int)cfd_width, 0.0, (double)cfd_width );
    //hSPE_rms = out_file->make<TH1D>("hSPE_rms", "Average SPE waveform;ticks from disc. fire-delay-deadtime", (int)3*cfd_deadtime, 0.0, (double)3*cfd_deadtime );
  }
  
  
}

void SPEcalibration::analyze(const art::Event& evt)
{
  // Run both routines
  irun    = (int)evt.run();
  isubrun = (int)evt.subRun();
  ievent    = (int)evt.event();
  //if ( !USE_PULSE_FINDER )
  anayzeLEDPulserMode( evt );    
  //else
  analyzePulseFindingMode( evt );

}


void SPEcalibration::anayzeLEDPulserMode(const art::Event& evt)
{
  // initialize data handles and services
  art::ServiceHandle<geo::UBOpReadoutMap> ub_PMT_channel_map;
  art::Handle< std::vector< raw::OpDetWaveform > > LogicHandle;
  art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;

  evt.getByLabel( OpDataModule, "OpdetBeamHighGain", wfHandle);
  evt.getByLabel( OpDataModule, "UnspecifiedLogic" , LogicHandle);
  std::vector<raw::OpDetWaveform> const& opwfms(*wfHandle);
  std::vector<raw::OpDetWaveform> const& logwfms(*LogicHandle);
  

  for(auto &wfm : logwfms)  {

    readout_ch = wfm.ChannelNumber();
    unsigned int c,s,f;
    ub_PMT_channel_map->GetCrateSlotFEMChFromReadoutChannel(readout_ch, c, s, f);
    slot = (int)s;
    ch = (int)f;

    if(slot == (int)LOGIC_SLOT && ch == (int)LOGIC_CH) {

      flasher_adcs.clear();
      for(auto &adc : wfm)
	flasher_adcs.push_back( (short)adc );
      break;
    }

  }

  LogicPulseFinder<unsigned short> LPF;
  std::vector<unsigned short> ttlpulses;
  ttlpulses = LPF.Get_TTL_Starts(flasher_adcs, LOGIC_THRESHOLD);


  for(auto &wfm : opwfms)  {

    readout_ch = wfm.ChannelNumber();
    unsigned int c,s,f;
    ub_PMT_channel_map->GetCrateSlotFEMChFromReadoutChannel(readout_ch, c, s, f);
    opchannel = (int)f;

    adcs.clear();
    for(auto &adc : wfm)
      adcs.push_back( (short)adc );


    //These four vectors will be filled in the next loop with baseline information, one baseline preceding and one                                          
    //following the pulse region                                                                                                                            

    std::vector<double> baseline_avgs;
    std::vector<double> baseline_vars;
    std::vector<double> baseline_avgs2;
    std::vector<double> baseline_vars2;

    for ( unsigned int ipulse=0; ipulse<ttlpulses.size(); ipulse++ ) {

      double x = 0.;
      double xx = 0.;
      double ticks = 0;
      for ( unsigned int tdc=ttlpulses.at(ipulse)+baseline_start; tdc<ttlpulses.at(ipulse)+baseline_end; tdc++ ) {
	if ( tdc>=adcs.size() )
	  break;
	x += (double)adcs.at( tdc );
	xx += ((double)adcs.at( tdc ))*((double)adcs.at(tdc));
	ticks += 1.0;
      }

      baseline_vars.push_back(sqrt( xx/ticks - (x/ticks)*(x/ticks) ));
      baseline_avgs.push_back(x / ticks);

      double x2 = 0.;
      double xx2 = 0.;
      double ticks2 = 0;

      for( unsigned int tdc=ttlpulses.at(ipulse) + baseline2_start; tdc<ttlpulses.at(ipulse)+baseline2_end; tdc++) {

	if (tdc>=adcs.size())
	  break;
	x2 += (double)adcs.at(tdc);
	xx2 += ((double)adcs.at(tdc))*((double)adcs.at(tdc));
	ticks2 += 1.0;

      }

      baseline_vars2.push_back(sqrt(xx2/ticks2 - (x2/ticks2)*(x2/ticks2)));
      baseline_avgs2.push_back(x2/ticks2);

    }

    //-----------------------------------------------------------------------------------------------------------                                           
    // Now that we have baseline information, we proceed to deal with the pulse region. Following we integrate                                              
    // and find the amplitude in the pulse region and use the baseline information to correct it                                                            

    for (unsigned int ipulse=0; ipulse<ttlpulses.size(); ipulse++ ) {
      unsigned int start = ttlpulses.at(ipulse)+pulse_start;
      unsigned int end   = ttlpulses.at(ipulse)+pulse_end;

      if ( end > adcs.size() )
	end = adcs.size();

      double ch_charge = 0.0;
      maxamp = -1;

      for (unsigned int i=start; i<end; i++ ) {
	
	ch_charge += adcs.at(i)-baseline_avgs[ipulse];
	
	if ( maxamp<adcs.at(i)-baseline_avgs[ipulse] )
	  maxamp = adcs.at(i)-baseline_avgs[ipulse];

      }

      // ------------------------------------------------------------------------------------------------------                
      // Fill the output tree with all the goodies we have calculated here                                                                                
      ave_baseline  = baseline_avgs[ipulse];                                                                                                              
      ave_baseline2 = baseline_avgs2[ipulse];                                                                                                             
      rms_baseline  = baseline_vars[ipulse];                                                                                                              
      rms_baseline2 = baseline_vars2[ipulse];                                                                                                             
      charge        = ch_charge;                                                                                                                          
      // ----------------------------------------------------------------                                                                                 
      fTouttree->Fill(); 
  
    }

  }

}

void SPEcalibration::analyzePulseFindingMode(const art::Event& evt)
{
  // initialize data handles and services
  art::ServiceHandle<geo::UBOpReadoutMap> ub_PMT_channel_map;
  art::ServiceHandle<util::TimeService> ts;
  art::Handle< std::vector< raw::OpDetWaveform > > LogicHandle;
  art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;

  evt.getByLabel( OpDataModule, "OpdetBeamHighGain", wfHandle);
  std::vector<raw::OpDetWaveform> const& opwfms(*wfHandle);

  event_maxamp = 0.0;
  event_nchfires->clear();
  event_nchfires->resize( 32, 0.0 );
  event_lastcosmicwin->clear();
  event_lastcosmicwin->resize( 32, 1.0 );
  event_lastsatcosmicwin->clear();
  event_lastsatcosmicwin->resize( 32, 1.0 );
  nsamples = 0;

  // get trigger time
  double trig_timestamp = ts->TriggerTime();
  
  // get channel maxamp and other event related info
  for(auto &wfm : opwfms)  {
    readout_ch = wfm.ChannelNumber();
    if ( readout_ch%100>=32 )
      continue;

    unsigned int c,s,f;
    ub_PMT_channel_map->GetCrateSlotFEMChFromReadoutChannel(readout_ch, c, s, f);
    slot = (int)s;
    ch = (int)f;
    if ( slot!=(int)SPE_SLOT )
      continue;
    if ( WINDOW_MINSIZE < wfm.size() ) {
      // beam window
      nsamples = (int)wfm.size();
    
      for (int i=0; i<(int)wfm.size(); i++) {
	if ( event_maxamp<(double)wfm.at(i) )
	  event_maxamp = (double)wfm.at(i);
      }
    }
    else {
      // cosmic window. we have a cosmic threshold (to prevent considering late or smaller light pulses)
      double dt_usec = wfm.TimeStamp()-trig_timestamp; // usec
      auto maxelem = std::max_element( wfm.begin(), wfm.end() );
      //std::cout << " " << ch << ", " << wfm.size() << ": " << dt_usec << " " << wfm.TimeStamp() << " " << trig_timestamp << " " << *maxelem << std::endl;
      if ( dt_usec>0 )
	continue;

      if ( (int)(*maxelem)>COSMIC_THRESHOLD ) {
	if ( event_lastcosmicwin->at(ch) > 0.0 || std::fabs(dt_usec)<std::fabs(event_lastcosmicwin->at(ch)) ) {
	  event_lastcosmicwin->at(ch) = dt_usec;
	}
	if ( (int)(*maxelem)>4090 ) {
	  // saturated channel
	  if ( event_lastsatcosmicwin->at(ch)>0.0 || (std::fabs(dt_usec)<std::fabs(event_lastsatcosmicwin->at(ch)) ) ) {
	    event_lastsatcosmicwin->at(ch) = dt_usec;
	  }
	}
      }//end of if cosmic threshold met
    }
  }
  
  
  for(auto &wfm : opwfms)  {
    readout_ch = wfm.ChannelNumber();

    if ( WINDOW_MINSIZE > wfm.size() )
      continue;
    nsamples = (int)wfm.size();

    unsigned int c,s,f;
    ub_PMT_channel_map->GetCrateSlotFEMChFromReadoutChannel(readout_ch, c, s, f);
    slot = (int)s;
    ch = (int)f;
    if ( slot!=(int)SPE_SLOT )
      continue;

    // convert waveform
    std::vector<double> data;
    for (int i=0; i<(int)wfm.size(); i++)
      data.push_back( (double)wfm.at(i) );

    // get pulse positions
    std::vector< int > t_fire;
    std::vector< int > amp_fire;
    std::vector< int > maxt_fire;
    std::vector< int > diff_fire;
    cpysubevent::runCFdiscriminatorCPP( t_fire, amp_fire, maxt_fire, diff_fire,
					data.data(), (int)cfd_delay, (int)cfd_threshold, (int)cfd_deadtime, (int)cfd_width, (int)data.size() );

    // we loop through the pulses.  
    if ( ch<32 )
      event_nchfires->at( ch ) = (int)t_fire.size();
    for ( int ifire=0; ifire<(int)t_fire.size(); ifire++ ) {

      int tick = (int)t_fire.at(ifire) - (int)cfd_delay;
      if ( tick<(int)2*cfd_deadtime || tick>(int)wfm.size()-(int)2*cfd_deadtime 
	   || ( ifire>0 && tick < t_fire.at(ifire-1)+(int)2*cfd_deadtime ) ) {
	// only use pulses 2 times the deadtime from the ends and 2 deadtimes away from the previous pulse
	continue;
      }
      
      // ok we like this one. baseline time. also save maximum amplitude of channel
      ave_baseline = 0.;
      rms_baseline = 0.;
      ave_baseline2 = 0.;
      rms_baseline2  = 0.;
      for (int i=0; i<(int)cfd_width; i++) {
	ave_baseline += (double)wfm.at( tick-cfd_width+i );
	rms_baseline  += (double)wfm.at( tick-cfd_width+i )*wfm.at( tick-cfd_width+i );
	ave_baseline2  += (double)wfm.at( tick+cfd_width+i );
	rms_baseline2 += (double)wfm.at( tick+cfd_width+i )*wfm.at( tick+cfd_width+i );
      }
      ave_baseline /= (double)cfd_width;
      ave_baseline2  /= (double)cfd_width;
      rms_baseline /= (double)cfd_width;
      rms_baseline2 /= (double)cfd_width;
      rms_baseline = sqrt( rms_baseline - ave_baseline*ave_baseline );
      rms_baseline2 = sqrt( rms_baseline2 - ave_baseline2*ave_baseline2 );

      // calculate integral
      maxamp = 0.0;
      charge = 0.0;
      double baseline_slope = (ave_baseline2-ave_baseline)/(double)cfd_width;
      for (int i=0; i<(int)cfd_width; i++) {
	double adc = wfm.at( tick + i ) - (baseline_slope*(i) + ave_baseline);
	charge += adc;
	if ( adc>maxamp )
	  maxamp = adc;
	if ( ch<32 && rms_baseline<AVESPE_RMSLIMIT && rms_baseline2<AVESPE_RMSLIMIT ) {
	  hSPE_ave[ch]->SetBinContent( i+1, (hSPE_ave[ch]->GetBinContent( i+1 )*hspe_nfills[ch] + adc)/double(hspe_nfills[ch]+1) );
	}
      }
      if ( ch<32 && rms_baseline<AVESPE_RMSLIMIT && rms_baseline2<AVESPE_RMSLIMIT )
	hspe_nfills[ch]++;

      // ok done for the pulse
      opchannel = readout_ch;
      fTpulsetree->Fill();
    }
  } //end of loop over waveforms

  fTevent->Fill();
}



DEFINE_ART_MODULE(SPEcalibration)
