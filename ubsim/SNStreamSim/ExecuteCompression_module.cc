////////////////////////////////////////////////////////////////////////
// Class:       ExecuteCompression
// Module Type: producer
// File:        ExecuteCompression_module.cc
//
// David Caratelli - davidc1@fnal.gov - July 13 2016
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// services etc...
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// data-products
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/Utilities/AssociationUtil.h"
//#include "lardata/ArtDataHelper/WireCreator.h"

// SN Compression
#include "ubsim/SNStreamSim/Fmwk/CompressionAlgoBase.h"
#include "ubsim/SNStreamSim/Algo/AlgorithmFactory.h"

// ROOT
#include "TVector3.h"
#include <TTree.h>
#include <TStopwatch.h>

// C++
#include <memory>
#include <iostream>
#include <utility>

class ExecuteCompression;

class ExecuteCompression : public art::EDProducer {
public:
  explicit ExecuteCompression(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ExecuteCompression(ExecuteCompression const &) = delete;
  ExecuteCompression(ExecuteCompression &&) = delete;
  ExecuteCompression & operator = (ExecuteCompression const &) = delete;
  ExecuteCompression & operator = (ExecuteCompression &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;


private:

  // debug (verbose) mode?
  bool _debug;

  // Tree
  TTree* _compress_chch_tree;
  TTree* _compress_evt_tree;
  
  // TTree variables
  double _compression;
  double _compressionU;
  double _compressionV;
  double _compressionY;
  int    _evt;
  // Tree for channel-by-channel compression
  double _ch_compression;
  int    _ch;
  int    _pl;
  // keep track of number of wires scanned per plane (to calculate compession)
  int _NplU, _NplV, _NplY;

  // timer to keep track of time-performance
  TStopwatch _evtwatch; // full event time
  TStopwatch _loopwatch;
  TStopwatch _watch;
  double _time_loop, _time_get, _time_algo, _time_calc, _time_swap;

  /// Compression Algorithm Object...performs compression
  std::unique_ptr< compress::CompressionAlgoBase > _compress_algo;

  // get the geometry
  art::ServiceHandle<geo::Geometry> _geom;

  const std::vector<std::pair< compress::tick, compress::tick> > ApplyCompression(const art::Ptr<raw::RawDigit> rawwf);
  void CalculateCompression(const std::vector<short> &beforeADCs,
			    const std::vector<std::pair< compress::tick, compress::tick> > &ranges,
			    int pl, int ch);

  void beginJob() override;
  void endJob() override;
    
};


void ExecuteCompression::beginJob()
{
  
  if (_debug) { std::cout << "begin job..." << std::endl; }

  art::ServiceHandle<art::TFileService> tfs;

  _compress_evt_tree = tfs->make<TTree>("_compress_evt_tree","SNCompression event tree");
  _compress_evt_tree->Branch("_evt",&_evt,"evt/I");
  _compress_evt_tree->Branch("_compression",&_compression,"compression/D");
  _compress_evt_tree->Branch("_compressionU",&_compressionU,"compressionU/D");
  _compress_evt_tree->Branch("_compressionV",&_compressionV,"compressionV/D");
  _compress_evt_tree->Branch("_compressionY",&_compressionY,"compressionY/D");

  _compress_chch_tree = tfs->make<TTree>("_compress_chch_tree","SNCompression channel tree");
  _compress_chch_tree->Branch("_ch_compression",&_ch_compression,"ch_compression/D");
  _compress_chch_tree->Branch("_ch",&_ch,"ch/I");
  _compress_chch_tree->Branch("_pl",&_pl,"pl/I");

  _compression = _compressionU = _compressionV = _compressionY = 0;
  _NplU = _NplV = _NplY = 0;

  return;
}

void ExecuteCompression::endJob()
{
  


}


ExecuteCompression::ExecuteCompression(fhicl::ParameterSet const & p)
  : _compress_algo(nullptr)
    // Initialize member data here.
{
  
  produces< std::vector< recob::Wire > >();

  _debug             = p.get<bool>       ("debug");
  
  if (_debug) { std::cout << "setting up default compression algo" << std::endl; }
  _compress_algo =  compress::AlgorithmFactory().MakeCompressionAlgo(p);
  
}

void ExecuteCompression::produce(art::Event & e)
{

  if (_debug)
  std::cout << "NEW EVENT" << std::endl;

  // produce OpFlash data-product to be filled within module
  std::unique_ptr< std::vector<recob::Wire> > wire_v(new std::vector<recob::Wire>);

  // load rawdigits
  art::Handle<std::vector<raw::RawDigit> > rawdigit_h;
  e.getByLabel("daq",rawdigit_h);

  // make sure rawdigits look good
  if(!rawdigit_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate RawDigit!"<<std::endl;
    throw std::exception();
  }


  _loopwatch.Start();

  for(size_t h=0; h < rawdigit_h->size(); h++){

    if (_debug) { std::cout << "new channel" << std::endl; }

    art::Ptr<raw::RawDigit> rawdigit(rawdigit_h, h);

    if (_debug) { std::cout << "applying compression..." << std::endl; }
    
    auto const& ranges = ApplyCompression(rawdigit);
    
    auto const& chan = rawdigit->Channel();

    if (_debug) { std::cout << "\t chan number : " << chan << std::endl; }
    
    lar::sparse_vector<float> wf_ROIs;
    
    //loop over new waveforms created
    for (size_t n=0; n < ranges.size(); n++){
      // prepare output waveform
      compress::tick t;
      float first_tick = (float)*(ranges[n].first);
      std::vector<float> out;
      for (t = ranges[n].first; t < ranges[n].second; t++)
	out.push_back( (float)*t - first_tick );
      
      unsigned int start_tick = ranges[n].first - _compress_algo->GetInputBegin();

      if (_debug) { std::cout << "\t adding range of size " << out.size() << std::endl; }
      wf_ROIs.add_range( start_tick, std::move(out) );
      
    }// for all saved ROIs
    
    _time_swap += _watch.RealTime();

    if (wf_ROIs.size() == 0) continue;

    recob::Wire wire(std::move(wf_ROIs), chan, _geom->View(chan) );

    //wire_v->emplace_back( recob::WireCreator(std::move(wf_ROIs),*rawdigit).move() );
    
    wire_v->emplace_back(wire);
    
  }// for all RawDigit vectors

  _compressionU /= _NplU;
  _compressionV /= _NplV;
  _compressionY /= _NplY;
  _compression  /= ( _NplU + _NplV + _NplY );
  _compress_evt_tree->Fill();
  _NplU = _NplV = _NplY = 0;
  _compressionU = _compressionV = _compressionY = 0;

  
  _time_loop += _loopwatch.RealTime();
  
  e.put(std::move(wire_v));
}

// function where compression is applied on a single wf
const std::vector<std::pair< compress::tick, compress::tick> > ExecuteCompression::ApplyCompression(const art::Ptr<raw::RawDigit> rawwf)
{
  
  //Check for empty waveforms!
  if(rawwf->ADCs().size()<1){
    std::cout << "Found 0-length waveform: Ch. " << rawwf->Channel() << std::endl;
  }//if wf size < 1

  auto const& ch   = rawwf->Channel();
  auto const& wids = _geom->ChannelToWire(ch);
  auto const& pl   = wids[0].Plane;

  // finally, apply compression:
  // *-------------------------*
  // 1) Convert tpc_data object to just the vector of shorts which make up the ADC ticks
  _watch.Start();
  const std::vector<short> ADCwaveformL = rawwf->ADCs();
  // cut size so that 3 blocks fit perfectly
  int nblocks = ADCwaveformL.size()/(3*64);
  std::vector<short>::const_iterator first = ADCwaveformL.begin();
  std::vector<short>::const_iterator last  = ADCwaveformL.begin()+(3*64*nblocks);
  std::vector<short> ADCwaveform(first,last);
  _time_get += _watch.RealTime();
  // 2) Now apply the compression algorithm. _compress_algo is an instance of CompressionAlgoBase
  _watch.Start();
  if (_debug) { std::cout << "\t applying algo-specific compression to waveform of size " << ADCwaveform.size() << " and 1st entry " << ADCwaveform.at(0) << "\t Nblocks : " << nblocks << std::endl; }
  _compress_algo->ApplyCompression(ADCwaveform,pl,ch);
  _time_algo += _watch.RealTime();
  // 3) Retrieve output ranges saved
  if (_debug) { std::cout << "\t getting output ranges" << std::endl; }
  auto const& ranges = _compress_algo->GetOutputRanges();
  if (_debug) { std::cout << "\t found " << ranges.size() << " ranges" << std::endl; }
  // 6) Study the Compression results for this channel
  _watch.Start();
  // 9) Calculate compression factor [ for now Ticks After / Ticks Before ]
  _watch.Start();
  if (_debug) { std::cout << "\t calculate compression" << std::endl; }
  CalculateCompression(ADCwaveform, ranges, pl, ch);
  _time_calc += _watch.RealTime();
  // 10) clear _InWF and _OutWF from compression algo object -> resetting algorithm for next time it is called
  _compress_algo->Reset();
  // return ranges
  return ranges;
}
 

void ExecuteCompression::CalculateCompression(const std::vector<short> &beforeADCs,
					      const std::vector<std::pair< compress::tick, compress::tick> > &ranges,
					      int pl, int ch){
  
  double inTicks = beforeADCs.size();
  double outTicks = 0;
  
  for (size_t n=0; n < ranges.size(); n++)
    outTicks += (ranges[n].second-ranges[n].first);
  
  if (pl==0){
    _compressionU += outTicks/inTicks;
    _NplU += 1;
  }
  else if (pl==1){
    _compressionV += outTicks/inTicks;
    _NplV += 1;
  }
  else if (pl==2){
    _compressionY += outTicks/inTicks;
    _NplY += 1;
  }
  else
    std::cout << "What plane? Error?" << std::endl;
  
  _ch_compression = outTicks/inTicks;
  
  _compression += outTicks/inTicks;
  
  _ch = ch;
  _pl = pl;
  _compress_chch_tree->Fill();
  
  return;
}

DEFINE_ART_MODULE(ExecuteCompression)
