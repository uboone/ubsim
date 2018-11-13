////////////////////////////////////////////////////////////////////////
// Class:       FakePhotonsAna
// Module Type: analyzer
// File:        FakePhotonsAna_module.cc
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include <iostream>
#include "lardataobj/RawData/OpDetWaveform.h"
#include <TTree.h>
#include <TFile.h>

class FakePhotonsAna;

class FakePhotonsAna : public art::EDAnalyzer {
public:
  explicit FakePhotonsAna(fhicl::ParameterSet const & p);
  virtual ~FakePhotonsAna();

  void analyze(art::Event const & evt) override;
  void endJob();

private:
  std::string _module_label;
  TTree* _tree;
  TFile* _file;
  int _threshold;
  int _store_wf;
  std::vector< std::vector<int> > _wf_v;
  std::vector< int > _max_v;
  std::vector< int > _ch_v;
  std::vector< int > _tp_v;
  std::vector< double > _ts_v;
};


FakePhotonsAna::FakePhotonsAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  _threshold = p.get<int>("Threshold",0);
  _store_wf = p.get<int>("StoreWaveform",0);
  _module_label = p.get<std::string>("ModuleLabel");
  std::string ana_fname = p.get<std::string>("AnaFile","ana.root");
  _file = TFile::Open(ana_fname.c_str(),"RECREATE");
  _tree = new TTree("tree","tree");
  if(_store_wf) _tree->Branch("wf_v",&_wf_v);
  _tree->Branch("max_v",&_max_v);
  _tree->Branch("ch_v",&_ch_v);
  _tree->Branch("tpeak_v",&_tp_v);
  _tree->Branch("tstart_v",&_ts_v);
}

FakePhotonsAna::~FakePhotonsAna() {
  // Clean up dynamic memory and other resources here.
}

void FakePhotonsAna::endJob() {
  if(_file) {
    _file->cd();
    _tree->Write();
    _file->Close();
  }
}

void FakePhotonsAna::analyze(art::Event const & evt) {

  _wf_v.clear();
  _max_v.clear();
  _ch_v.clear();
  _tp_v.clear();
  _ts_v.clear();

  art::Handle< std::vector<raw::OpDetWaveform> > wf_h;
  evt.getByLabel(_module_label,wf_h);
  if(!wf_h.isValid()) {
    std::cout << Form("Did not find any OpDetWaveform a prodcuer: %s",_module_label.c_str()) << std::endl;
    return;
  }

  _wf_v.reserve(wf_h->size());
  _max_v.reserve(wf_h->size());
  _ts_v.reserve(wf_h->size());
  _tp_v.reserve(wf_h->size());
  _ch_v.reserve(wf_h->size());

  std::vector<int> out_wf;
  for(size_t i=0; i<wf_h->size(); ++i) {
    art::Ptr<::raw::OpDetWaveform> wf_ptr(wf_h,i);
    auto const& wf = (*wf_ptr);
    int max = 0;
    int tpeak = 0;
    for(int j=0; j<(int)(wf.size()); ++j) {
      if(wf[j] <= max) continue;
      max = wf[j];
      tpeak = j;
    }
    if(max<_threshold) continue;
    if(_store_wf) {
      out_wf.resize(wf.size());
      for(size_t j=0; j<wf.size(); ++j) out_wf[j]=(int)(wf[j]);
      _wf_v.push_back(out_wf);
    }
    _max_v.push_back(max);
    _ch_v.push_back(wf_ptr->ChannelNumber());
    _tp_v.push_back(tpeak);
    _ts_v.push_back(wf_ptr->TimeStamp());
  }
  _tree->Fill();
  
}

DEFINE_ART_MODULE(FakePhotonsAna)
