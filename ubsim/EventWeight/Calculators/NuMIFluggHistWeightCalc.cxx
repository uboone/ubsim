#include "larsim/EventWeight/Base/WeightCalcCreator.h"
#include "larsim/EventWeight/Base/WeightCalc.h"

#include <iostream>

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "nutools/RandomUtils/NuRandomService.h"

#include "CLHEP/Random/RandGaussQ.h"

#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "TFile.h"
#include "TH1F.h"

namespace evwgh {
  class NuMIFluggHistWeightCalc : public WeightCalc
  {
  public:
    NuMIFluggHistWeightCalc() = default;
    void Configure(fhicl::ParameterSet const& pset,
                   CLHEP::HepRandomEngine& engine);
    std::vector<std::vector<double> > GetWeight(art::Event & e);

  private:
    std::vector<double> fWeightArray{};
    int fNmultisims{};
    std::string fMode{};
    std::string fHorn{};
    std::string fGenieModuleLabel{};

    //         numu, numubar, nue, nuebar
    //         |  50MeV bins
    //         |  |
    double fRW[4][200];

    DECLARE_WEIGHTCALC(NuMIFluggHistWeightCalc)
  };

  void NuMIFluggHistWeightCalc::Configure(fhicl::ParameterSet const& p,
                                     CLHEP::HepRandomEngine& engine)
  {
    //global config
    fGenieModuleLabel= p.get< std::string > ("genie_module_label");

    fhicl::ParameterSet const &pset=p.get<fhicl::ParameterSet> (GetName());
    //calc config
    fNmultisims = pset.get<int>("number_of_multisims");
    fMode       = pset.get<std::string>("mode");
    fHorn       = pset.get<std::string>("horn", "FHC");
    std::string dataInput = pset.get< std::string >("rw_hist_file");

    // cet::search_path sp("FW_SEARCH_PATH");
    // std::string rwfile = sp.find_file(dataInput);
    std::string rwfile = dataInput;

    std::string ntype[] = {"numu", "numubar", "nue", "nuebar"};

    TFile frw(Form("%s",rwfile.c_str()));
    for (int intyp=0;intyp<4;intyp++) {
      TH1F* hrw = dynamic_cast<TH1F*> (frw.Get(Form("ratio_%s_%s",ntype[intyp].c_str(), fHorn.c_str())));
      for (int ibin=0;ibin<200;ibin++) {
        fRW[intyp][ibin]= hrw->GetBinContent(ibin+1);
      }
    }
    frw.Close();

    fWeightArray.resize(fNmultisims);

    if (fMode.find("multisim") != std::string::npos )
      for (double& weight : fWeightArray) weight = CLHEP::RandGaussQ::shoot(&engine, 0, 1.);
    else
      for (double& weight : fWeightArray) weight = 1.;
  }

  std::vector<std::vector<double> > NuMIFluggHistWeightCalc::GetWeight(art::Event & e)
  {
    //calculate weight(s) here
    std::vector<std::vector<double> > weight;

    // * MC flux information
    art::Handle< std::vector<simb::MCFlux> > mcfluxListHandle;
    std::vector<art::Ptr<simb::MCFlux> > fluxlist;
    if (e.getByLabel(fGenieModuleLabel,mcfluxListHandle))
      art::fill_ptr_vector(fluxlist, mcfluxListHandle);
    else{ return weight;}

    // * MC truth information
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    if (e.getByLabel(fGenieModuleLabel,mctruthListHandle))
      art::fill_ptr_vector(mclist, mctruthListHandle);
    else{return weight;}


    weight.resize(mclist.size());
    for (unsigned int inu=0;inu<mclist.size();inu++) {
      weight[inu].resize(fNmultisims);

      int ntype=-9999;
      int bin=-9999;

      if ( fluxlist[inu]->fntype==14 ) ntype=0;
      else if ( fluxlist[inu]->fntype==-14 ) ntype=1;
      else if ( fluxlist[inu]->fntype==12 ) ntype=2;
      else if ( fluxlist[inu]->fntype==-12 ) ntype=3;
      else {
	      throw cet::exception(__FUNCTION__) << GetName()<<"::Unknown ntype "<<fluxlist[0]->fntype<< std::endl;
      }

      double enu=mclist[inu]->GetNeutrino().Nu().E();
      bin=enu/0.05;
      for (int i=0;i<fNmultisims;i++) {
        double test = bin < 200 ? fRW[ntype][bin]*fWeightArray[i] : 1;
        // cap weights at +/- 100%
        if(test >= 2.){test = 2.;}
        if(test <= 0.){test = 1.e-10;}
        // Guards against inifinite weights
        if(std::isfinite(test)){ weight[inu][i] = test;}
        else{weight[inu][i] = 1;}

      }
    }
    return weight;
  }
  REGISTER_WEIGHTCALC(NuMIFluggHistWeightCalc)
}
