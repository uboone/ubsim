#include "larsim/EventWeight/Base/WeightCalcCreator.h"
#include "larsim/EventWeight/Base/WeightCalc.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Random/RandGaussQ.h"

#include "MakeReweight.h"
#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"
#include "nutools/EventGeneratorBase/GENIE/MCTruthAndFriendsItr.h"
#include "cetlib/search_path.h"

#include "TSystem.h"
#include "TFile.h"
#include "TH2F.h"

namespace evwgh {

template <typename T>
class HistReader
{
public:
  HistReader(TFile* fin, std::string key, double limit = -1.)
    : fFile(fin), fLimit(limit)
  {
      fHist = (T*)fFile->Get(key.c_str());
  }
  template <typename... dtypes>
  double operator()(dtypes... values) const {
      int binidx = fHist->FindBin(static_cast<double>(values)...);
      if(fLimit < 0.)
        return (double)(fHist->GetBinContent(binidx));
      else
        return (fHist->GetBinContent(binidx) < fLimit) ? (double)fHist->GetBinContent(binidx) : fLimit;
  }
private:
  T* fHist;
  TFile* fFile;
  double fLimit;
};

class UBOLDPPFXWeightCalc : public WeightCalc
  {
  public:
    UBOLDPPFXWeightCalc();
    void Configure(fhicl::ParameterSet const& p,
                   CLHEP::HepRandomEngine& engine) override;

    std::vector<std::vector<double> > GetWeight(art::Event & e) override;
  private:
    std::string fGenieModuleLabel;

    std::vector<std::string>  fInputLabels;
    std::string fPPFXMode;
    std::string fMode;
    std::string fHorn;
    std::string fTarget;
    int fSeed;
    int fVerbose;
    NeutrinoFluxReweight::MakeReweight* fPPFXrw;

    // reweighting
    std::string fOldFluxReweightFile;
    // order is numu, nue, numubar, nuebar
    std::vector<HistReader<TH2D>> fHists;
    double GetOldFluxReweight(const simb::MCFlux* mcflux);

    DECLARE_WEIGHTCALC(UBOLDPPFXWeightCalc)
  };

UBOLDPPFXWeightCalc::UBOLDPPFXWeightCalc()
{
}

double UBOLDPPFXWeightCalc::GetOldFluxReweight(const simb::MCFlux* mcflux)
{

  // get our variables of interest
  double Enu    = mcflux->fnenergyn;
  double L      = mcflux->fdk2gen + mcflux->fgen2vtx;
  int    nupdg  = mcflux->fntype;
  double nuvx   = mcflux->fvx;
  double nuvy   = mcflux->fvy;
  double nuvz   = mcflux->fvz;

  // no weight for bad volume regions
  if((std::abs(nuvx) > 58.42 || ((nuvy > 35.56 && nuvz > 350) || (nuvy > 65 && nuvz < 350))) && nuvz < 4569.9)
    return 0.;
  // get histogram based on order
  int idx = -1;
  if(nupdg == 14)
    idx = 0;
  if(nupdg == 12)
    idx = 1;
  if(nupdg == -14)
    idx = 2;
  if(nupdg == -12)
    idx = 3;

  if(idx < 0) return 1.;
  // return the weight
  return (fHists.at(idx))(Enu, L);
}

void UBOLDPPFXWeightCalc::Configure(fhicl::ParameterSet const& p,
                                      CLHEP::HepRandomEngine& engine)
{
  //get configuration for this function
  fhicl::ParameterSet const &pset=p.get<fhicl::ParameterSet> (GetName());
  fGenieModuleLabel = p.get<std::string> ("genie_module_label");

  //ppfx setup
  fInputLabels            = pset.get<std::vector<std::string>>("input_labels");
  fPPFXMode               = pset.get<std::string>("ppfx_mode");
  fOldFluxReweightFile    = pset.get<std::string>("oldflux_reweight");
  fVerbose                = pset.get<int>("verbose");
  fMode                   = pset.get<std::string>("mode");
  fHorn	                  = pset.get<std::string>("horn_curr");
  fTarget 	              = pset.get<std::string>("target_config");
  fSeed                   = pset.get<int>("random_seed");

  gSystem->Setenv("MODE", fPPFXMode.c_str());

  mf::LogWarning("evwgh::UBOLDPPFXWeightCalc.cxx")                   <<
    "\nUSE ONLY WHEN WORKING FROM OLD FLUX FILES AND"                <<
    "\nRUNNING WITH MODERN PPFX (TUNED TO G4.10.4)"                  <<
    "\nTHIS REWEIGHTS THE OLD FLUX PREDICTION TO THE NEW G4 VERSION" <<
    "\nAND THEN CALCULATES THE MODERN PPFX WEIGHT";


  fPPFXrw = NeutrinoFluxReweight::MakeReweight::getInstance();
  std::cout<<"PPFX instance "<<fPPFXrw<<std::endl;

  std::string inputOptions;                                 // Full path.
  std::string mode_file = "inputs_" + fPPFXMode + ".xml";   // Just file name.
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(mode_file, inputOptions);

  // file to add additional weights into ppfx weights
  TFile* reweightFile = new TFile(fOldFluxReweightFile.c_str(), "read");
  assert(reweightFile->IsOpen());
  // store the relevant reweight histograms split by flavor
  std::string horn_key = (fHorn == "-200i") ? "rhc" : "fhc";
  fHists.emplace_back(reweightFile, "ratio_cv_"+horn_key+"_enu_length_numu");
  fHists.emplace_back(reweightFile, "ratio_cv_"+horn_key+"_enu_length_nue");
  fHists.emplace_back(reweightFile, "ratio_cv_"+horn_key+"_enu_length_numubar");
  fHists.emplace_back(reweightFile, "ratio_cv_"+horn_key+"_enu_length_nuebar");

  if (fSeed != -1) fPPFXrw->setBaseSeed(fSeed); // Set the random seed
  std::cout << "is PPFX setup : " << fPPFXrw->AlreadyInitialized() << std::endl;
  std::cout << "Setting PPFX, inputs: " << inputOptions << std::endl;
  std::cout << "Setting Horn Current Configuration to: " << fHorn << std::endl;
  std::cout << "Setting Target Configuration to: " << fTarget << std::endl;
  if(!(fPPFXrw->AlreadyInitialized())){
    fPPFXrw->SetOptions(inputOptions);
  }
  std::cout << "PPFX just set with mode: " << fPPFXMode << std::endl;
}

std::vector<std::vector<double> > UBOLDPPFXWeightCalc::GetWeight(art::Event & e)
{
  std::vector<std::vector<double> > weight;
  evgb::MCTruthAndFriendsItr mcitr(e,fInputLabels);
  //calculate weight(s) here

  int nmctruth=0, nmatched=0;
  bool flag = true;
  int  ievt = -1;
  std::vector<art::Handle<std::vector<bsim::Dk2Nu>>> dk2nus2;
  e.getManyByType(dk2nus2);
  for (size_t dk2=0; dk2 < dk2nus2.size(); ++dk2) {
    art::Handle<std::vector<bsim::Dk2Nu>> dk2nus = dk2nus2[dk2];
  }
  while ( ( flag = mcitr.Next() ) ) {
    std::string label                  = mcitr.GetLabel();
    const simb::MCTruth*     pmctruth  = mcitr.GetMCTruth();
    // const simb::GTruth*  pgtruth  = mcitr.GetGTruth();
    const simb::MCFlux*      pmcflux   = mcitr.GetMCFlux();
    const bsim::Dk2Nu*       pdk2nu    = mcitr.GetDk2Nu();
    //not-used//const bsim::NuChoice*    pnuchoice = mcitr.GetNuChoice();
    // art::Ptr<simb::MCTruth>  mctruthptr = mcitr.GetMCTruthPtr();

    ++ievt;
    ++nmctruth;
    if ( fVerbose > 0 ) {
      std::cout << "FluxWeightCalculator [" << std::setw(4) << ievt << "] "
        << " label \"" << label << "\" MCTruth@ " << pmctruth
        << " Dk2Nu@ " << pdk2nu << std::endl;
    }
    if ( ! pdk2nu ) continue;
    ++nmatched;

    double oldflux_rwgt = GetOldFluxReweight(pmcflux);

    //RWH//bsim::Dk2Nu*  tmp_dk2nu  = fTmpDK2NUConversorAlg.construct_toy_dk2nu( &mctruth,&mcflux);
    //RWH//bsim::DkMeta* tmp_dkmeta = fTmpDK2NUConversorAlg.construct_toy_dkmeta(&mctruth,&mcflux);
    //RWH// those appear to have been memory leaks
    //RWH// herein begins the replacment for the "construct_toy_dkmeta"

    static bsim::DkMeta dkmeta_obj;        //RWH// create this on stack (destroyed when out-of-scope)  ... or static
    dkmeta_obj.tgtcfg  = fTarget;
    dkmeta_obj.horncfg = fHorn;
    (dkmeta_obj.vintnames).push_back("Index_Tar_In_Ancestry");
    (dkmeta_obj.vintnames).push_back("Playlist");
    bsim::DkMeta* tmp_dkmeta = &dkmeta_obj;
    //RWH// herein ends this block

    // sigh ....
    //RWH// this is the signature in NeutrinoFluxReweight::MakeReweight :
    //      //! create an interaction chain from the new dk2nu(dkmeta) format
    //      void calculateWeights(bsim::Dk2Nu* nu, bsim::DkMeta* meta);
    //RWH// these _should_ be "const <class>*" because we don't need to change them
    //RWH// and the pointers we get out of the ART record are going to be const.
    bsim::Dk2Nu* tmp_dk2nu = const_cast<bsim::Dk2Nu*>(pdk2nu);  // remove const-ness
    try {
      fPPFXrw->calculateWeights(tmp_dk2nu,tmp_dkmeta);
    } catch (...) {
      std::cerr<<"Failed to calcualte weight"<<std::endl;
      continue;
    }
    //Get weights:
    if (fMode=="reweight") {
      double ppfx_cv_wgt = fPPFXrw->GetCVWeight();
      std::vector<double> wvec = {ppfx_cv_wgt*oldflux_rwgt};
      weight.push_back(wvec);
    } else {
      std::vector<double> vmipppion    = fPPFXrw->GetWeights("MIPPNumiPionYields");
      std::vector<double> vmippkaon    = fPPFXrw->GetWeights("MIPPNumiKaonYields");
      std::vector<double> vtgtatt      = fPPFXrw->GetWeights("TargetAttenuation");
      std::vector<double> vabsorp      = fPPFXrw->GetWeights("TotalAbsorption");
      std::vector<double> vttpcpion    = fPPFXrw->GetWeights("ThinTargetpCPion");
      std::vector<double> vttpckaon    = fPPFXrw->GetWeights("ThinTargetpCKaon");
      std::vector<double> vttpcnucleon = fPPFXrw->GetWeights("ThinTargetpCNucleon");
      std::vector<double> vttncpion    = fPPFXrw->GetWeights("ThinTargetnCPion");
      std::vector<double> vttnucleona  = fPPFXrw->GetWeights("ThinTargetnucleonA");
      std::vector<double> vttmesoninc  = fPPFXrw->GetWeights("ThinTargetMesonIncident");
      std::vector<double> vothers      = fPPFXrw->GetWeights("Other");

      std::vector<double> tmp_vhptot;
      for(unsigned int iuniv=0;iuniv<vtgtatt.size();iuniv++){
        tmp_vhptot.push_back(float(vmipppion[iuniv]*vmippkaon[iuniv]*vtgtatt[iuniv]*vabsorp[iuniv]*vttpcpion[iuniv]*
                                   vttpckaon[iuniv]*vttpcnucleon[iuniv]*vttncpion[iuniv]*vttnucleona[iuniv]*
                                   vttmesoninc[iuniv]*vothers[iuniv])*oldflux_rwgt);
      }
      weight.push_back(tmp_vhptot);
    }
  }

  return weight;
}
REGISTER_WEIGHTCALC(UBOLDPPFXWeightCalc)
}
