/**
 * Histogram-based reweighting calculator.
 *
 * Reweights MC based on an externally provided histogram, e.g. with a
 * cross section ratio between two models. Currently limited to a TH2F
 * binned in q0 vs. q3, and somewhat ambitiously named.
 *
 * The calculator is configured with a few FHICL parameters:
 *
 *     genie_module_label   string  Name of MCTruth producer
 *     mode                 string  "multisim", or any other value for unisims
 *     number_of_multisims  int     Number of multisim universes
 *
 *     rw_hist_file         string  Name of the histogram ROOT file
 *     rw_hist_object       string  Name of the histogram object in the file
 *     sigma                double  Scale of shifts (alternative model is 1 sigma)
 *     norm_scale           double  Additional normalization scale factor
 *
 *     event_filter         string  Event type ("ccqe" or "ccmec")
 *
 * The file is located based on the FW_SEARCH_PATH environment variable.
 *
 * The "event_filter" parameter is used to apply the reweighting only to
 * specific classes of events. Currently this includes CCQE and CCMEC.
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>, 2018/04/30
 */

#include <iostream>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Persistency/Provenance/ModuleContext.h"
#include "nurandom/RandomUtils/NuRandomService.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "TFile.h"
#include "TH2F.h"
#include "larsim/EventWeight/Base/WeightCalcCreator.h"
#include "larsim/EventWeight/Base/WeightCalc.h"

#include "CLHEP/Random/RandGaussQ.h"

namespace evwgh {

class HistogramWeightWeightCalc : public WeightCalc {
public:
  /** Constructor. */
  HistogramWeightWeightCalc();

  /** Configure the calculator based on FHICL parameters. */
  void Configure(fhicl::ParameterSet const& p,
                 CLHEP::HepRandomEngine& engine);

  /** Compute the list of weights. */
  std::vector<std::vector<double> > GetWeight(art::Event& e);

private:
  CLHEP::RandGaussQ* fGaussRandom;  //!< Random number generator
  std::vector<double> fWeightArray;  //!< Random numbers for weights
  double fNormScale;  //!< Normalization scaling
  int fNmultisims;  //!< Number of multisim universes
  std::string fGenieModuleLabel;  //!< Module name for MCTruth
  std::string fMode;  //!< Multisim vs. unisim mode
  std::string fEventFilter;  //!< Event filter string (see notes in header)
  bool fOneSided;  //!< One-sided, use upper half of a standard normal
  TH2F* fRWHist;  //!< Reweighting histogram

  DECLARE_WEIGHTCALC(HistogramWeightWeightCalc)
};


HistogramWeightWeightCalc::HistogramWeightWeightCalc() : fRWHist(NULL) {}


void HistogramWeightWeightCalc::Configure(fhicl::ParameterSet const& p,
                                          CLHEP::HepRandomEngine& engine)
{
  fhicl::ParameterSet const& pset = p.get<fhicl::ParameterSet>(GetName());

  // Global config
  fGenieModuleLabel = p.get<std::string>("genie_module_label");
  fGaussRandom = new CLHEP::RandGaussQ(engine);

  // Load weight histogram from a file
  double sigma = pset.get<double>("sigma");
  std::string histFileName = pset.get<std::string>("rw_hist_file");
  std::string histObjectName = pset.get<std::string>("rw_hist_object");

  cet::search_path sp("FW_SEARCH_PATH");
  std::string histFilePath = sp.find_file(histFileName);
  TFile histFile(Form("%s", histFilePath.c_str()));
  assert(histFile.IsOpen());

  TH2F* h = dynamic_cast<TH2F*>(histFile.Get(histObjectName.c_str()));
  assert(h);
  fRWHist = dynamic_cast<TH2F*>(h->Clone("rwhist"));
  fRWHist->SetDirectory(NULL);

  histFile.Close();

  // Calculator config
  fNmultisims = pset.get<int>("number_of_multisims");
  fMode = pset.get<std::string>("mode");
  fEventFilter = pset.get<std::string>("event_filter");
  fOneSided = pset.get<bool>("one_sided", true);
  fNormScale = pset.get<double>("norm_scale");
  assert(fEventFilter == "ccqe" || fEventFilter == "ccmec");

  // Initialize the list of random numbers for each universe
  fWeightArray.resize(fNmultisims);
  for (int i=0; i<fNmultisims; i++) {
    if (fMode == "multisim") {
      fWeightArray[i] = fGaussRandom->fire(0, sigma);

      // One sided uncertainty: upper half of standard normal
      if (fOneSided) {
        fWeightArray[i] = std::abs(fWeightArray[i]);
      }
    }
    else {
      assert(fNmultisims == 1);
      fWeightArray[i] = sigma;
    }
  }
}


std::vector<std::vector<double> >
HistogramWeightWeightCalc::GetWeight(art::Event& e) {
  std::vector<std::vector<double> > weight;

  // MC truth information
  auto mctruthListHandle = e.getHandle<std::vector<simb::MCTruth> >(fGenieModuleLabel);
  if (!mctruthListHandle) {
    return weight;
  }
  std::vector<simb::MCTruth> const& mclist = *mctruthListHandle;

  weight.resize(mclist.size());

  // Loop through MCTruth neutrinos
  for (unsigned int inu=0; inu<mclist.size(); inu++) {
    weight[inu].resize(fNmultisims, 1.0);

    // Truth-level event filtering and kinematics
    simb::MCNeutrino nu = mclist[inu].GetNeutrino();
    simb::MCParticle lep = nu.Lepton();

    int mode = nu.Mode();
    int ccnc = nu.CCNC();
    double q0 = nu.Nu().E() - lep.E();
    double q3 = (nu.Nu().Momentum().Vect() - lep.Momentum().Vect()).Mag();

    // Loop through multisim universes
    for (int i=0; i<fNmultisims; i++) {
      if ((fEventFilter == "ccqe"  && ccnc == simb::kCC && mode == simb::kQE ) ||
          (fEventFilter == "ccmec" && ccnc == simb::kCC && mode == simb::kMEC)) {
        int xbin = fRWHist->GetXaxis()->FindBin(q3);
        int ybin = fRWHist->GetYaxis()->FindBin(q0);
        double w = fRWHist->GetBinContent(xbin, ybin);
        weight[inu][i] = 1.0 - (1.0 - fNormScale * w) * fWeightArray[i];
        weight[inu][i] = std::max(0.0, weight[inu][i]);

        //std::cout << "Histogram weight: "
        //          << "type = " << fEventFilter << ", "
        //          << "q0 = " << q0 << ", "
        //          << "q3 = " << q3 << ", "
        //          << "w = " << weight[inu][i] << std::endl;
      }
      else {
        weight[inu][i] = 1.0; 
      }
    }
  }

  return weight;
}

REGISTER_WEIGHTCALC(HistogramWeightWeightCalc)

}  // namespace evwgh
