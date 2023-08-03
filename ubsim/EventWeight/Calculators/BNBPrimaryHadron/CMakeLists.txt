cet_make_library(
  SOURCE
  PrimaryHadronFeynmanScalingWeightCalc.cxx
  PrimaryHadronNormalizationWeightCalc.cxx
  PrimaryHadronSWCentralSplineVariationWeightCalc.cxx
  PrimaryHadronSanfordWangWeightCalc.cxx
  LIBRARIES
  PUBLIC
  larsim::EventWeight_Base
  art::Framework_Principal
  nusimdata::SimulationBase
  ROOT::Hist
)

install_headers()
install_fhicl()
install_source()


