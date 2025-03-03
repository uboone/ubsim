cet_build_plugin(G4RWManagerService art::service
  geant4reweight::ReweightBaseLib
  fhiclcpp::fhiclcpp
)

cet_make_library(
  SOURCE
  FluxHistWeightCalc.cxx
  FluxUnisimWeightCalc.cxx
  Geant4WeightCalc.cxx
  HistogramWeightWeightCalc.cxx
  LEESignalElectronWeightCalc.cxx
  ReinteractionWeightCalc.cxx
  SCCWeightCalc.cxx
  UBGenieWeightCalc.cxx
  #UBOLDPPFXCVWeightCalc.cxx
  #UBOLDPPFXWeightCalc.cxx
  UBPPFXCVWeightCalc.cxx
  UBPPFXWeightCalc.cxx
  LIBRARIES
  PUBLIC
  ubsim_EventWeight_Calculators_G4RWManagerService_service
  larsim::EventWeight_Base
  larcore::Geometry_Geometry_service
  nusimdata::SimulationBase
  nugen::EventGeneratorBase_GENIE
  Geant4::G4event
  ppfx::ppfx
  dk2nu::Tree
  art_root_io::TFileService_service
  fhiclcpp::fhiclcpp
  geant4reweight::ReweightBaseLib
  GENIE::GRwClc
  GENIE::GFwInt
  GENIE::GRwFwk
  GENIE::GFwMsg
  GENIE::GFwUtl
  log4cpp::log4cpp
  ROOT::Hist
  ROOT::Tree
)

install_headers()
install_fhicl()
install_source()

add_subdirectory(BNBPrimaryHadron)
add_subdirectory(Splines)
