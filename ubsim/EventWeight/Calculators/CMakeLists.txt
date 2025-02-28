cet_build_plugin(G4RWManagerService art::service
  geant4reweight::ReweightBaseLib
  fhiclcpp::fhiclcpp
)

cet_make_library(
  SOURCE
  FluxHistWeightCalc.cxx
  FluxUnisimWeightCalc.cxx
  #Geant4WeightCalc.cxx
  HistogramWeightWeightCalc.cxx
  LEESignalElectronWeightCalc.cxx
  ReinteractionWeightCalc.cxx
  SCCWeightCalc.cxx
  UBGenieWeightCalc.cxx
  #UBOLDPPFXCVWeightCalc.cxx
  #UBOLDPPFXWeightCalc.cxx
  UBPPFXCVWeightCalc.cxx
  UBPPFXWeightCalc.cxx
  Geant4WeightCalc.cxx
  LIBRARIES
  PUBLIC
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
  ROOT::MathCore
  geant4reweight::ReweightBaseLib
  larcorealg::Geometry
  ubsim_EventWeight_Calculators_G4RWManagerService_service
)

# get_cmake_property(_variableNames VARIABLES)
# list (SORT _variableNames)
# foreach (_variableName ${_variableNames})
#     message(STATUS "${_variableName}=${${_variableName}}")
# endforeach()

message("LARCORE: ${larcore_INCLUDE_DIR}")

target_include_directories(
  ubsim_EventWeight_Calculators
  PUBLIC
  ${larcore_INCLUDE_DIR}
  ${larcorealg_INCLUDE_DIR}  
  $ENV{RANGE_INC}
  ${geant4reweight_INCLUDE_DIR}
)

install_headers()
install_fhicl()
install_source()

add_subdirectory(BNBPrimaryHadron)
add_subdirectory(Splines)
