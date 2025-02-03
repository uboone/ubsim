cet_make_library(
  SOURCE
  FluxHistWeightCalc.cxx
  FluxUnisimWeightCalc.cxx
  HistogramWeightWeightCalc.cxx
  LEESignalElectronWeightCalc.cxx
  ReinteractionWeightCalc.cxx
  SCCWeightCalc.cxx
  UBPPFXCVWeightCalc.cxx
  UBPPFXWeightCalc.cxx
  Geant4WeightCalc.cxx
  LIBRARIES
  PUBLIC
  larsim::EventWeight_Base
  nusimdata::SimulationBase
  nugen::EventGeneratorBase_GENIE
  Geant4::G4event
  ppfx::ppfx
  dk2nu::Tree
  art_root_io::TFileService_service
  ROOT::Hist
  ROOT::Tree
  ROOT::MathCore
  geant4reweight::ReweightBaseLib
  larcorealg::Geometry
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
