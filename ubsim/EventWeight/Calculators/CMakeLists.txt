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
)

install_headers()
install_fhicl()
install_source()

add_subdirectory(BNBPrimaryHadron)
add_subdirectory(Splines)
