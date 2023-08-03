cet_build_plugin(
  RawDigitSimulator art::EDProducer
  LIBRARIES
  PRIVATE
  ubevt::Utilities_SignalShapingServiceMicroBooNE_service
  larevt::ElectronicsCalibService
  lardata::DetectorClocksService
  lardata::Utilities_LArFFT_service
  larcore::Geometry_Geometry_service
  lardataobj::Simulation
  art_root_io::TFileService_service
  nurandom::RandomUtils_NuRandomService_service
  ROOT::MathCore
)

cet_build_plugin(
  SimWireMicroBooNEAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardataobj::RawData
)

cet_build_plugin(
  SimWireMicroBooNE art::EDProducer
  LIBRARIES
  PRIVATE
  ubevt::Utilities_SignalShapingServiceMicroBooNE_service
  larevt::DetPedestalService
  lardata::Utilities_LArFFT_service
  larcore::Geometry_Geometry_service
  lardataobj::Simulation
  lardataobj::RawData
  art_root_io::TFileService_service
  nurandom::RandomUtils_NuRandomService_service
  ROOT::MathCore
)

# install_headers()
install_fhicl()
install_source()

add_subdirectory(MicroBooNEWireCell)
