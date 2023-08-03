cet_build_plugin(
  UBOpticalChConfig art::service
  LIBRARIES
  PUBLIC
  ubcore::Geometry_UBOpReadoutMap_service
  larcore::Geometry_Geometry_service
  lardata::LArPropertiesService
)

cet_make_library(
  SOURCE
  RandomServer.cxx
  SimpleChConfig.cxx
  UBADCBase.cxx
  UBLogicPulseADC.cxx
  UBOpticalADC.cxx
  UBOpticalException.cxx
  WFAlgoAnalyticalSPE.cxx
  WFAlgoArray.cxx
  WFAlgoBase.cxx
  WFAlgoDigitizedSPE.cxx
  WFAlgoPedestal.cxx
  WFAlgoSPEBase.cxx
  WFAlgoUtilities.cxx
  LIBRARIES
  PUBLIC
  ubsim::OpticalDetectorSim_UBOpticalChConfig_service
  lardataalg::DetectorInfo
  larcore::Geometry_Geometry_service
  lardata::DetectorClocksService
  ROOT::MathCore
)

cet_build_plugin(
  FakePhotonsAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardataobj::RawData
  ROOT::Tree
)

cet_build_plugin(
  FakePhotons art::EDProducer
  LIBRARIES
  PRIVATE
  larcore::Geometry_Geometry_service
  lardataobj::Simulation
  art_root_io::TFileService_service
)

cet_build_plugin(
  LYSimPhotonScaling art::EDProducer
  LIBRARIES
  PRIVATE
  larevt::CalibrationDBI_IOVData
  larcore::Geometry_Geometry_service
  lardataobj::Simulation
  nurandom::RandomUtils_NuRandomService_service
)

cet_build_plugin(
  OpticalDRAMReadout art::EDProducer
  LIBRARIES
  PRIVATE
  ubcore::Geometry_UBOpReadoutMap_service
  lardata::DetectorClocksService
  larcore::ServiceUtil
  lardataobj::headers
  art_root_io::TFileService_service
)

cet_build_plugin(
  OpticalFEM art::EDProducer
  LIBRARIES
  PRIVATE
  ubsim::OpticalDetectorSim_UBOpticalChConfig_service
  ubcore::Geometry_UBOpReadoutMap_service
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::headers
  art_root_io::TFileService_service
)

cet_build_plugin(
  UBFlasherMC art::EDProducer
  LIBRARIES
  PRIVATE
  ubsim::OpticalDetectorSim
  art_root_io::TFileService_service
)

cet_build_plugin(
  UBOpticalADCSim art::EDProducer
  LIBRARIES
  PRIVATE
  ubsim::OpticalDetectorSim
  art_root_io::TFileService_service
)

install_headers()
install_fhicl()
install_source()
