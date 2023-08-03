cet_make_library(
  SOURCE
  FluxReader.cxx
  FluxReaderBNB.cxx
  FluxReaderNuMI.cxx
  GenKinematics.cxx
  LIBRARIES
  PUBLIC
  larcore::Geometry_Geometry_service
  art::Framework_Services_Registry
  ifdh_art::IFDH_service
  fhiclcpp::fhiclcpp
  dk2nu::Tree
  GENIE::GFwParDat
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  HeavyNeutralLeptonGenFromBNBFlux art::EDProducer
  LIBRARIES
  PRIVATE
  ubsim::EventGenerator_HeavyNeutralLeptonGen
  lardata::Utilities
  art_root_io::TFileService_service
  nurandom::RandomUtils_NuRandomService_service
)

cet_build_plugin(
  HeavyNeutralLeptonGenFromNuMIFlux art::EDProducer
  LIBRARIES
  PRIVATE
  ubsim::EventGenerator_HeavyNeutralLeptonGen
  lardata::Utilities
  art_root_io::TFileService_service
  nurandom::RandomUtils_NuRandomService_service
)

add_subdirectory(job)

install_headers()
install_fhicl()
install_source()
