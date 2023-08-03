cet_make_library(
  SOURCE
  FluxReader.cxx
  FluxReaderBNB.cxx
  FluxReaderNuMI.cxx
  GenKinematics.cxx
  LIBRARIES
  PUBLIC
  larcore::Geometry_Geometry_service
  ifdh_art::IFDH_service
  art::Framework_Services_Registry
  dk2nu::Tree
  GENIE::GFwParDat
  fhiclcpp::fhiclcpp
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  HiggsPortalScalarGenFromBNBFlux art::EDProducer
  LIBRARIES
  PRIVATE
  ubsim::EventGenerator_HiggsPortalScalarGen
  lardata::Utilities
  art_root_io::TFileService_service
  nurandom::RandomUtils_NuRandomService_service
)

cet_build_plugin(
  HiggsPortalScalarGenFromNuMIFlux art::EDProducer
  LIBRARIES
  PRIVATE
  ubsim::EventGenerator_HiggsPortalScalarGen
  lardata::Utilities
  art_root_io::TFileService_service
  nurandom::RandomUtils_NuRandomService_service
)

add_subdirectory(job)

install_headers()
install_fhicl()
install_source()

