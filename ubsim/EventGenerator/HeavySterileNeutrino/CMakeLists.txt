cet_build_plugin(
  HSNGen art::EDProducer
  LIBRARIES
  PRIVATE
  ubsim::EventGenerator_HeavySterileNeutrino_Helpers
  larcore::Geometry_Geometry_service
  art_root_io::TFileService_service
  nurandom::RandomUtils_NuRandomService_service
  nusimdata::SimulationBase
  nutools::EventGeneratorBase_CRY
  GENIE::GFwParDat
  ROOT::Tree
  ROOT::MathCore
  SQLite::SQLite3
)

add_subdirectory(Helpers)
add_subdirectory(DataObjects)

install_headers()
install_fhicl()
install_source()

