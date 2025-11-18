cet_build_plugin(
  TwoBodyDecayGen art::EDProducer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  art_root_io::TFileService_service
  nurandom::RandomUtils_NuRandomService_service
  nusimdata::SimulationBase
  nugen::EventGeneratorBase
  ROOT::EG
)

install_headers()
install_fhicl()
install_source()
