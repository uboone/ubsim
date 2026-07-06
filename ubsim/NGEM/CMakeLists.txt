cet_add_compiler_flags( CXX -DHIDE_GENIE_LOG_XXX )

cet_build_plugin(
  GENIEGenNGEM art::EDProducer
  LIBRARIES
  PRIVATE
  lardata::AssociationUtil
  lardataalg::MCDumpers
  larcore::Geometry_Geometry_service
  lardataobj::Simulation
  nurandom::RandomUtils_NuRandomService_service
  nugen::EventGeneratorBase
  nugen::EventGeneratorBase_GENIE
  art_root_io::TFileService_service
  dk2nu::Genie
  dk2nu::Tree
  GENIE::GFwEG
  ROOT::Tree
  ROOT::EG
)

install_headers()
install_fhicl()
install_source()

