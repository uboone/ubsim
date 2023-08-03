cet_build_plugin(
  GeneratorPositionFilter art::EDFilter
  LIBRARIES
  PRIVATE
  nusimdata::SimulationBase
)

cet_build_plugin(
  MultiPartVertex art::EDProducer
  LIBRARIES
  PRIVATE
  larcore::Geometry_Geometry_service
  nurandom::RandomUtils_NuRandomService_service
  nusimdata::SimulationBase
  ROOT::EG
)

add_subdirectory(GENIE)
add_subdirectory(HeavySterileNeutrino)
add_subdirectory(HiggsPortalScalarGen)
add_subdirectory(HeavyNeutralLeptonGen)

install_fhicl()
install_source()
