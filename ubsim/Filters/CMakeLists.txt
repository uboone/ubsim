cet_build_plugin(
  CorsikaCosmicGen art::EDProducer
  LIBRARIES
  PRIVATE
  larcore::Geometry_Geometry_service
  nusimdata::SimulationBase
  nurandom::RandomUtils_NuRandomService_service
  nutools::EventGeneratorBase_CRY
  CRY::CRY
  ROOT::EG
)

cet_build_plugin(
  EtaFilter art::EDFilter
  LIBRARIES
  PRIVATE
  nusimdata::SimulationBase
)

cet_build_plugin(
  HighEnergyNuMINuMuCCFilter art::EDFilter
  LIBRARIES
  PRIVATE
  lardataobj::MCBase
  nusimdata::SimulationBase
)

cet_build_plugin(
  NeutrinoTruthFilter art::EDFilter
  LIBRARIES
  PRIVATE
  nusimdata::SimulationBase
)

cet_build_plugin(
  NuMIKDARFilter art::EDFilter
  LIBRARIES
  PRIVATE
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  nusimdata::SimulationBase
)

cet_build_plugin(
  NueBackgroundFilter art::EDFilter
  LIBRARIES
  PRIVATE
  lardataobj::MCBase
  nusimdata::SimulationBase
)

cet_build_plugin(
  StopMuFilter art::EDFilter
  LIBRARIES
  PRIVATE
  larevt::SpaceChargeService
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  nusimdata::SimulationBase
)

install_headers()
install_fhicl()
install_source()

add_subdirectory(jobs)
