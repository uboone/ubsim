cet_build_plugin(
  UBPhotonLibraryPropagation art::EDProducer
  LIBRARIES
  PRIVATE
  larsim::IonizationScintillation
  larsim::PhotonPropagation_PhotonVisibilityService_service
  larsim::PhotonPropagation
  larevt::SpaceChargeService
  lardata::LArPropertiesService
  larcore::Geometry_Geometry_service
  nurandom::RandomUtils_NuRandomService_service
)

install_fhicl()
install_source()
