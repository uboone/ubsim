include_directories ( $ENV{GEANT4_FQ_DIR}/include )

art_make(
  LIBRARY_NAME ubsim_EventWeight_Calculators
  LIB_LIBRARIES
    larsim_EventWeight_Base
    lardata_Utilities
    nugen::NuReweight_art
    nugen::NuReweight
    larsim_Simulation
    nurandom::RandomUtils_NuRandomService_service
    nusimdata::SimulationBase
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    art::Framework_Services_Optional_RandomNumberGenerator_service
    art::Persistency_Common canvas::canvas canvas::canvas
    art::Persistency_Provenance canvas::canvas canvas::canvas
    art::Utilities canvas::canvas canvas::canvas
    cetlib::cetlib cetlib_except
    fhiclcpp::fhiclcpp
    CLHEP::CLHEP

    ${G4_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()

add_subdirectory(BNBPrimaryHadron)

