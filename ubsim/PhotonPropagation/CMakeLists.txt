art_make(
    EXCLUDE UBPhotonLibraryPropagation_module.cc
    MODULE_LIBRARIES
        larsim_PhotonPropagation
        larsim_PhotonPropagation_PhotonVisibilityService_service
        larsim_IonizationScintillation
        lardataobj_Simulation
        larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
        art::Framework_Core
        art::Framework_Principal
        art::Framework_Services_Registry
        
        art::Framework_Services_Optional_RandomNumberGenerator_service
        nurandom::RandomUtils_NuRandomService_service
        art::Persistency_Common
        art::Persistency_Provenance
        art::Utilities
        canvas::canvas
        messagefacility::MF_MessageLogger
        fhiclcpp::fhiclcpp
        cetlib::cetlib cetlib_except
        CLHEP::CLHEP
        ${ROOT_BASIC_LIB_LIST}
)

install_fhicl()
install_source()
