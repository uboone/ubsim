art_make(
    EXCLUDE UBPhotonLibraryPropagation_module.cc
    MODULE_LIBRARIES
        larsim_PhotonPropagation
        larsim_PhotonPropagation_PhotonVisibilityService_service
        larsim_IonizationScintillation
        lardataobj_Simulation
        larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_PRINCIPAL}
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
        nurandom_RandomUtils_NuRandomService_service
        art_Persistency_Common
        art_Persistency_Provenance
        art_Utilities
        canvas
        ${MF_MESSAGELOGGER}
        ${FHICLCPP}
        cetlib cetlib_except
        ${CLHEP}
        ${ROOT_BASIC_LIB_LIST}
)

install_fhicl()
install_source()
