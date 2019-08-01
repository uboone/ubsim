art_make(
    MODULE_LIBRARIES
        lardataobj_Simulation
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_PRINCIPAL}
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
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
