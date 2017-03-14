# canvas_utlities

art_make_library( LIBRARY_NAME CRT
  SOURCE
    CRTData.cc
    CRTGeoObjectSorter.cc
    CRTChannelMapAlg.cc
  LIBRARIES
    larcore_Geometry
    ${ROOT_BASIC_LIB_LIST}
    ${CETLIB}
    ${FHICLCPP}
)

simple_plugin( CRTDaqSim module
    larcore_Geometry
    CRT
    art_Utilities 
    canvas_Utilities
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${MF_MESSAGELOGGER}
    ${ROOT_BASIC_LIB_LIST}
    ${CLHEP}
    ${CETLIB}
    ${bernfebdaq-core_Overlays}
    ${ROOT_BASIC_LIB_LIST}
    ${ART_UTILITIES}
    ${ARTDAQ_CORE}
    ${FHICLCPP}
    ${MF_MESSAGELOGGER}
    ${CETLIB}
    ${GALLERY}  
  BASENAME_ONLY
)

simple_plugin( CRTDetSim module
    CRT
    art_Utilities 
    lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
    nutools_RandomUtils_NuRandomService_service
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
    ${MF_MESSAGELOGGER}
    ${ROOT_BASIC_LIB_LIST}
    ${CLHEP}
    ${CETLIB}
  BASENAME_ONLY
)

simple_plugin( CRTGeometryHelper service
    larcore_Geometry
    CRT
    art_Utilities 
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${MF_MESSAGELOGGER}
    ${ROOT_BASIC_LIB_LIST}
    ${CETLIB}
  BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
