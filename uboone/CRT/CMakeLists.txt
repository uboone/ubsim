# canvas_utlities
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{NUTOOLS_INC} )

art_make( EXCLUDE
    CRTGeoObjectSorter.cc
    CRTChannelMapAlg.cc
    CRTAuxFunctions.cc
    CRTDaqSim_module.cc
  BASENAME_ONLY
  LIBRARY_NAME CRTData
  LIB_LIBRARIES 
      larcore_Geometry
      ${ROOT_BASIC_LIB_LIST}
      ${CETLIB}
      ${FHICLCPP}
  MODULE_LIBRARIES 
)   



art_make_library( LIBRARY_NAME CRT
  SOURCE
    CRTData.cc
    CRTGeoObjectSorter.cc
    CRTChannelMapAlg.cc
    CRTAuxFunctions.cc
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
    bernfebdaq_core_overlays
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

simple_plugin(CRTCalibration "module"
        CRT 
            bernfebdaq_core_overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(CRTCoincidence "module"
            CRT
              bernfebdaq_core_overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

install_headers()
install_fhicl()
install_source()
