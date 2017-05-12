# canvas_utlities
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{NUTOOLS_INC} )

art_make( EXCLUDE
    CRTAuxFunctions.cc
    CRTCalibration_module.cc
    CRTChannelMapAlg.cc
    CRTCoincidence_module.cc
    CRTHitProducer_module.cc	
    CRTDaqSim_module.cc
    CRTDetSim_module.cc
    CRTGeometryHelper_service.cc
    CRTMerger_module.cc
    CRTTriggerFilter_module.cc
    CRTGeoObjectSorter.cc
    CRTChannelMapAlg.cc
    CRTAuxFunctions.cc
    CRTDaqSim_module.cc
    CRTMerger_module.cc
  BASENAME_ONLY
  LIBRARY_NAME 
    CRTSimData
  LIB_LIBRARIES
    ${ART_FRAMEWORK_CORE} 

  MODULE_LIBRARIES
    ${CETLIB}
    ${ROOT_GEOM}
    ${ROOT_XMLIO}
    ${ROOT_GDML}
    ${ROOT_BASIC_LIB_LIST}
    art_Persistency_Common canvas_Persistency_Common
    art_Persistency_Provenance canvas_Persistency_Provenance
    art_Utilities canvas_Utilities
)

art_make_library( LIBRARY_NAME CRT
  SOURCE
    CRTGeoObjectSorter.cc
    CRTChannelMapAlg.cc
    CRTAuxFunctions.cc
  LIBRARIES
    CRTSimData
    larcore_Geometry
    ${ROOT_BASIC_LIB_LIST}
    ${CETLIB}
    ${FHICLCPP}
)

simple_plugin( CRTDaqSim module
    larcore_Geometry
    CRTSimData
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
    CRTSimData
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


simple_plugin(CRTHitProducer "module"
            CRT
              bernfebdaq_core_overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin( CRTMerger module
  bernfebdaq_core_overlays
  ${ROOT_BASIC_LIB_LIST}
  ${ART_UTILITIES}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
  ${ARTDAQ_CORE}
  ${FHICLCPP}
  ${MF_MESSAGELOGGER}
  ${MF_UTILITIES}
  ${CETLIB}
  ${GALLERY}
  CRTSimData
)


install_headers()
install_fhicl()
install_source()
