include_directories ( $ENV{GEANT4_FQ_DIR}/include )

# For UBGenieWeightCalc
include_directories ( $ENV{GENIE_INC} )
include_directories ( $ENV{GENIE_INC}/GENIE )
include_directories ( $ENV{LOG4CPP_INC} )

# Build against extra reweighting functionality in GENIE v3.0.4 MicroBooNE patch 02
add_definitions(-DGENIE_UB_PATCH)

art_make(
  LIBRARY_NAME ubsim_EventWeight_Calculators
  LIB_LIBRARIES
    larsim_EventWeight_Base
    lardata_Utilities
    nutools_EventGeneratorBase_GENIE
    nutools_NuReweight_art
    nutools_NuReweight
    larsim_Simulation
    larcorealg_Geometry
    larcore_Geometry_Geometry_service
    nutools_RandomUtils_NuRandomService_service
    nusimdata_SimulationBase
    ReweightBaseLib
    PropBaseLib
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
    art_Persistency_Common canvas canvas
    art_Persistency_Provenance canvas canvas
    art_Utilities canvas canvas
    cetlib cetlib_except
    ${FHICLCPP}
    ${CLHEP}
    ${ROOT_BASIC_LIB_LIST}
    ${Boost_SYSTEM_LIBRARY}
    ${G4_LIB_LIST}
    ${GENIE_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()

add_subdirectory(BNBPrimaryHadron)
add_subdirectory(Splines)
