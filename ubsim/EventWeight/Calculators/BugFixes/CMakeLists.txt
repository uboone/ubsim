include_directories ( $ENV{LIBXML2_INC} )
include_directories ( $ENV{GENIE}/src )

# Needed to build SplineWeightCalc. Note that the GENIE ones need to be changed
# if building against v2
# TODO: add a flag to determine the right GENIE version, as in current versions of nutools
cet_find_library( XML2 NAMES xml2 PATHS $ENV{LIBXML2_FQ_DIR}/lib NO_DEFAULT_PATH )
cet_find_library( GFWGHEP NAMES GFwGHEP PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GFWINT NAMES GFwInt PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GFWNTP NAMES GFwNtp PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GFWNUM NAMES GFwNum PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GFWUTL NAMES GFwUtl PATHS ENV GENIE_LIB NO_DEFAULT_PATH )

art_make(
  LIBRARY_NAME ubsim_EventWeight_Calculators_BugFixes
  LIB_LIBRARIES
    ${XML2}
    ${GHEP}
    ${GINTERACTION}
    ${GNTUPLE}
    ${GUTILS}
    larsim_EventWeight_Base
    lardata_Utilities
    nutools_NuReweight_art   
    nutools_NuReweight       
    larsim_Simulation
    nutools_RandomUtils_NuRandomService_service
    nutools_EventGeneratorBase_GENIE
    nusimdata_SimulationBase
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
    art_Persistency_Common canvas canvas
    art_Persistency_Provenance canvas canvas
    art_Utilities canvas canvas
    cetlib cetlib_except
    ${FHICLCPP}
    ${CLHEP}
    ${ROOT_BASIC_LIB_LIST}
    ${Boost_SYSTEM_LIBRARY}
    ${G4_LIB_LIST}
    ${MF_MESSAGELOGGER}
    ${GFWGHEP}
    ${GFWINT}
    ${GFWNTP}
    ${GFWNUM}
    ${GFWUTL}
) 

install_headers()
install_fhicl()
install_source()
