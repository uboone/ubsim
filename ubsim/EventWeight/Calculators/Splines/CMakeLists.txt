include_directories ( $ENV{LIBXML2_INC} )
include_directories ( $ENV{GENIE}/src )

# Needed to build SplineWeightCalc. Note that the GENIE ones need to be changed
# if building against v2
# TODO: add a flag to determine the right GENIE version, as in current versions of nutools
#cet_find_library( XML2 NAMES xml2 PATHS $ENV{LIBXML2_FQ_DIR}/lib NO_DEFAULT_PATH )
#cet_find_library( GFWGHEP NAMES GFwGHEP PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
#cet_find_library( GFWINT NAMES GFwInt PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
#cet_find_library( GFWNTP NAMES GFwNtp PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
#cet_find_library( GFWNUM NAMES GFwNum PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
#cet_find_library( GFWUTL NAMES GFwUtl PATHS ENV GENIE_LIB NO_DEFAULT_PATH )

art_make(
  LIBRARY_NAME ubsim_EventWeight_Calculators_Splines
  LIB_LIBRARIES
    xml2
    larsim_EventWeight_Base
    lardata_Utilities
    nugen::NuReweight_art   
    nugen::NuReweight       
    larsim_Simulation
    nurandom::RandomUtils_NuRandomService_service
    nugen::EventGeneratorBase_GENIE
    nusimdata::SimulationBase
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    art::Framework_Services_Optional_RandomNumberGenerator_service
    art::Persistency_Common canvas canvas
    art::Persistency_Provenance canvas canvas
    art::Utilities canvas canvas
    cetlib::cetlib
    cetlib_except::cetlib_except
    fhiclcpp::fhiclcpp
    CLHEP::CLHEP
    
    messagefacility::MF_MessageLogger
    GENIE::GFwGHEP
    GENIE::GFwInt
    GENIE::GFwNtp
    GENIE::GFwNum
    GENIE::GFwUtl
) 

install_headers()
install_fhicl()
install_source()
