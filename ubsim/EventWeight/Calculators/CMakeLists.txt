include_directories ( $ENV{GEANT4_FQ_DIR}/include )

# For UBGenieWeightCalc
include_directories ( $ENV{GENIE_INC} )
include_directories ( $ENV{GENIE_INC}/GENIE )
include_directories ( $ENV{LOG4CPP_INC} )
include_directories ( $ENV{PPFX_DIR}/include )
include_directories ( $ENV{DK2NU_INC} )

art_make(
  LIBRARY_NAME ubsim_EventWeight_Calculators
  EXCLUDE UBGenieWeightCalc.cxx Geant4WeightCalc.cxx
  LIB_LIBRARIES
    larsim::EventWeight_Base
    lardata::Utilities
    nugen::NuReweight_art
    nugen::NuReweight
    nugen::EventGeneratorBase_GENIE
    larsim::Simulation
    nurandom::RandomUtils_NuRandomService_service
    nusimdata::SimulationBase
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    art::Framework_Services_Optional_RandomNumberGenerator_service
    art::Persistency_Common
    art::Persistency_Provenance
    art::Utilities
    canvas::canvas
    cetlib::cetlib
    cetlib_except::cetlib_except
    fhiclcpp::fhiclcpp
    ppfx::ppfx
    dk2nu::Tree
    dk2nu::Genie
    CLHEP::CLHEP

    log4cpp::log4cpp
    messagefacility::MF_MessageLogger
)

install_headers()
install_fhicl()
install_source()

add_subdirectory(BNBPrimaryHadron)
add_subdirectory(Splines)
