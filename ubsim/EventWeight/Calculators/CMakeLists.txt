include_directories ( $ENV{GEANT4_FQ_DIR}/include )

# For UBGenieWeightCalc
include_directories ( $ENV{GENIE_INC} )
include_directories ( $ENV{GENIE_INC}/GENIE )
include_directories ( $ENV{LOG4CPP_INC} )

# Determine the major version number of GENIE based on the GENIE_VERSION
# variable set by the genie ups product. Define the GENIE_PRE_R3 preprocessor
# macro (needed to set the right GENIE header locations, etc.) if needed.
#
# The top-level nutools CMakeLists.txt does this as well, but I don't think
# that we can access the result of that version check in the ubsim build.
# GENIE v3 library list
set (GENIE_LIB_LIST GFwMsg
                    GFwReg
                    GFwAlg
                    GFwInt
                    GFwGHEP
                    GFwNum
                    gsl       # FWNUM relies on GSL
                    GFwUtl
                    GFwParDat
                    GFwEG
                    GFwNtp
                    GPhXSIg
                    GPhPDF
                    GPhNuclSt
                    GPhCmn
                    GPhDcy
                    GPhHadTransp
                    GPhHadnz
                    GPhDeEx
                    GPhAMNGXS
                    GPhAMNGEG
                    GPhChmXS
                    GPhCohXS
                    GPhCohEG
                    GPhDISXS
                    GPhDISEG
                    GPhDfrcXS
                    GPhDFRCEG
                    GPhGlwResXS
                    GPhGlwResEG
                    GPhIBDXS
                    GPhIBDEG
                    GPhMNucXS
                    GPhMNucEG
                    GPhMEL
                    GPhNuElXS
                    GPhNuElEG
                    GPhQELXS
                    GPhQELEG
                    GPhResXS
                    GPhResEG
                    GPhStrXS
                    GPhStrEG
                    GPhNDcy
                    GTlGeo
                    GTlFlx
                    GRwFwk
                    GRwIO
                    GRwClc)

# Build against extra reweighting functionality in GENIE v3.0.4 MicroBooNE patch 02
add_definitions(-DGENIE_UB_PATCH)

art_make(
  LIBRARY_NAME ubsim_EventWeight_Calculators
  EXCLUDE UBGenieWeightCalc.cxx Geant4WeightCalc.cxx
  LIB_LIBRARIES
    larsim_EventWeight_Base
    lardata_Utilities
    nugen::NuReweight_art
    nugen::NuReweight
    larsim_Simulation
    nurandom::RandomUtils_NuRandomService_service
    nusimdata::SimulationBase
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    art::Framework_Services_Optional_RandomNumberGenerator_service
    art::Persistency_Common canvas::canvas canvas::canvas
    art::Persistency_Provenance canvas::canvas canvas::canvas
    art::Utilities canvas::canvas canvas::canvas
    cetlib::cetlib cetlib_except
    fhiclcpp::fhiclcpp
    CLHEP::CLHEP

    ${G4_LIB_LIST}
    log4cpp
    ${MF_MESSAGELOGGER}
)

install_headers()
install_fhicl()
install_source()

add_subdirectory(BNBPrimaryHadron)
add_subdirectory(Splines)
