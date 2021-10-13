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
    nugen_EventGeneratorBase_GENIE
    nugen_NuReweight_art   
    nugen_NuReweight       
    larsim_Simulation
    larcorealg_Geometry
    larcore_Geometry_Geometry_service
    nurandom_RandomUtils_NuRandomService_service
    nusimdata_SimulationBase
    #ReweightBaseLib
    #PropBaseLib
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
    ${GENIE_LIB_LIST}
    ${FHICLCPP}
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
