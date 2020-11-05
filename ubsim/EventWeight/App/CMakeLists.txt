set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} )

include_directories ( $ENV{GENIE_INC} )
include_directories ( $ENV{GENIE_INC}/GENIE )
include_directories ( $ENV{LOG4CPP_INC} )

# Determine the major version number of GENIE based on the GENIE_VERSION
# variable set by the genie ups product. Define the GENIE_PRE_R3 preprocessor
# macro (needed to set the right GENIE header locations, etc.) if needed.
#
# The top-level nutools CMakeLists.txt does this as well, but I don't think
# that we can access the result of that version check in the ubsim build.

parse_ups_version( ${GENIE_VERSION} )
if( ${VMAJ} LESS 3 )

  add_definitions(-DGENIE_PRE_R3)

  # GENIE v2 library list
  set(GENIE_LIB_LIST  ${GALGORITHM}
                    ${GBARYONRESONANCE}
                    ${GBASE}
                    ${GBODEKYANG}
                    ${GCHARM}
                    ${GCOH}
                    ${GDFRC}
                    ${GDIS}
                    ${GCROSSSECTIONS}
                    ${GDECAY}
                    ${GELAS}
                    ${GELFF}
                    ${GHEP}
                    ${GEVGCORE}
                    ${GEVGMODULES}
                    ${GEVGDRIVERS}
                    ${GGIBUU}
                    ${GHADRONTRANSP}
                    ${GFRAGMENTATION}
                    ${GINTERACTION}
                    ${GLLEWELLYNSMITH}
                    ${GMEC}
                    ${GMESSENGER}
                    ${GNUGAMMA}
                    ${GNUE}
                    ${GNTUPLE}
                    ${GNUCLEAR}
                    ${GNUMERICAL}
                    ${GQPM}
                    ${GPDG}
                    ${GPDF}
                    ${GQEL}
                    ${GRES}
                    ${GREGISTRY}
                    ${GREINSEHGAL}
                    ${GUTILS}
                    ${GGEO}
                    ${GFLUXDRIVERS}
                    ${GMUELOSS}
                    ${GREWEIGHT}
                    ${GNUCLEONDECAY} )

  #MESSAGE("--UBSIM-- GENIE_VERSION=${VMAJ} ${VMIN} -DGENIE_PRE_R3")
else()

  # GENIE v3 library list
  set (GENIE_LIB_LIST ${GFWMSG}
                      ${GFWREG}
                      ${GFWALG}
                      ${GFWINT}
                      ${GFWGHEP}
                      ${GFWNUM}
                      ${GSL}  # FWNUM relies on GSL
                      ${GFWUTL}
                      ${GFWPARDAT}
                      ${GFWEG}
                      ${GFWNTP}
                      ${GPHXSIG}
                      ${GPHPDF}
                      ${GPHNUCLST}
                      ${GPHCMN}
                      ${GPHDCY}
                      ${GPHHADTRANSP}
                      ${GPHHADNZ}
                      ${GPHDEEX}
                      ${GPHAMNGXS}
                      ${GPHAMNGEG}
                      ${GPHCHMXS}
                      ${GPHCOHXS}
                      ${GPHCOHEG}
                      ${GPHDISXS}
                      ${GPHDISEG}
                      ${GPHDFRCXS}
                      ${GPHDFRCEG}
                      ${GPHGLWRESXS}
                      ${GPHGLWRESEG}
                      ${GPHIBDXS}
                      ${GPHIBDEG}
                      ${GPHMNUCXS}
                      ${GPHMNUCEG}
                      ${GPHMEL}
                      ${GPHNUELXS}
                      ${GPHNUELEG}
                      ${GPHQELXS}
                      ${GPHQELEG}
                      ${GPHRESXS}
                      ${GPHRESEG}
                      ${GPHSTRXS}
                      ${GPHSTREG}
                      ${GPHNDCY}
                      ${GTLGEO}
                      ${GTLFLX}
                      ${GRWFWK}
                      ${GRWIO}
                      ${GRWCLC})

  #MESSAGE("--UBSIM-- GENIE_VERSION=${VMAJ} ${VMIN} not defining -DGENIE_PRE_R3")
endif()

# Build against extra reweighting functionality in GENIE v3.0.4 MicroBooNE patch 02
add_definitions(-DGENIE_UB_PATCH)

art_make(LIBRARY_NAME     EventWeightSysYC
         LIB_LIBRARIES    larsim_EventWeight_Base
                          larsim_EventWeight_Calculators
                          UBEventWeight service
                          larsim_EventWeight_Base
                          larsim_EventWeight_Calculators
                          ubsim_EventWeight_Calculators
                          ubsim_EventWeight_Calculators_BNBPrimaryHadron
                          ${Boost_SYSTEM_LIBRARY}
                          larsim_EventWeight_Base
                          larsim_EventWeight_Calculators
                          ubsim_EventWeight_Calculators
                          ubsim_EventWeight_Calculators_BNBPrimaryHadron
                          ubsim_EventWeight_Calculators_Splines
                          ReweightBaseLib
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
                          ${GENIE_LIB_LIST}
                          ${FHICLCPP}
                          ${CLHEP}
                          ${ROOT_BASIC_LIB_LIST}
                          ${Boost_SYSTEM_LIBRARY}
                          ${G4_LIB_LIST}
                          ${LOG4CPP}
                          ${MF_MESSAGELOGGER}
         MODULE_LIBRARIES larsim_EventWeight_Base
                          larsim_EventWeight_Calculators
                          #UBEventWeight service
                          larsim_EventWeight_Base
                          larsim_EventWeight_Calculators
                          ubsim_EventWeight_Calculators
                          ubsim_EventWeight_Calculators_BNBPrimaryHadron
                          ${Boost_SYSTEM_LIBRARY}
                          larsim_EventWeight_Base
                          larsim_EventWeight_Calculators
                          ubsim_EventWeight_Calculators
                          ubsim_EventWeight_Calculators_BNBPrimaryHadron
                          ubsim_EventWeight_Calculators_Splines
                          ${MF_MESSAGELOGGER})

#simple_plugin(
#  UBEventWeight service
#  larsim_EventWeight_Base
#  larsim_EventWeight_Calculators
#  ubsim_EventWeight_Calculators
#  ubsim_EventWeight_Calculators_BNBPrimaryHadron
#  ${Boost_SYSTEM_LIBRARY}
#  larsim_EventWeight_Base
#  larsim_EventWeight_Calculators
#  ubsim_EventWeight_Calculators
#  ubsim_EventWeight_Calculators_BNBPrimaryHadron
#  ubsim_EventWeight_Calculators_Splines
#  ${MF_MESSAGELOGGER}
#)

install_headers()
install_fhicl()
install_source()

