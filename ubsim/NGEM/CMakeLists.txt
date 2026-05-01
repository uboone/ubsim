
cet_add_compiler_flags( CXX -DHIDE_GENIE_LOG_XXX )

include_directories ( $ENV{LOG4CPP_INC} )
include_directories ( $ENV{DK2NUDATA_INC} )
include_directories ( $ENV{DK2NUGENIE_INC} )

# genie
cet_find_library( DK2NUDATA  NAMES dk2nuTree  PATHS ENV DK2NUDATA_LIB )
cet_find_library( DK2NUGENIE NAMES dk2nuGenie PATHS ENV DK2NUGENIE_LIB )

# Determine the major version number of GENIE based on the GENIE_VERSION
# variable set by the genie ups product. Define the GENIE_PRE_R3 preprocessor
# macro (needed to set the right GENIE header locations, etc.) if needed.
#
# The top-level nutools CMakeLists.txt does this as well, but I don't think
# that we can access the result of that version check in the larsim build.
parse_ups_version( ${GENIE_VERSION} )
include_directories ( $ENV{GENIE_INC}/GENIE )
if( ${VMAJ} LESS 3 )

  add_definitions(-DGENIE_PRE_R3)

  # GENIE v2 library list
  set(GENIE_LIB_LIST  ${GBARYONRESONANCE}
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
                    ${GEVGCORE}
                    ${GEVGMODULES}
                    ${GGIBUU}
                    ${GHADRONTRANSP}
                    ${GFRAGMENTATION}
                    ${GLLEWELLYNSMITH}
                    ${GMEC}
                    ${GNUGAMMA}
                    ${GNUE}
                    ${GNUCLEAR}
                    ${GNUMERICAL}
                    ${GQPM}
                    ${GPDG}
                    ${GPDF}
                    ${GQEL}
                    ${GRES}
                    ${GREGISTRY}
                    ${GREINSEHGAL}
                    ${GGEO}
                    ${GMUELOSS}
                    ${GREWEIGHT}
                    ${GNUCLEONDECAY} )

  #MESSAGE("--LARSIM-- GENIE_VERSION=${VMAJ} ${VMIN} -DGENIE_PRE_R3")
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

  #MESSAGE("--LARSIM-- GENIE_VERSION=${VMAJ} ${VMIN} not defining -DGENIE_PRE_R3")
endif()

#MESSAGE("--LARSIM-- GENIE ${VMAJ} GENIE_LIB_LIST=${GENIE_LIB_LIST}")

art_make( 
          MODULE_LIBRARIES 
                        lardataalg_MCDumpers
                        larsim_Simulation lardataobj_Simulation
                        nusimdata_SimulationBase
                        larcoreobj_SummaryData
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        nutools_RandomUtils_NuRandomService_service
                        nutools_EventGeneratorBase_GENIE
			#nurandom_RandomUtils_NuRandomService_service
      			#nugen_EventGeneratorBase_GENIE			
			${DK2NUDATA}
                        ${DK2NUGENIE}
			#dk2nu_Tree
      			#dk2nu_Genie
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        art_Persistency_Common
                        art_Persistency_Provenance
                        art_Utilities
			#art_root_io_TFileService_service
                        canvas
                        ${MF_MESSAGELOGGER}                        
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${CLHEP}
                        ${GENIE_LIB_LIST}
                        ${ROOT_EGPYTHIA6}    # FIXME!!! - resolving genie run time reference
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_EG}
                        ${LOG4CPP}
			#Generator_Framework_GHEP
        )

install_headers()
install_fhicl()
install_source()

