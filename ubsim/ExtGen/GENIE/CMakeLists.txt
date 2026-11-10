# GENIE algorithms that let GENIE's flux and geometry drivers choose
# interactions using the total cross sections of an external generator
# (GiBUU, ...). This is a plain GENIE plugin library, not an art module:
# GENIE instantiates the classes by name through their ROOT dictionary.

cet_add_compiler_flags( CXX -DHIDE_GENIE_LOG_XXX )

art_make_library( LIBRARY_NAME ubsim_ExtGen_GENIE
                  SOURCE ExtGenPXSec.cxx
                         ExtGenInteractionListGenerator.cxx
                  LIBRARIES
                    log4cpp::log4cpp
                    GENIE::GFwUtl
                    GENIE::GFwParDat
                    GENIE::GFwAlg
                    GENIE::GFwInt
                    GENIE::GFwEG
                    GENIE::GFwMsg
                    GENIE::GFwReg
                    ROOT::Tree
                    ROOT::EG
                )

# Produces libubsim_ExtGen_GENIE_dict.so and its .rootmap, which is how
# GENIE (gROOT->GetClass) finds the classes at run time.
build_dictionary( NO_CHECK_CLASS_VERSION
                  DICTIONARY_LIBRARIES
                    ubsim::ExtGen_GENIE
                )

install_headers()
install_source()
