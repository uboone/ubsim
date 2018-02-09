#include_directories ( $ENV{GENIE_INC}/GENIE )

art_make(LIBRARY_NAME  uboone_EventWeight_Calculators_BNBPrimaryHadron
         LIB_LIBRARIES
                       larsim_EventWeight_Base
                       lardata_Utilities
                       ${ART_FRAMEWORK_CORE}
                       ${ART_FRAMEWORK_PRINCIPAL}
                       ${ART_FRAMEWORK_SERVICES_REGISTRY}
                       ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                       ${CLHEP}
                       ${ROOT_BASIC_LIB_LIST}
                       ${Boost_SYSTEM_LIBRARY}
        )

install_headers()
install_fhicl()
install_source()


