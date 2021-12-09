
art_make(
  LIBRARY_NAME ubsim_EventWeight_Calculators_MultiVariateRNG
  LIB_LIBRARIES
    ${FHICLCPP}
    ${CLHEP}
    ${ROOT_BASIC_LIB_LIST}
    ${Boost_SYSTEM_LIBRARY}
    cetlib cetlib_except
)

install_headers()
install_source()
