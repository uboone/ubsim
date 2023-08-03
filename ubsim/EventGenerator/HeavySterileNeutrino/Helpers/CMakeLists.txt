cet_make_library(
  SOURCE
  Helper.cxx
  LIBRARIES
  PUBLIC
  ubsim::EventGenerator_HeavySterileNeutrino_DataObjects
)

install_headers()
install_fhicl()
install_source()

