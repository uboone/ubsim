cet_make_library(
  SOURCE
  AlgorithmFactory.cxx
  MicrobooneFirmware.cxx
  LIBRARIES
  PUBLIC
  ubsim::SNStreamSim_Fmwk
)

install_headers()
install_fhicl()
install_source()

