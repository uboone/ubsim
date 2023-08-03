cet_make_library(
  SOURCE
  Channel.cxx
  Flux.cxx
  FourMomentum.cxx
  Observables.cxx
  SterileNeutrino.cxx
  LIBRARIES
  PUBLIC
  CLHEP::Random
)

install_headers()
install_fhicl()
install_source()
