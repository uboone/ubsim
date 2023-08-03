cet_make_library(
  SOURCE
  CompressionAlgoBase.cxx
  LIBRARIES
  PUBLIC
  fhiclcpp::fhiclcpp
)

install_headers()
install_fhicl()
install_source()

