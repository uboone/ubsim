cet_build_plugin(
  UBEventWeight art::service
  LIBRARIES
  PRIVATE
  larsim::EventWeight_Base
)

install_headers()
install_fhicl()
install_source()
