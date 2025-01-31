cet_build_plugin(
  UBEventWeight art::service
  LIBRARIES
  PRIVATE
  larsim::EventWeight_Base
  larsim::EventWeight_Calculators
  ubsim::EventWeight_Calculators
  ubsim::EventWeight_Calculators_BNBPrimaryHadron
  ubsim::EventWeight_Calculators_Splines
)

install_headers()
install_fhicl()
install_source()
