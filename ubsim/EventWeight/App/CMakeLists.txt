simple_plugin(
  UBEventWeight service
  larsim_EventWeight_Base
  larsim_EventWeight_Calculators
  ubsim_EventWeight_Calculators
  ubsim_EventWeight_Calculators_BNBPrimaryHadron
  larsim_EventWeight_Base
  larsim_EventWeight_Calculators
  ubsim_EventWeight_Calculators
  ubsim_EventWeight_Calculators_BNBPrimaryHadron
  ubsim_EventWeight_Calculators_Splines
  ${MF_MESSAGELOGGER}
)

install_headers()
install_fhicl()
install_source()

