#include_directories ( $ENV{GENIE_INC}/GENIE )

art_make(LIBRARY_NAME  ubsim_EventWeight_Calculators_BNBPrimaryHadron
         LIB_LIBRARIES
                       larsim_EventWeight_Base
                       lardata_Utilities
                       cetlib::cetlib
                       cetlib_except::cetlib_except
                       art::Framework_Core
                       art::Framework_Principal
                       art::Framework_Services_Registry
                       art_root_io::tfile_support ROOT::Core
                       art::Framework_Services_Optional_RandomNumberGenerator_service
                       CLHEP::CLHEP

                       nugen::NuReweight_art
                       nugen::NuReweight
                       larsim_Simulation
                       nurandom::RandomUtils_NuRandomService_service
                       nusimdata::SimulationBase
                       art::Persistency_Common canvas::canvas canvas::canvas
                       art::Persistency_Provenance canvas::canvas canvas::canvas
                       art::Utilities canvas::canvas canvas::canvas
                       fhiclcpp::fhiclcpp
        )

install_headers()
install_fhicl()
install_source()


