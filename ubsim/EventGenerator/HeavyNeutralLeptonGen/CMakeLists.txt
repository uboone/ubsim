find_ups_product( ifdhc )
find_ups_product( ifdh_art )

# ifdhc
cet_find_library( IFDH NAMES ifdh PATHS ENV IFDHC_LIB NO_DEFAULT_PATH )
include_directories( $ENV{IFDHC_FQ_DIR}/inc )
# ifdhc_art  # as a service
cet_find_library( IFDH_SERVICE NAMES IFDH_service PATHS ENV IFDH_ART_LIB NO_DEFAULT_PATH )
include_directories( $ENV{IFDH_ART_FQ_DIR}/inc )

#include_directories(${IFDH_ART_INC})
include_directories(${IFDHC_FQ_DIR}/inc)

find_ups_product(genie)

art_make(BASENAME_ONLY
          LIBRARY_NAME uboonecode_HeavyNeutralLeptonGen
          LIB_LIBRARIES
            larcoreobj_SummaryData
            larcorealg_Geometry
            larcore_Geometry_Geometry_service
            nurandom::RandomUtils_NuRandomService_service
            lardataobj_Simulation
            nusimdata::SimulationBase
            art::Framework_Core
            art::Framework_Principal
            art::Framework_Services_Registry

            art::Framework_Services_Optional_RandomNumberGenerator_service
            art::Persistency_Common
            art::Persistency_Provenance
            art::Utilities
            canvas::canvas
            cetlib::cetlib
            cetlib_except::cetlib_except
            messagefacility::MF_MessageLogger

            fhiclcpp::fhiclcpp
            CLHEP::CLHEP
            SQLite::SQLite3
            ROOT::EG
            ROOT::Tree
            ROOT::Geom
            ifdh::ifdh

            dk2nu::Tree
            GENIE::GFwMsg
            GENIE::GFwParDat
            log4cpp::log4cpp
      MODULE_LIBRARIES
            uboonecode_HeavyNeutralLeptonGen
            larcoreobj_SummaryData
            larcorealg_Geometry
            larcore_Geometry_Geometry_service
            nurandom::RandomUtils_NuRandomService_service
            lardataobj_Simulation
            nusimdata::SimulationBase
            art::Framework_Core
            art::Framework_Principal
            art::Framework_Services_Registry
            art_root_io::tfile_support
            art_root_io::TFileService_service
            art::Persistency_Common
            art::Persistency_Provenance
            art::Utilities
            canvas::canvas
            cetlib::cetlib
            cetlib_except::cetlib_except
            messagefacility::MF_MessageLogger

            fhiclcpp::fhiclcpp
            CLHEP::CLHEP
            SQLite::SQLite3
            ROOT::Core
            ROOT::EG
            ROOT::Geom
            ifdh::ifdh
            dk2nu::Tree
            GENIE::GFwMsg
            GENIE::GFwParDat
            log4cpp::log4cpp
            )

add_subdirectory(job)

install_headers()
install_fhicl()
install_source()

