if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
cet_add_compiler_flags(CXX -Wno-return-type-c-linkage )
endif()

art_make( BASENAME_ONLY
          EXCLUDE SimWireMicroBooNE_module.cc
          MODULE_LIBRARIES larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
			lardata_Utilities_LArFFT_service
                        lardata_Utilities
                        larevt_Filters
                        lardataalg_DetectorInfo
                        lardataobj_RawData
                        ubevt::SignalShapingServiceMicroBooNE_service
			larevt_CalibrationDBI_Providers
                        nurandom::RandomUtils_NuRandomService_service
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
                        art::Framework_Services_Optional_RandomNumberGenerator_service
                        art_root_io::TFileService_service
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities
                        canvas::canvas
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib
                        cetlib_except::cetlib_except
                        CLHEP::CLHEP
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml
        )

simple_plugin(SimWireMicroBooNE "module"
          BASENAME_ONLY
              larcorealg_Geometry
              larcore_Geometry_Geometry_service
              larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
	      lardata_Utilities_LArFFT_service
              lardata_Utilities
              larevt_Filters
              lardataalg_DetectorInfo
              lardataobj_RawData
              lardata_DetectorInfoServices_DetectorClocksServiceStandard_service # BUG note that we link to a specific service implementation
              ubevt::SignalShapingServiceMicroBooNE_service
              larevt_CalibrationDBI_Providers
              nurandom::RandomUtils_NuRandomService_service
              art::Framework_Core
              art::Framework_Principal
              art::Framework_Services_Registry
              art_root_io::tfile_support ROOT::Core
              art::Framework_Services_Optional_RandomNumberGenerator_service
              art_root_io::TFileService_service
              art::Persistency_Common
              art::Persistency_Provenance
              art::Utilities
              canvas::canvas
              messagefacility::MF_MessageLogger
              fhiclcpp::fhiclcpp
              cetlib::cetlib
              cetlib_except::cetlib_except
              CLHEP::CLHEP
              ROOT::Geom
              ROOT::XMLIO
              ROOT::Gdml
        )

# install_headers()
install_fhicl()
install_source()

add_subdirectory(MicroBooNEWireCell)
