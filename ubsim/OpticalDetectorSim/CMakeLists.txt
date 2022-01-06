art_make( EXCLUDE UBOpticalChConfig_service.cc OpticalDRAMReadout_module.cc
	  #NO_PLUGINS
	  BASENAME_ONLY
	  LIBRARY_NAME  uboonecode_uboone_OpticalDetectorSim
          LIB_LIBRARIES UBOpticalChConfig_service
			larcoreobj_SummaryData
			lardataobj_RawData
			lardata_Utilities
			larcorealg_Geometry
                        nurandom::RandomUtils_NuRandomService_service
                        nusimdata::SimulationBase
                        art::Framework_Core
			art::Framework_Principal
			art::Framework_Services_Registry
			art_root_io::tfile_support ROOT::Core
			art_root_io::TFileService_service
                        art::Framework_Services_Optional_RandomNumberGenerator_service
                        art::Persistency_Common
			art::Utilities
                        canvas::canvas
			messagefacility::MF_MessageLogger
	                fhiclcpp::fhiclcpp
			cetlib::cetlib
                        cetlib_except::cetlib_except

          MODULE_LIBRARIES larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
			lardataobj_OpticalDetectorData
			larcoreobj_SummaryData
	       		ubevt_Database
			UBOpticalChConfig_service
			ubcore::UBOpReadoutMap_service
			uboonecode_uboone_OpticalDetectorSim
			lardata_Utilities
			ubcore::UBooNEGeometryHelper_service
			ubcore_Geometry
                        nusimdata::SimulationBase
                        art::Framework_Core
			art::Framework_Principal
			art::Framework_Services_Registry
			art_root_io::tfile_support ROOT::Core
			art_root_io::TFileService_service
                        art::Framework_Services_Optional_RandomNumberGenerator_service
                        art::Persistency_Common
                        art::Persistency_Provenance
                        canvas::canvas
                        art::Utilities
			messagefacility::MF_MessageLogger
	                fhiclcpp::fhiclcpp
                        cetlib::cetlib
                        cetlib_except::cetlib_except
                        ROOT::Tree
         )

simple_plugin(UBOpticalChConfig "service"
	      ubcore_Geometry
	      ubcore::UBOpReadoutMap_service
	      larcore_Geometry_Geometry_service
	      larcorealg_Geometry
              art::Utilities
              canvas::canvas
	      messagefacility::MF_MessageLogger
	      fhiclcpp::fhiclcpp
              cetlib::cetlib
              cetlib_except::cetlib_except

              BASENAME_ONLY
)

simple_plugin(OpticalDRAMReadout "module"
                  larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
                  lardata_DetectorInfoServices_DetectorClocksServiceStandard_service # BUG note that we link to a specific service implementation
                  lardataobj_OpticalDetectorData
                  UBOpticalChConfig_service
                  ubcore::UBOpReadoutMap_service
                  uboonecode_uboone_OpticalDetectorSim
                  lardata_Utilities
                  ubcore::UBooNEGeometryHelper_service
                  ubcore_Geometry
                  nusimdata::SimulationBase
                  art::Framework_Core
                  art::Framework_Principal
                  art::Framework_Services_Registry
                  art_root_io::tfile_support ROOT::Core
                  art_root_io::TFileService_service
                  art::Framework_Services_Optional_RandomNumberGenerator_service
                  art::Persistency_Common
                  art::Persistency_Provenance
                  art::Utilities
                  canvas::canvas
                  messagefacility::MF_MessageLogger
                  fhiclcpp::fhiclcpp
                  cetlib::cetlib
                  cetlib_except::cetlib_except

              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
