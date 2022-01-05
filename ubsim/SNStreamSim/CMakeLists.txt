
art_make(
	  LIB_LIBRARIES

			art::Framework_Core
			art::Utilities canvas::canvas
			messagefacility::MF_MessageLogger
			cetlib::cetlib
                        cetlib_except::cetlib_except

          MODULE_LIBRARIES
			ubsim_SNStreamSim_Fmwk
			ubsim_SNStreamSim_Algo
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			larsim_MCCheater_ParticleInventoryService_service
			larsim_MCCheater_BackTrackerService_service
			larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
			lardataobj_RecoBase
			lardataobj_AnalysisBase
			lardata_RecoObjects
			lardata_Utilities
			larevt_Filters
			nusimdata::SimulationBase
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
                        ROOT::Tree
        )

add_subdirectory(Fmwk)
add_subdirectory(Algo)

install_headers()
install_fhicl()
install_source()

