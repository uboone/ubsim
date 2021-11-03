
art_make(
	  LIB_LIBRARIES

			${ART_FRAMEWORK_CORE}
			art_Utilities canvas
			${MF_MESSAGELOGGER}

			cetlib cetlib_except

          MODULE_LIBRARIES
			ubsim_SNStreamSim_Fmwk
			ubsim_SNStreamSim_Algo
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			larsim_MCCheater_ParticleInventoryService_service
			larsim_MCCheater_BackTrackerService_service
			larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
			lardataobj_RecoBase
			lardataobj_AnalysisBase
			lardata_RecoObjects
			lardata_Utilities
			larevt_Filters
			nusimdata::SimulationBase
			${LARDATA_LIB}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
                        canvas
			${MF_MESSAGELOGGER}
	                ${FHICLCPP}
			cetlib cetlib_except
			${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ROOT::Tree
        )

add_subdirectory(Fmwk)
add_subdirectory(Algo)

install_headers()
install_fhicl()
install_source()

