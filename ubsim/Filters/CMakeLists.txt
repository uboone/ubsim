if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
cet_add_compiler_flags(CXX -Wno-return-type-c-linkage )
endif()

#cet_find_library( CRY NAMES CRY PATHS ENV CRY_LIB NO_DEFAULT_PATH )
include_directories ( $ENV{CRYHOME}/src )

art_make( BASENAME_ONLY

          MODULE_LIBRARIES larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
                        lardata_Utilities
			lardataobj_RecoBase
			lardataobj_AnalysisBase
                        larevt_Filters
			larcoreobj_SummaryData
			larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
                        lardataobj_RawData
                        nurandom_RandomUtils_NuRandomService_service
			nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
			${LARRECO_LIB}
			${LARDATA_LIB}
                        ${MF_MESSAGELOGGER}                        
                        ${FHICLCPP}
                        cetlib cetlib_except
                        CLHEP::CLHEP
			${ROOT_PHYSICS}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${ROOT_BASIC_LIB_LIST}
			${ROOT_CORE}
			${ROOT_EG}
			CRY
        )

install_headers()
install_fhicl()
install_source()

add_subdirectory(jobs)
