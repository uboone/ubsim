if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
cet_add_compiler_flags(CXX -Wno-return-type-c-linkage )
endif()

art_make( BASENAME_ONLY
          EXCLUDE SimWireMicroBooNE_module.cc
          MODULE_LIBRARIES larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
			lardata_Utilities_LArFFT_service
                        lardata_Utilities
                        larevt_Filters
                        lardataalg_DetectorInfo
                        lardataobj_RawData
                        SignalShapingServiceMicroBooNE_service
			larevt_CalibrationDBI_Providers
                        nurandom_RandomUtils_NuRandomService_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}

                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}

        )

simple_plugin(SimWireMicroBooNE "module"
          BASENAME_ONLY
              larcorealg_Geometry
              larcore_Geometry_Geometry_service
              larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
	      lardata_Utilities_LArFFT_service
              lardata_Utilities
              larevt_Filters
              lardataalg_DetectorInfo
              lardataobj_RawData
              lardata_DetectorInfoServices_DetectorClocksServiceStandard_service # BUG note that we link to a specific service implementation
              SignalShapingServiceMicroBooNE_service
              larevt_CalibrationDBI_Providers
              nurandom_RandomUtils_NuRandomService_service
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              art_Persistency_Common canvas
              art_Persistency_Provenance canvas
              art_Utilities canvas
              ${MF_MESSAGELOGGER}

              ${FHICLCPP}
              cetlib cetlib_except
              ${CLHEP}
              ${ROOT_GEOM}
              ${ROOT_XMLIO}
              ${ROOT_GDML}

        )

# install_headers()
install_fhicl()
install_source()
