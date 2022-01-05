if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
cet_add_compiler_flags(CXX -Wno-return-type-c-linkage )
endif()

#cet_find_library( CRY NAMES CRY PATHS ENV CRY_LIB NO_DEFAULT_PATH )
include_directories ( $ENV{CRYHOME}/src )

art_make( BASENAME_ONLY

          MODULE_LIBRARIES larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
                        lardata_Utilities
			lardataobj_RecoBase
			lardataobj_AnalysisBase
                        larevt_Filters
			larcoreobj_SummaryData
			larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
                        lardataobj_RawData
                        nurandom::RandomUtils_NuRandomService_service
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
                        art::Framework_Services_Optional_RandomNumberGenerator_service
                        art_root_io::TFileService_service
                        art::Persistency_Common canvas::canvas
                        art::Persistency_Provenance canvas::canvas
                        art::Utilities canvas::canvas
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib
                        cetlib_except::cetlib_except
                        CLHEP::CLHEP
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml
                        ROOT::EG
                        CRY::CRY

        )

install_headers()
install_fhicl()
install_source()

add_subdirectory(jobs)
