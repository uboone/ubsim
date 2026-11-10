# art modules of the ExtGen chain, and the library ubsim_ExtGen_Modules
# (ExtGenIO: request/served-event files and frame rotation) they share.

art_make_library( LIBRARY_NAME ubsim_ExtGen_Modules
                  SOURCE ExtGenIO.cxx
                  LIBRARIES
                    cetlib::cetlib
                    cetlib_except::cetlib_except
                    ROOT::Tree
                    ROOT::EG
                )

cet_build_plugin( ExtGenClosureCheck art::EDAnalyzer
                  LIBRARIES
                  PRIVATE
                    larcoreobj::SummaryData
                    lardataobj::Simulation
                    nusimdata::SimulationBase
                )

cet_build_plugin( ExtGenMixCheck art::EDAnalyzer
                  LIBRARIES
                  PRIVATE
                    larcoreobj::SummaryData
                    lardataobj::Simulation
                    nusimdata::SimulationBase
                )

cet_build_plugin( ExtGenRequestDumper art::EDAnalyzer
                  LIBRARIES
                  PRIVATE
                    larcoreobj::SummaryData
                    lardataobj::Simulation
                    nusimdata::SimulationBase
                )

cet_build_plugin( ExtGenTruthServer art::EDAnalyzer
                  LIBRARIES
                  PRIVATE
                    ubsim::ExtGen_Modules
                    larcoreobj::SummaryData
                    lardataobj::Simulation
                    nusimdata::SimulationBase
                )

cet_build_plugin( ExtGenMerge art::EDProducer
                  LIBRARIES
                  PRIVATE
                    ubsim::ExtGen_Modules
                    lardata::Utilities
                    larcoreobj::SummaryData
                    lardataobj::Simulation
                    nusimdata::SimulationBase
                    nurandom::RandomUtils_NuRandomService_service
                    dk2nu::Tree
                )

install_headers()
install_source()
