cet_build_plugin(
  ExecuteCompression art::EDProducer
  LIBRARIES
  PRIVATE
  ubsim::SNStreamSim_Algo
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

add_subdirectory(Fmwk)
add_subdirectory(Algo)

install_headers()
install_fhicl()
install_source()

