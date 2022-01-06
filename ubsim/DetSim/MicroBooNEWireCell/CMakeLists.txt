# Build the libWireCellubsim WCT plugin library.

include_directories($ENV{EIGEN_INC})

link_directories( $ENV{WIRECELL_LIB})

file(GLOB wcls_components_sources "*.cxx")

# This is what's needed for building the WireCellubsim WCT component
# plugin library.  
art_make_library(
  LIBRARY_NAME WireCellubsim
  SOURCE ${wcls_components_sources}
  LIBRARIES
        art::Framework_Core
        art::Framework_Principal
        art::Framework_Services_Registry
        art::Persistency_Common art_Persistency_Provenance
        art::Utilities
        canvas::canvas cetlib
        cetlib_except::cetlib_except
        
	jsoncpp
	lardata_Utilities
	lardataobj_RawData
	lardataobj_RecoBase
	lardataobj_Simulation
        larcorealg_Geometry
	larevt_CalibrationDBI_IOVData
	larevt_CalibrationDBI_Providers
        #${WIRECELL_LIBS}
	WireCellApps
	WireCellSigProc
	WireCellIface
	WireCellUtil
	WireCellGen
	WireCellAux
        fhiclcpp::fhiclcpp
        cetlib::cetlib
        cetlib_except::cetlib_except
)
