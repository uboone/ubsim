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
        art_Framework_Core art_Framework_Principal art_Framework_Services_Registry
        art_Persistency_Common art_Persistency_Provenance
        art_Utilities
        canvas cetlib cetlib_except
        ${Boost_SYSTEM_LIBRARY}
        ${ROOT_BASIC_LIB_LIST}
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
        ${FHICLCPP}
        cetlib cetlib_except
)
