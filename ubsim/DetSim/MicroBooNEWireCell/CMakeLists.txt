cet_make_library(
  LIBRARY_NAME WireCellubsim       # Wire cell toolkit needs this name.
  SOURCE
  ReweightedDepoTransform.cxx
  UbDrifter.cxx
  LIBRARIES
  PUBLIC
  ubevt::Database_TPCEnergyCalib
  ubobj::DataOverlay
  larwirecell::IArtEventVisitor
  larevt::ElectronicsCalibService
  lardata::DetectorPropertiesService
  larcore::ServiceUtil
  WireCell::Gen
  art::Framework_Services_Registry
  art::Framework_Principal
  ROOT::Hist
)
