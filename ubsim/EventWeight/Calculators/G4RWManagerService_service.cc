#include "G4RWManagerService.h"

evwgh::G4RWManagerService::G4RWManagerService(
    fhicl::ParameterSet const& pset) :
    fMaterial(pset.get<fhicl::ParameterSet>("Material")),
    fRWManager(new G4ReweightManager({fMaterial})) {}

evwgh::G4RWManagerService::G4RWManagerService(
    fhicl::ParameterSet const& pset,
    art::ActivityRegistry&) : G4RWManagerService(pset) {}

G4ReweightManager * evwgh::G4RWManagerService::GetManager() const {
  return fRWManager;
}

const fhicl::ParameterSet & evwgh::G4RWManagerService::GetMaterial() const {
  return fMaterial;
}

DEFINE_ART_SERVICE(evwgh::G4RWManagerService)