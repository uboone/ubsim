#ifndef G4RWManagerService_h
#define G4RWManagerService_h

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "geant4reweight/src/ReweightBase/G4ReweightManager.hh"

namespace evwgh {
  class G4RWManagerService;
}

class evwgh::G4RWManagerService {

public:

  G4RWManagerService(fhicl::ParameterSet const& pset);
  G4RWManagerService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  G4ReweightManager * GetManager() const;
  const fhicl::ParameterSet & GetMaterial() const;
private:
  fhicl::ParameterSet fMaterial;  
  G4ReweightManager * fRWManager;
};

DECLARE_ART_SERVICE(evwgh::G4RWManagerService, LEGACY)

#endif