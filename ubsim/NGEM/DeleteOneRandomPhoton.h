#ifndef EVGEN_DELETE_ONE_RANDOM_PHOTON_H
#define EVGEN_DELETE_ONE_RANDOM_PHOTON_H

#include <vector>
#include <iostream>

// ROOT
#include "TRandom3.h"

// LArSoft
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

namespace evgen {

  inline void DeleteOneRandomPhoton(simb::MCTruth& originalMCTruth, simb::MCTruth& newMCTruth) {
    TRandom3 randomGen;
    randomGen.SetSeed(0);
    std::vector<int> photon_indices;
    for (int i = 0; i < originalMCTruth.NParticles(); ++i) {
      int pdgCode = originalMCTruth.GetParticle(i).PdgCode();
      int statusCode = originalMCTruth.GetParticle(i).StatusCode();
      if (pdgCode == 22 && statusCode == 1) {
        photon_indices.push_back(i);
      }
    }
    int num_photons = photon_indices.size();
    std::cout << "Number of primary photons: " << num_photons << std::endl;
    if (num_photons == 0) {
      std::cout << "No photons to delete" << std::endl;
      for (int i = 0; i < originalMCTruth.NParticles(); ++i) {
        const simb::MCParticle& particle = originalMCTruth.GetParticle(i);
	simb::MCParticle non_const_particle = simb::MCParticle(particle);
        newMCTruth.Add(non_const_particle);
      }
    } else {
      int index_to_delete = photon_indices[randomGen.Integer(num_photons)];
      for (int i = 0; i < originalMCTruth.NParticles(); ++i) {
        if (i == index_to_delete) {
	  std::cout << "Deleted photon at index: " << index_to_delete << std::endl;
        } else { // copying over particles that aren't getting deleted
          const simb::MCParticle& particle = originalMCTruth.GetParticle(i);
	  simb::MCParticle non_const_particle = simb::MCParticle(particle);
          newMCTruth.Add(non_const_particle);
        }
      }
    }
    simb::MCNeutrino neutrino = originalMCTruth.GetNeutrino();
    int CCNC = neutrino.CCNC();
    int mode = neutrino.Mode();
    int interactionType = neutrino.InteractionType();
    int target = neutrino.Target();
    int nucleon = neutrino.HitNuc();
    int quark = neutrino.HitQuark();
    double w = neutrino.W();
    double x = neutrino.X();
    double y = neutrino.Y();
    double qsqr = neutrino.QSqr();
    newMCTruth.SetNeutrino(CCNC, mode, interactionType, target, nucleon, quark, w, x, y, qsqr);
    newMCTruth.SetOrigin(originalMCTruth.Origin());
  }

  // In-place convenience wrapper
  inline void DeleteOneRandomPhoton(simb::MCTruth& mcTruth) {
    simb::MCTruth newTruth;
    DeleteOneRandomPhoton(mcTruth, newTruth);
    mcTruth = newTruth;
  }

} // namespace evgen

#endif // EVGEN_DELETE_ONE_RANDOM_PHOTON_H
