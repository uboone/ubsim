#ifndef EVGEN_MANUALLYDECAYPI0S_TO_TWO_PHOTONS_H
#define EVGEN_MANUALLYDECAYPI0S_TO_TWO_PHOTONS_H

#include <vector>
#include <iostream>
#include <cmath>

// ROOT
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"

// LArSoft
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

namespace evgen {

  // note that this does not consider rarer decays of pi0s
  inline void ManuallyDecayPi0sToTwoPhotons(simb::MCTruth& originalMCTruth, simb::MCTruth& newMCTruth) {
    for (int i = 0; i < originalMCTruth.NParticles(); ++i) {
      TRandom3 randomGen;
      randomGen.SetSeed(0);

      int pdgCode = originalMCTruth.GetParticle(i).PdgCode();
      int statusCode = originalMCTruth.GetParticle(i).StatusCode();
      if (pdgCode == 111 && statusCode == 1) { // pi0, GENIE status code 1 = kIstStableFinalState
	std::cout << "Decaying pi0 at index " << i << std::endl;
        const simb::MCParticle& pi0 = originalMCTruth.GetParticle(i);

        TLorentzVector pi0_position(pi0.Vx(), pi0.Vy(), pi0.Vz(), pi0.T());
        TLorentzVector pi0_momentum = pi0.Momentum();

        // 3 = kIStDecayedState from GENIE
	simb::MCParticle newPi0(pi0.TrackId(), pi0.PdgCode(), pi0.Process(), pi0.Mother(), pi0.Mass(), 3);
        newPi0.AddTrajectoryPoint(pi0_position, pi0_momentum);

        // track_id, pdg, process, mother, mass, status_code
	simb::MCParticle gamma1(979797971, 22, "primary", pi0.TrackId(), 0.0, 1);
	simb::MCParticle gamma2(979797972, 22, "primary", pi0.TrackId(), 0.0, 1);

        // GENIE vertices relative to nucleus
        gamma1.SetGvtx(pi0.Gvx(), pi0.Gvy(), pi0.Gvz(), pi0.Gvt());
        gamma2.SetGvtx(pi0.Gvx(), pi0.Gvy(), pi0.Gvz(), pi0.Gvt());

        TVector3 pi0_boost_vector = pi0_momentum.BoostVector();

        // Calculate decay momenta
        double pi0_mass = pi0.Mass();
        double theta = randomGen.Uniform(0, 2 * M_PI);
        double phi = acos(1 - 2 * randomGen.Uniform(0, 1));
        double pE = pi0_mass / 2;
        double px = pE * sin(phi) * cos(theta);
        double py = pE * sin(phi) * sin(theta);
        double pz = pE * cos(phi);
        TLorentzVector rest_frame_momentum_1(px, py, pz, pE);
        TLorentzVector rest_frame_momentum_2(-px, -py, -pz, pE);

        // transform the momenta to the lab frame using a Lorentz boost
        TLorentzVector lab_frame_momentum_1 = rest_frame_momentum_1;
        TLorentzVector lab_frame_momentum_2 = rest_frame_momentum_2;
        lab_frame_momentum_1.Boost(pi0_boost_vector);
        lab_frame_momentum_2.Boost(pi0_boost_vector);

        gamma1.AddTrajectoryPoint(pi0_position, lab_frame_momentum_1);
        gamma2.AddTrajectoryPoint(pi0_position, lab_frame_momentum_2);

	std::cout << "Conservation check:" << std::endl;
        TLorentzVector sum = lab_frame_momentum_1 + lab_frame_momentum_2;
	std::cout << "4-momentum difference (should be ~0): "
                  << (pi0_momentum - sum).Px() << ", "
                  << (pi0_momentum - sum).Py() << ", "
                  << (pi0_momentum - sum).Pz() << ", "
                  << (pi0_momentum - sum).E() << std::endl;

        newMCTruth.Add(newPi0);
        newMCTruth.Add(gamma1);
        newMCTruth.Add(gamma2);
      } else {
        // Copy over non-pi0 particles
        const simb::MCParticle& particle = originalMCTruth.GetParticle(i);
	simb::MCParticle non_const_particle = simb::MCParticle(particle);
        newMCTruth.Add(non_const_particle);
      }
    }

    // Copy over the neutrino information
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
  inline void ManuallyDecayPi0sToTwoPhotons(simb::MCTruth& mcTruth) {
    simb::MCTruth newTruth;
    ManuallyDecayPi0sToTwoPhotons(mcTruth, newTruth);
    mcTruth = newTruth;
  }

} // namespace evgen

#endif // EVGEN_MANUALLYDECAYPI0S_TO_TWO_PHOTONS_H
