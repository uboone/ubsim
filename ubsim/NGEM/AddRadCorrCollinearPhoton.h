#ifndef EVGEN_ADD_RADCORR_COLLINEAR_PHOTON_H
#define EVGEN_ADD_RADCORR_COLLINEAR_PHOTON_H

#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>

// ROOT
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"

// LArSoft
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "RadCorrCollinearPhoton.h"

namespace evgen {

  // Adds one QED radiative-correction photon, collinear with the outgoing
  // charged lepton of a CC interaction, following Tomalak et al.,
  // Phys. Rev. D 106, 093006 (2022), Eqs. (35)-(36).
  //
  // The GENIE final-state lepton is treated as the tree-level lepton "jet":
  //   E_tree = E_lep(GENIE),  E_gamma = (1 - x) E_tree,  E_lep' = x E_tree,
  // and (x, eta) are sampled from dj/deta in the region
  //   E_gamma > deltaE,  E_lep' > m_lep,  lepton-photon opening angle < maxAngleDeg.
  // The radiated lepton and photon are placed so that their summed momentum
  // stays along the original lepton direction, with the sampled opening angle
  // between them and a uniformly random azimuth around that axis. Energy is
  // conserved exactly, but the lepton + photon system has invariant mass
  //   M^2 = m_lep^2 + 2 E_gamma (E_lep' - p_lep' cos(opening angle)) > m_lep^2,
  // so its 3-momentum magnitude is smaller than the tree-level lepton's. This is
  // negligible for small opening angles, but reaches tens of MeV for wide-angle
  // photons near maxAngleDeg. The missing momentum is not given to the
  // (unmodified) hadronic system.
  //
  // Particle list: the GENIE lepton is kept at its index with its tree-level kinematics,
  // but with status 3 (decayed, as for pi0s in ManuallyDecayPi0sToTwoPhotons), and the
  // radiated lepton and the photon are appended as its status-1 daughters. Only status-1
  // particles are tracked by Geant4. MCNeutrino::Lepton() stays the tree-level lepton, so
  // GENIE event reweighting (which rebuilds the GENIE event with the final-state lepton
  // momentum from MCNeutrino::Lepton()) sees exactly the original GENIE event.
  //
  // Emission modes:
  //   forceEmission = false: a photon is added with probability P(E_tree), where P is
  //     the integral of j over the allowed region. The resulting sample is unweighted.
  //   forceEmission = true:  a photon is always added (when P > 0), and the photon's
  //     MCParticle::Weight() is set to P, the weight to apply to the event.
  //
  // verbose = true also prints a line for every interaction, not just when a photon is added.
  //
  // Returns P (0 if the event is not eligible). radiated is set true if a photon was added.
  inline double AddRadCorrCollinearPhoton(int leptonPdg, double deltaE, double maxAngleDeg,
                                          bool forceEmission, bool verbose, TRandom3& randomGen,
                                          simb::MCTruth& originalMCTruth, simb::MCTruth& newMCTruth,
                                          bool& radiated) {
    radiated = false;

    // find the outgoing CC lepton
    int lepton_index = -1;
    if (originalMCTruth.NeutrinoSet() && originalMCTruth.GetNeutrino().CCNC() == simb::kCC) {
      const int lepton_track_id = originalMCTruth.GetNeutrino().Lepton().TrackId();
      for (int i = 0; i < originalMCTruth.NParticles(); ++i) {
        const simb::MCParticle& part = originalMCTruth.GetParticle(i);
        if (part.TrackId() == lepton_track_id && std::abs(part.PdgCode()) == leptonPdg && part.StatusCode() == 1) {
          lepton_index = i;
          break;
        }
      }
    }

    const double maxAngleRad = maxAngleDeg * TMath::DegToRad();
    double prob = 0.0;
    if (lepton_index >= 0) {
      const simb::MCParticle& lepton = originalMCTruth.GetParticle(lepton_index);
      prob = radcorr::RadiationProbability(lepton.E(), lepton.Mass(), deltaE, maxAngleRad);
      if (verbose) std::cout << "Rad corr: lepton pdg " << lepton.PdgCode() << " at index " << lepton_index
                << ", E_tree = " << lepton.E() << " GeV, emission probability = " << prob << std::endl;
    } else if (verbose) {
      std::cout << "Rad corr: no CC lepton with |pdg| " << leptonPdg << " found, not adding a photon" << std::endl;
    }

    // sample the photon energy fraction x and angle variable eta
    double x_sample = 0.0;
    double eta = 0.0;
    if (prob > 0.0 && (forceEmission || randomGen.Uniform(0.0, 1.0) < prob)) {
      const simb::MCParticle& lepton = originalMCTruth.GetParticle(lepton_index);
      radiated = radcorr::SampleXEta(lepton.E(), lepton.Mass(), deltaE, maxAngleRad, randomGen, x_sample, eta);
    }

    if (!radiated) {
      for (int i = 0; i < originalMCTruth.NParticles(); ++i) {
        const simb::MCParticle& particle = originalMCTruth.GetParticle(i);
        simb::MCParticle non_const_particle = simb::MCParticle(particle);
        newMCTruth.Add(non_const_particle);
      }
    } else {
      const simb::MCParticle& lepton = originalMCTruth.GetParticle(lepton_index);
      const double E_tree = lepton.E();
      const double m_lep = lepton.Mass();
      const double x = x_sample;
      const double opening_angle = eta * m_lep / E_tree;

      const double E_lep = x * E_tree;
      const double p_lep = std::sqrt(std::max(0.0, E_lep * E_lep - m_lep * m_lep));
      const double E_gamma = (1.0 - x) * E_tree;

      // angles of the lepton (a) and photon (b) from the jet axis, on opposite sides,
      // with a + b = opening angle and balanced transverse momentum
      const double a = std::atan2(E_gamma * std::sin(opening_angle), p_lep + E_gamma * std::cos(opening_angle));
      const double b = opening_angle - a;

      const TVector3 axis = lepton.Momentum().Vect().Unit();
      const TVector3 e1 = axis.Orthogonal().Unit();
      const TVector3 e2 = axis.Cross(e1);
      const double phi = randomGen.Uniform(0.0, TMath::TwoPi());
      const TVector3 transverse = std::cos(phi) * e1 + std::sin(phi) * e2;

      const TVector3 lep_dir = std::cos(a) * axis + std::sin(a) * transverse;
      const TVector3 gamma_dir = std::cos(b) * axis - std::sin(b) * transverse;
      const TLorentzVector lep_momentum(p_lep * lep_dir, E_lep);
      const TLorentzVector gamma_momentum(E_gamma * gamma_dir, E_gamma);

      int max_track_id = 0;
      for (int i = 0; i < originalMCTruth.NParticles(); ++i) {
        max_track_id = std::max(max_track_id, originalMCTruth.GetParticle(i).TrackId());
      }
      const int radiated_lepton_track_id = max_track_id + 1;
      const int gamma_track_id = max_track_id + 2;

      for (int i = 0; i < originalMCTruth.NParticles(); ++i) {
        const simb::MCParticle& particle = originalMCTruth.GetParticle(i);
        if (i == lepton_index) {
          // tree-level lepton, unchanged except for status 3 (decayed) and its two new daughters
          // track_id, pdg, process, mother, mass, status_code
          simb::MCParticle treeLepton(particle.TrackId(), particle.PdgCode(), particle.Process(), particle.Mother(), particle.Mass(), 3);
          treeLepton.SetPolarization(particle.Polarization());
          treeLepton.SetGvtx(particle.Gvx(), particle.Gvy(), particle.Gvz(), particle.Gvt());
          treeLepton.SetRescatter(particle.Rescatter());
          treeLepton.SetWeight(particle.Weight());
          treeLepton.SetEndProcess(particle.EndProcess());
          for (int i_d = 0; i_d < particle.NumberDaughters(); ++i_d) treeLepton.AddDaughter(particle.Daughter(i_d));
          treeLepton.AddDaughter(radiated_lepton_track_id);
          treeLepton.AddDaughter(gamma_track_id);
          treeLepton.AddTrajectoryPoint(particle.Position(), particle.Momentum());
          newMCTruth.Add(treeLepton);
        } else {
          simb::MCParticle non_const_particle = simb::MCParticle(particle);
          newMCTruth.Add(non_const_particle);
        }
      }

      // the radiated lepton and the photon start at the interaction vertex
      simb::MCParticle radiatedLepton(radiated_lepton_track_id, lepton.PdgCode(), "primary", lepton.TrackId(), m_lep, 1);
      radiatedLepton.SetPolarization(lepton.Polarization());
      radiatedLepton.SetGvtx(lepton.Gvx(), lepton.Gvy(), lepton.Gvz(), lepton.Gvt());
      radiatedLepton.AddTrajectoryPoint(lepton.Position(), lep_momentum);
      newMCTruth.Add(radiatedLepton);

      simb::MCParticle gamma(gamma_track_id, 22, "primary", lepton.TrackId(), 0.0, 1);
      gamma.SetGvtx(lepton.Gvx(), lepton.Gvy(), lepton.Gvz(), lepton.Gvt());
      gamma.SetWeight(forceEmission ? prob : 1.0);
      gamma.AddTrajectoryPoint(lepton.Position(), gamma_momentum);
      newMCTruth.Add(gamma);

      std::cout << "Rad corr: added photon, E_tree = " << E_tree << " GeV, emission probability = " << prob << ", x = " << x << ", eta = " << eta
                << ", E_gamma = " << E_gamma << " GeV, opening angle = " << opening_angle * TMath::RadToDeg() << " deg" << std::endl;
      TLorentzVector diff = lepton.Momentum() - lep_momentum - gamma_momentum;
      std::cout << "Rad corr: (tree lepton - radiated lepton - photon) 4-momentum: ("
                << diff.Px() << ", " << diff.Py() << ", " << diff.Pz() << ", " << diff.E() << ")" << std::endl;
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

    return prob;
  }

  // In-place convenience wrapper
  inline double AddRadCorrCollinearPhoton(int leptonPdg, double deltaE, double maxAngleDeg,
                                          bool forceEmission, bool verbose, TRandom3& randomGen,
                                          simb::MCTruth& mcTruth, bool& radiated) {
    simb::MCTruth newTruth;
    double prob = AddRadCorrCollinearPhoton(leptonPdg, deltaE, maxAngleDeg, forceEmission, verbose, randomGen, mcTruth, newTruth, radiated);
    mcTruth = newTruth;
    return prob;
  }

} // namespace evgen

#endif // EVGEN_ADD_RADCORR_COLLINEAR_PHOTON_H
