#ifndef EVGEN_MODIFYPARTICLE_H
#define EVGEN_MODIFYPARTICLE_H

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

  // Utility: sample from PDF defined by bin edges and per-bin probabilities (weights)
  inline double sampleFromBinned(const std::vector<double>& binEdges,
                                 const std::vector<double>& binProbs,
                                 TRandom3& randomGen) {
    const std::size_t numBins = binProbs.size();

    if (binEdges.size() != numBins + 1) {
      throw std::runtime_error("sampleFromBinned: Number of bin edges must be equal to number of weights + 1");
    }
    if (numBins == 0) {
      throw std::runtime_error("sampleFromBinned: Number of bins must be greater than 0");
    }

    // check that the bin edges are in ascending order
    for (std::size_t i = 0; i < numBins; ++i) {
      if (binEdges[i] >= binEdges[i + 1]) {
        throw std::runtime_error("sampleFromBinned: Bin edges must be in ascending order");
      }
    }

    // check that the bin weights are non-negative
    for (std::size_t i = 0; i < numBins; ++i) {
      if (binProbs[i] < 0.0) {
        throw std::runtime_error("sampleFromBinned: Bin weight must be greater than or equal to 0");
      }
    }

    // Build cumulative probabilities (weights are per-bin probabilities; uniform within each bin)
    std::vector<double> cumulative;
    cumulative.reserve(numBins + 1);
    cumulative.push_back(0.0);
    for (std::size_t i = 0; i < numBins; ++i) {
      const double prob = binProbs[i];
      cumulative.push_back(cumulative.back() + prob);
    }
    const double totalProb = cumulative.back();
    const double tol = 1e-6;
    if (std::abs(totalProb - 1.0) > tol) {
      throw std::runtime_error("sampleFromBinned: Sum of bin probabilities must be 1 within tolerance, but is " + std::to_string(totalProb));
    }

    // Sample probability then locate bin
    const double target = randomGen.Uniform(0.0, totalProb);
    std::size_t binIndex = 0;
    for (std::size_t i = 0; i < numBins; ++i) {
      if (target < cumulative[i + 1]) { binIndex = i; break; }
    }

    // Uniform within the selected bin
    const double a = binEdges[binIndex];
    const double b = binEdges[binIndex + 1];
    return a + randomGen.Uniform(0.0, 1.0) * (b - a);
  }

  // modify particle kinematics
  inline void ModifyParticle(int modPdg, std::string modVar, 
                const std::vector<double>& varBinEdges, const std::vector<double>& varBinProbs, int seed
                simb::MCTruth& originalMCTruth, simb::MCTruth& newMCTruth) {
    for (int i = 0; i < originalMCTruth.NParticles(); ++i) {
      TRandom3 randomGen;
      randomGen.SetSeed(seed);

      int pdgCode = originalMCTruth.GetParticle(i).PdgCode();
      int statusCode = originalMCTruth.GetParticle(i).StatusCode();
      if (pdgCode == modPdg && statusCode == 1) { // GENIE status code 1 = kIstStableFinalState
	      std::cout << "Modifying " << modVar << " of particle with pdg code " << modPdg << " at index " << i << std::endl;
        const simb::MCParticle& part = originalMCTruth.GetParticle(i);

        TLorentzVector part_position(part.Vx(), part.Vy(), part.Vz(), part.T());
        TLorentzVector part_momentum = part.Momentum();

        // GENIE vertices relative to nucleus
        //Gvtx(part.Gvx(), part.Gvy(), part.Gvz(), part.Gvt());

        TVector3 part_boost_vector = part_momentum.BoostVector();

        // Calculate decay momenta
        double part_mass = part.Mass();
        //double theta = randomGen.Uniform(0, 2 * M_PI);
        //double phi = acos(1 - 2 * randomGen.Uniform(0, 1));
        double pE = part_momentum.E(); //part_mass / 2;
        double px = part_momentum.Px(); //pE * sin(phi) * cos(theta);
        double py = part_momentum.Py(); //pE * sin(phi) * sin(theta);
        double pz = part_momentum.Pz(); //pE * cos(phi);
        //TLorentzVector rest_frame_momentum(px, py, pz, pE);

        double T = part_position.T();
        double x = part_position.Vx();
        double y = part_position.Vy();
        double z = part_position.Vz();

        // transform the momenta to the lab frame using a Lorentz boost
        //TLorentzVector lab_frame_momentum = rest_frame_momentum;
        //lab_frame_momentum.Boost(part_boost_vector);

        auto new_part_mass = part.Mass();
        auto new_part_position = part_position;
        auto new_part_momentum = part_momentum;
        
        auto new_part_var = sampleFromBinned(varBinEdges, varBinProbs, randomGen);

        if (modVar == "mass"){
          new_part_mass = new_part_var;
        } else if (modVar == "E") {
          new_part_momentum = TLorentzVector mom(px, py, pz, new_part_var);
        } else if (modVar == "px") {
          new_part_momentum = TLorentzVector mom(new_part_var, py, pz, pE);
        } else if (modVar == "py") {
          new_part_momentum = TLorentzVector mom(px, new_part_var, pz, pE);
        } else if (modVar == "pz") {
          new_part_momentum = TLorentzVector mom(px, py, new_part_var, pE);
        } else if (modVar == "theta") {
          double p = new_part_momentum.P();
          double phi = new_part_momentum.Phi();
          double new_px = p * sin(new_part_var) * cos(phi);
          double new_py = p * sin(new_part_var) * sin(phi);
          double new_pz = p * cos(new_part_var);
          new_part_momentum = TLorentzVector mom(new_px, new_py, new_pz, pE);
        } else if (modVar == "phi") {
          double p = new_part_momentum.P();
          double theta = new_part_momentum.Theta();
          double new_px = p * sin(theta) * cos(new_part_var);
          double new_py = p * sin(theta) * sin(new_part_var);
          double new_pz = p * cos(theta);
          new_part_momentum = TLorentzVector mom(new_px, new_py, new_pz, pE);
        } else if (modVar == "p") {
          double theta = new_part_momentum.Theta();
          double phi = new_part_momentum.Phi();
          double new_px = new_part_var * sin(theta) * cos(phi);
          double new_py = new_part_var * sin(theta) * sin(phi);
          double new_pz = new_part_var * cos(theta);
          new_part_momentum = TLorentzVector mom(new_px, new_py, new_pz, pE);
        } else if ( modVar == "x") {
          new_part_position = TLorentzVector(new_part_var, y, z, T);
        } else if ( modVar == "y") {
          new_part_position = TLorentzVector(x, new_part_var, z, T);
        } else if ( modVar == "z") {
          new_part_position = TLorentzVector(x, y, new_part_var, T);
        } else if ( modVar == "T") {
          new_part_position = TLorentzVector(x, y, z, new_part_var);
        } else {
          throw std::runtime_error("ModifyParticle: Invalid modVar " + modVar);
        }

        // make new Part particle with modifications
	      simb::MCParticle newPart(part.TrackId(), part.PdgCode(), part.Process(), part.Mother(), new_part_mass, part.StatusCode());
        newPart.AddTrajectoryPoint(new_part_position, new_part_momentum);

        newMCTruth.Add(newPart);
      } else {
        // Copy over other particles
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
  inline void ModifyParticle(simb::MCTruth& mcTruth) {
    simb::MCTruth newTruth;
    ModifyParticle(mcTruth, newTruth);
    mcTruth = newTruth;
  }

} // namespace evgen

#endif // EVGEN_MODIFYPARTICLE_H
