#ifndef EVGEN_GENERATE_ISOTROPIC_SINGLE_PHOTON_H
#define EVGEN_GENERATE_ISOTROPIC_SINGLE_PHOTON_H

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

// ROOT
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TMath.h"

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

  // Replace the MCTruth particle list with a single photon using
  // binned distributions for energy and polar angle (cosTheta).
  inline void GenerateIsotropicSinglePhotonWithBins(
						    const std::vector<double>& energyBinEdges,
						    const std::vector<double>& energyBinProbs,
						    const std::vector<double>& cosThetaBinEdges,
						    const std::vector<double>& cosThetaBinProbs,
						    simb::MCTruth& originalMCTruth,
						    simb::MCTruth& newMCTruth) {
    TRandom3 randomGen;
    randomGen.SetSeed(0);

    const double E = sampleFromBinned(energyBinEdges, energyBinProbs, randomGen);

    const double cosTheta = sampleFromBinned(cosThetaBinEdges, cosThetaBinProbs, randomGen);
    const double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
    const double phi = randomGen.Uniform(0.0, TMath::TwoPi());
    const double px = E * sinTheta * std::cos(phi);
    const double py = E * sinTheta * std::sin(phi);
    const double pz = E * cosTheta;
    TLorentzVector momentum(px, py, pz, E);

    const simb::MCNeutrino neutrino = originalMCTruth.GetNeutrino();
    const simb::MCParticle nu = neutrino.Nu();
    TLorentzVector position(nu.Vx(), nu.Vy(), nu.Vz(), nu.T());

    simb::MCParticle gamma(1, 22, "primary", 0, 0.0, 1);
    gamma.AddTrajectoryPoint(position, momentum);
    gamma.SetGvtx(nu.Gvx(), nu.Gvy(), nu.Gvz(), nu.Gvt());

    newMCTruth.Add(gamma);

    // Preserve neutrino/meta information
    if (originalMCTruth.NeutrinoSet()) {
      simb::MCNeutrino neutrino = originalMCTruth.GetNeutrino();
      newMCTruth.SetNeutrino(
			     neutrino.CCNC(), neutrino.Mode(), neutrino.InteractionType(),
			     neutrino.Target(), neutrino.HitNuc(), neutrino.HitQuark(),
			     neutrino.W(), neutrino.X(), neutrino.Y(), neutrino.QSqr()
			     );
    }
    newMCTruth.SetOrigin(originalMCTruth.Origin());
  }

  // In-place wrapper requiring explicit bins
  inline void GenerateIsotropicSinglePhotonWithBins(
						    const std::vector<double>& energyBinEdges,
						    const std::vector<double>& energybinProbs,
						    const std::vector<double>& cosThetaBinEdges,
						    const std::vector<double>& cosThetabinProbs,
						    simb::MCTruth& mcTruth) {
    simb::MCTruth newTruth;
    GenerateIsotropicSinglePhotonWithBins(
					  energyBinEdges, energybinProbs, cosThetaBinEdges, cosThetabinProbs,
					  mcTruth, newTruth
					  );
    mcTruth = newTruth;
  }

} // namespace evgen

#endif // EVGEN_GENERATE_ISOTROPIC_SINGLE_PHOTON_H
