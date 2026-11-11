#ifndef EVGEN_RADCORR_COLLINEAR_PHOTON_H
#define EVGEN_RADCORR_COLLINEAR_PHOTON_H

// Sampling of an energetic collinear photon radiated from a charged-current
// final-state lepton, using the one-loop collinear function of
//   Tomalak et al., Phys. Rev. D 106, 093006 (2022), Eqs. (35)-(36)
//
// Variables (all energies in GeV):
//   E_tree = E_lep + E_gamma          (tree-level lepton "jet" energy)
//   x      = E_lep / E_tree           (so E_gamma = (1 - x) E_tree)
//   eta    = dtheta * E_tree / m_lep  (dtheta = lepton-photon opening angle [rad])
//
// j(x, eta) = dsigma^gamma / (dx dsigma_LO) for photons inside a cone of
// size eta, so d^2 sigma^gamma / (dx deta dsigma_LO) = dj/deta.

#include <iostream>
#include <cmath>
#include <algorithm>

// ROOT
#include "TRandom3.h"

namespace evgen {
namespace radcorr {

  constexpr double kAlphaQED = 1.0 / 137.036;

  // Eq. (35): j(mu/m_l, x, eta)
  inline double CollinearJ(double x, double eta) {
    if (x <= 0.0 || x >= 1.0 || eta <= 0.0) return 0.0;
    const double xe2 = x * x * eta * eta;
    return (kAlphaQED / M_PI) * (0.5 * (1.0 + x * x) / (1.0 - x) * std::log1p(xe2)
                                 - x / (1.0 - x) * xe2 / (1.0 + xe2));
  }

  // dj/deta, the double-differential density in (x, eta)
  inline double CollinearDJDEta(double x, double eta) {
    if (x <= 0.0 || x >= 1.0 || eta <= 0.0) return 0.0;
    const double xe2 = x * x * eta * eta;
    return (kAlphaQED / M_PI) * (x * x * eta) / ((1.0 - x) * (1.0 + xe2))
      * ((1.0 + x * x) - 2.0 * x / (1.0 + xe2));
  }

  // Allowed region in x (the photon must have E_gamma > deltaE, and the lepton E_lep > m_lep)
  inline double XMin(double eTree, double mLep) { return mLep / eTree; }
  inline double XMax(double eTree, double deltaE) { return 1.0 - deltaE / eTree; }

  // Total probability of radiating a photon in the allowed region:
  //   P = integral of dj/deta over x in [xMin, xMax], eta in [0, etaMax]
  //     = integral of j(x, etaMax) over x in [xMin, xMax]   (since j(x, 0) = 0)
  // computed with a simple midpoint sum.
  inline double RadiationProbability(double eTree, double mLep, double deltaE, double maxAngleRad) {
    const double xMin = XMin(eTree, mLep);
    const double xMax = XMax(eTree, deltaE);
    const double etaMax = maxAngleRad * eTree / mLep;
    if (xMin >= xMax || etaMax <= 0.0) return 0.0;
    const int nSteps = 100000;
    const double dx = (xMax - xMin) / nSteps;
    double prob = 0.0;
    for (int i = 0; i < nSteps; ++i) prob += CollinearJ(xMin + (i + 0.5) * dx, etaMax) * dx;
    return prob;
  }

  // Sample (x, eta) from dj/deta in the allowed region with rejection sampling:
  // draw (x, eta) uniformly in the box, and accept with probability dj/deta / fMax.
  // Returns false if the region is empty.
  inline bool SampleXEta(double eTree, double mLep, double deltaE, double maxAngleRad,
                         TRandom3& randomGen, double& x, double& eta) {
    const double xMin = XMin(eTree, mLep);
    const double xMax = XMax(eTree, deltaE);
    const double etaMax = maxAngleRad * eTree / mLep;
    if (xMin >= xMax || etaMax <= 0.0) return false;

    // maximum of the density, from a scan on a grid (the maximum is near x = xMax,
    // eta ~ 1), with a safety margin
    const int nGrid = 200;
    double fMax = 0.0;
    for (int i = 0; i <= nGrid; ++i) {
      for (int k = 0; k <= nGrid; ++k) {
        fMax = std::max(fMax, CollinearDJDEta(xMin + (xMax - xMin) * i / nGrid, etaMax * k / nGrid));
      }
    }
    fMax *= 1.2;
    if (fMax <= 0.0) return false;

    for (long iTry = 0; iTry < 100000000; ++iTry) {
      x = randomGen.Uniform(xMin, xMax);
      eta = randomGen.Uniform(0.0, etaMax);
      const double f = CollinearDJDEta(x, eta);
      if (f > fMax) {
        std::cout << "Warning: SampleXEta: density " << f << " above the assumed maximum " << fMax << std::endl;
      }
      if (randomGen.Uniform(0.0, fMax) < f) return true;
    }
    std::cout << "Warning: SampleXEta: no (x, eta) accepted" << std::endl;
    return false;
  }

} // namespace radcorr
} // namespace evgen

#endif // EVGEN_RADCORR_COLLINEAR_PHOTON_H
