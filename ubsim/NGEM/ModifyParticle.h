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

#include "GenerateIsotropicSinglePhoton.h"

namespace evgen {

  // modify particle kinematics
  inline void ModifyParticle(int modPdg, std::string modVar, 
    const std::vector<double>& varBinEdges, const std::vector<double>& varBinProbs, 
    int seed, simb::MCTruth& originalMCTruth, simb::MCTruth& newMCTruth) {

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
        double pE = part_momentum.E(); 
        //double px = part_momentum.Px(); 
        //double py = part_momentum.Py(); 
        //double pz = part_momentum.Pz(); 
        //TLorentzVector rest_frame_momentum(px, py, pz, pE);

        double T = part.T();
        double x = part.Vx();
        double y = part.Vy();
        double z = part.Vz();

        // transform the momenta to the lab frame using a Lorentz boost
        //TLorentzVector lab_frame_momentum = rest_frame_momentum;
        //lab_frame_momentum.Boost(part_boost_vector);

        auto new_part_mass = part_mass;
        auto new_part_position = part_position;
        auto new_part_momentum = part_momentum;
        
        auto new_part_var = sampleFromBinned(varBinEdges, varBinProbs, randomGen);

        if (modVar == "mass"){
          new_part_mass = new_part_var;
          double new_pE = pE;
          // Compute new |p| from E^2 = p^2 + m^2
          double new_p2 = new_pE*new_pE - new_part_mass*new_part_mass;
          if (new_p2 < 0) new_p2 = 0; // protect against numerical issues / invalid input
          double new_p = std::sqrt(new_p2);

          // Keep original direction
          TVector3 dir = part_momentum.Vect().Unit();

          double new_px = new_p * dir.X();
          double new_py = new_p * dir.Y();
          double new_pz = new_p * dir.Z();
          TLorentzVector mom(new_px, new_py, new_pz, new_pE);
          new_part_momentum = mom;
        } else if (modVar == "E") { 
          double new_pE = new_part_var;
          // Compute new |p| from E^2 = p^2 + m^2
          double new_p2 = new_pE*new_pE - part_mass*part_mass;
          if (new_p2 < 0) new_p2 = 0; // protect against numerical issues / invalid input
          double new_p = std::sqrt(new_p2);

          // Keep original direction
          TVector3 dir = part_momentum.Vect().Unit();

          double new_px = new_p * dir.X();
          double new_py = new_p * dir.Y();
          double new_pz = new_p * dir.Z();
          TLorentzVector mom(new_px, new_py, new_pz, new_pE);
          new_part_momentum = mom;
        } else if (modVar == "theta") {
          double p = new_part_momentum.P();
          double phi = new_part_momentum.Phi();
          double new_px = p * sin(new_part_var) * cos(phi);
          double new_py = p * sin(new_part_var) * sin(phi);
          double new_pz = p * cos(new_part_var);
          TLorentzVector mom(new_px, new_py, new_pz, pE);
          new_part_momentum = mom;
        } else if (modVar == "phi") {
          double p = new_part_momentum.P();
          double theta = new_part_momentum.Theta();
          double new_px = p * sin(theta) * cos(new_part_var);
          double new_py = p * sin(theta) * sin(new_part_var);
          double new_pz = p * cos(theta);
          TLorentzVector mom(new_px, new_py, new_pz, pE);
          new_part_momentum = mom;
        } else if (modVar == "p") {
          double new_p = new_part_var;
          // Compute new |p| from E^2 = p^2 + m^2
          double new_pE = std::sqrt(new_p*new_p + part_mass*part_mass);
          if (new_pE < 0) new_pE = 0; // protect against numerical issues / invalid input

          // Keep original direction
          TVector3 dir = part_momentum.Vect().Unit();

          double new_px = new_p * dir.X();
          double new_py = new_p * dir.Y();
          double new_pz = new_p * dir.Z();
          TLorentzVector mom(new_px, new_py, new_pz, new_pE);
          new_part_momentum = mom;
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
  inline void ModifyParticle(int modPdg, std::string modVar, 
    const std::vector<double>& varBinEdges, const std::vector<double>& varBinProbs, 
    int seed, simb::MCTruth& mcTruth) {
    simb::MCTruth newTruth;
    ModifyParticle(modPdg, modVar, varBinEdges, varBinProbs, seed, mcTruth, newTruth);
    mcTruth = newTruth;
  }

} // namespace evgen

#endif // EVGEN_MODIFYPARTICLE_H
