/**
 * \class evwgh::PhotoNuclearWeightCalc
 * \brief PhotoNuclear interaction event rewighting
 * \author Y. Jwa <yj2429@columbia.edu>, 2019/12
 *
 * \Largily adopting ReinteractionWeightCalc written by A. Mastbaum <mastbaum@uchicago.edu>
 *
 * Reweight events based on photonuclear interaction possibilities.
 */

#include <map>
#include <string>
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "Geant4/G4LossTableManager.hh"
#include "Geant4/G4ParticleTable.hh"
#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4MaterialCutsCouple.hh"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Provenance/ModuleContext.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "larsim/EventWeight/Base/WeightCalcCreator.h"
#include "larsim/EventWeight/Base/WeightCalc.h"

namespace evwgh {

class PhotoNuclearWeightCalc : public WeightCalc {
public:
  PhotoNuclearWeightCalc() {}

  void Configure(fhicl::ParameterSet const& p, CLHEP::HepRandomEngine& engine);

  std::vector<std::vector<double> > GetWeight(art::Event& e);

  /**
   * \struct ParticleDef
   * \brief A reweightable particle definition
   */
  class ParticleDef {
  public:
    ParticleDef() {}
    ParticleDef(std::string _name, std::string objname,
                int _pdg, float _sigma, TFile* probFile)
        : name(_name), pdg(_pdg), par_sigma(_sigma) {
      pint = dynamic_cast<TH1D*>(probFile->Get(objname.c_str()));
      assert(pint);

      std::cout << "objname : " << objname.c_str() << std::endl;

      // Reconstitute the cross section vs. KE
      char hname[100];
      snprintf(hname, 100, "_xs_%s", name.c_str());

      std::cout << "hname : " << hname << std::endl; 

      xs = dynamic_cast<TH1D*>(pint->Clone(hname));
      assert(xs);

      for (int j=1; j<xs->GetNbinsX()+1; j++) {
        float p1 = pint->GetBinContent(j);
        float p2 = pint->GetBinContent(j-1);
        float v = 0;

        if (p1 > p2 && p1 < 1) {
          v = -1.0 * log((1.0 - p1) / (1.0 - p2));
        }

        xs->SetBinContent(j, v);
      }
    }

    std::string name;  //!< String name
    int pdg;  //!< PDG code
    float par_sigma;  //!< Variation sigma set by user
    TH1D* pint;  //!< Interaction probability as a function of KE
    TH1D* xs;  //!< Derived effective cross section
    std::vector<double> sigmas;  //!< Sigmas for universes
  };

private:
  std::string fMCParticleProducer;  //!< Label for the MCParticle producer
  std::string fMCTruthProducer;  //!< Label for the MCTruth producer
  CLHEP::RandGaussQ* fGaussRandom;  //!< Random number generator
  TFile* fProbFile;  //!< File with interaction probabilities, uncertainties
  std::map<int, ParticleDef> fParticles;  //!< Particles to reweight
  unsigned fNsims;  //!< Number of multisims
  float fXSUncertainty;  //!< Flat cross section uncertainty

  DECLARE_WEIGHTCALC(PhotoNuclearWeightCalc)
};


void PhotoNuclearWeightCalc::Configure(fhicl::ParameterSet const& p,
                                        CLHEP::HepRandomEngine& engine)
{
  // Get configuration for this function
  fhicl::ParameterSet const& pset = p.get<fhicl::ParameterSet>(GetName());
  fMCParticleProducer = pset.get<std::string>("MCParticleProducer", "largeant");
  fMCTruthProducer = pset.get<std::string>("MCTruthProducer", "generator");
  std::vector<std::string> pars = pset.get< std::vector<std::string> >("parameter_list");//"gamma"	
  std::vector<float> sigmas = pset.get<std::vector<float> >("parameter_sigma");	
  std::string mode = pset.get<std::string>("mode");//"multisim"
  fXSUncertainty = pset.get<float>("xs_uncertainty", 0.3);
  std::string probFileName = pset.get<std::string>("ProbFileName", "photonuprob.root");//systematics/reint/interaction_probabilities.root");//what is this file?
  fNsims = pset.get<int> ("number_of_multisims", 0);

  std::cout << ":: 1 :: here? " << std::endl;

  // Prepare random generator
  fGaussRandom = new CLHEP::RandGaussQ(engine);

  std::cout << ":: 2 :: here? "<< std::endl;

  // Load interaction probabilities
  cet::search_path sp("FW_SEARCH_PATH");

  std::cout << ":: 3 :: here? "<< std::endl;

  std::string probFilePath = sp.find_file(probFileName);

  std::cout << ":: 4 :: here? "<< std::endl;

  fProbFile = TFile::Open(probFilePath.c_str());

  std::cout << ":: 5 :: here? "<< std::endl;

  assert(fProbFile && fProbFile->IsOpen());

  std::cout << ":: 6 :: here? "<< std::endl;

  // Build parameter list
  for (size_t i=0; i<pars.size(); i++) {

    std::cout << ":: 7 :: here? #pars? " << pars.size() << " , iterator? " << i << " , #sigmas? "  << sigmas.size() << std::endl;

    if (pars[i] == "gamma") {

      std::cout << ":: 8 :: here? "<< std::endl;

      fParticles[22] = ParticleDef("gamma",   "prob_nuc_e",  22, sigmas[i], fProbFile);
    }
    else {

      std::cout << ":: 8.1 :: here? "<< std::endl;

      std::cerr << "Unknown particle type: " << pars[i] << std::endl;
      assert(false);
    }
  };

  std::cout << ":: 9 :: here? "<< std::endl;

  // Set up universes
  for (auto& it : fParticles) {

    std::cout << ":: 10 :: here? "<< std::endl;

    if (mode == "pm1sigma") {

      std::cout << ":: 11 :: here? "<< std::endl;

      // pm1sigma mode: 0 = +1sigma, 1 = -1sigma
      it.second.sigmas.push_back( 1.0);
      it.second.sigmas.push_back(-1.0);
      fNsims = 2;
    }
    else if (mode == "multisim") {

      std::cout << ":: 12 :: here? "<< std::endl;

      // multisim mode: Scale factors sampled within the given uncertainty
      for (unsigned j=0; j<fNsims; j++) {
        double r = fGaussRandom->fire(0.0, 1.0);
        it.second.sigmas.push_back(it.second.par_sigma * r);
      }
    }
    else {
      // Anything else mode: Set scale to user-defined scale factor
      it.second.sigmas.push_back(it.second.par_sigma);
      fNsims = 1;
    }
  }
}

  

std::vector<std::vector<double> >
PhotoNuclearWeightCalc::GetWeight(art::Event& e) {

  std::cout << ":: 13 :: here? "<< std::endl;

  // Get MCParticles for each MCTruth in this event
  art::Handle<std::vector<simb::MCTruth> > truthHandle;
  e.getByLabel(fMCTruthProducer, truthHandle);
  const art::FindManyP<simb::MCParticle> truthParticles(truthHandle, e, fMCParticleProducer);
  assert(truthParticles.isValid());

  // Initialize the vector of event weights
  std::vector<std::vector<double> > weight(truthHandle->size());

  // Loop over sets of MCTruth-associated particles
  for (size_t itruth=0; itruth<truthParticles.size(); itruth++) {

    // Initialize weight vector for this MCTruth
    weight[itruth].resize(fNsims, 1.0);

    // Loop over MCParticles in the event
    auto const& mcparticles = truthParticles.at(itruth);

    for (size_t i=0; i<mcparticles.size(); i++) {
      const simb::MCParticle& p = *mcparticles.at(i);
      int pdg = p.PdgCode();
      double ke = p.E() - p.Mass();
      std::string endProc = p.EndProcess();

      bool photonu_absorbed = (endProc.find("Photonuclear") != std::string::npos);
      //      bool interaced = (endProc.find("Inelastic") != std::string::npos);

      // Reweight particles under consideration
      if (fParticles.find(pdg) != fParticles.end()) {
        ParticleDef& def = fParticles[pdg];
        int kebin = def.pint->FindBin(ke);

        // Loop through universes
        for (size_t j=0; j<weight[0].size(); j++) {
          // Integrate a modified cross section to find a survival probability
          float sprob = 1.0;
          for (int k=1; k<kebin+1; k++) {
            float wbin = 1.0 + fXSUncertainty * def.sigmas[j];
            float xs = wbin * def.xs->GetBinContent(k);
            sprob *= exp(-1.0 * xs);
          }

          float w;
          if (photonu_absorbed) {
            w = (1.0 - sprob) / def.pint->GetBinContent(kebin);
          }
          else {
            w = sprob / (1.0 - def.pint->GetBinContent(kebin));
          }

          // Total weight is the product of track weights in the event
          weight[itruth][j] *= std::max((float)0.0, w);
        }
      }
    }
  }

  return weight;
}

REGISTER_WEIGHTCALC(PhotoNuclearWeightCalc)

}  // namespace evwgh
