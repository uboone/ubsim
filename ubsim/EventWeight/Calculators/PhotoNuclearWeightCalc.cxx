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
#include "TRandom3.h"
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
  TRandom3 *r3=new TRandom3();

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
  fXSUncertainty = pset.get<float>("xs_uncertainty", 1.0);
  std::string probFileName = pset.get<std::string>("ProbFileName", "photonuprob.root");//systematics/reint/interaction_probabilities.root");//what is this file?
  fNsims = pset.get<int> ("number_of_multisims", 0);

  // Prepare random generator
  fGaussRandom = new CLHEP::RandGaussQ(engine);
  // Load interaction probabilities
  cet::search_path sp("FW_SEARCH_PATH");
  std::string probFilePath = sp.find_file(probFileName);
  fProbFile = TFile::Open(probFilePath.c_str());
  assert(fProbFile && fProbFile->IsOpen());

  // Build parameter list
  for (size_t i=0; i<pars.size(); i++) {
    if (pars[i] == "gamma") {
      fParticles[22] = ParticleDef("gamma",   "prob_nuc_e",  22, sigmas[i], fProbFile);
    }
    else {
      std::cerr << "Unknown particle type: " << pars[i] << std::endl;
      assert(false);
    }
  };

  // Set up universes
  for (auto& it : fParticles) {
    if (mode == "pm1sigma") {
      // pm1sigma mode: 0 = +1sigma, 1 = -1sigma
      it.second.sigmas.push_back( 1.0);
      it.second.sigmas.push_back(-1.0);
      fNsims = 2;
    }
    else if (mode == "multisim") {
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

  // Get MCParticles for each MCTruth in this event
  art::Handle<std::vector<simb::MCTruth> > truthHandle;
  e.getByLabel(fMCTruthProducer, truthHandle);
  const art::FindManyP<simb::MCParticle> truthParticles(truthHandle, e, fMCParticleProducer);
  assert(truthParticles.isValid());

  std::cout << "Event Number : " << e.id() << std::endl;

  // Initialize the vector of event weights
  std::vector<std::vector<double> > weight(truthHandle->size());

  // Loop over sets of MCTruth-associated particles
  for (size_t itruth=0; itruth<truthParticles.size(); itruth++) {

    // Initialize weight vector for this MCTruth
    weight[itruth].resize(fNsims, 1.0);

    // Loop over MCParticles in the event
    auto const& mcparticles = truthParticles.at(itruth);

    int photon_count = 0;
    int primary_photon_count = 0;
    std::vector<int> pi0_trackIDs;

    // collecting pi0s in mcparticles
    for (size_t i=0; i<mcparticles.size(); i++) {
      const simb::MCParticle& p = *mcparticles.at(i);
      int pdg = p.PdgCode();
      if (pdg ==111){
	std::cout <<"pi0 :: EndProcess : " << p.EndProcess() << " , Statuscode : " << p.StatusCode() << " , TrackID : " << p.TrackId() << " , P : " << p.P() << std::endl;
	pi0_trackIDs.push_back(p.TrackId());
      }
    }
    std::cout << "size of pi0 set : " << pi0_trackIDs.size() << std::endl;

    for (size_t i=0; i<mcparticles.size(); i++) {

      const simb::MCParticle& p = *mcparticles.at(i);
      int pdg = p.PdgCode();

      bool is_from_pi0=false;

      if (pdg == 22){
	//	std::cout <<itruth<<"-th truth, "<< photon_count<<"-th photon :: EndProcess : " << p.EndProcess() << " , Mother: " << p.Mother() << " , Statuscode : " << p.StatusCode() << " , TrackID : " << p.TrackId() << " , P : " << p.P() << std::endl;
	photon_count++;

	if (std::find(pi0_trackIDs.begin(), pi0_trackIDs.end(),p.Mother())!=pi0_trackIDs.end()){
	  is_from_pi0 = true;
	  primary_photon_count++;
	  std::cout <<itruth<<"-th truth, "<< primary_photon_count<<"-th primary photon :: EndProcess : " << p.EndProcess() << " , Mother: " << p.Mother() << " , Statuscode : " << p.StatusCode() << " , TrackID : " << p.TrackId() << " , P : " << p.P() << " , primary? :" <<is_from_pi0<<std::endl;
	}
      }
      
      if (!is_from_pi0) continue;

      double ke = p.E() - p.Mass();
      std::string endProc = p.EndProcess();

      bool photonu_absorbed = (endProc.find("photonNuclear") != std::string::npos);
      //      bool interaced = (endProc.find("Inelastic") != std::string::npos);

      // Reweight particles under consideration
      if (fParticles.find(pdg) != fParticles.end()) {
        ParticleDef& def = fParticles[pdg];
        int kebin = def.pint->FindBin(ke);

        // Loop through universes
        for (size_t j=0; j<weight[0].size(); j++) {
          // Integrate a modified cross section to find a survival probability
          float sprob = 1.0;

	  /*
          for (int k=1; k<kebin+1; k++) {
            float wbin = 1.0 + fXSUncertainty * def.sigmas[j];
            float xs = wbin * def.xs->GetBinContent(k);
            sprob *= exp(-1.0 * xs);
          }
	  */

	  //just use prob?
	  float wbin = 1.0 + fXSUncertainty * def.sigmas[j];                                                                             
	  if (wbin<0.) {
	    wbin = 0.;
	    std::cout << "setting wbin to be non-negative, wbin is set to :"<< wbin << std::endl; 
	  }                   
	  //float xs = wbin * def.xs->GetBinContent(kebin);     
	  sprob = 1- wbin*def.pint->GetBinContent(kebin); //exp(-1.0 * xs); 

          float w;
          if (photonu_absorbed) {
            w = (1.0 - sprob) / def.pint->GetBinContent(kebin);
	    std::cout << "photonu end process! , nominal sprob : " << 1-def.pint->GetBinContent(kebin) << " , sprob= "<<sprob<<" , when r= "<< def.sigmas[j] << " , weight  : "<< w <<std::endl; 
          }
          else {
            w = sprob / (1.0 - def.pint->GetBinContent(kebin));
	    std::cout << "other end process !, nominal sprob : " << 1-def.pint->GetBinContent(kebin) << " , sprob= "<<sprob<<" , when r= "<< def.sigmas[j]<< " , weight: " << w << std::endl;
          }
	  

          // Total weight is the product of track weights in the event
	  // Maybe we should just do it  for primary photons
          weight[itruth][j] *= std::max((float)0.0, w);

	  //std::cout << "itruth : " << itruth << " , j-th univ : " << j << " , weight : " << std::max((float)0.0, w) << std::endl;

        }
      }
    }

    std::cout << "photon count : " << photon_count << " , primary photon count : " << primary_photon_count<<std::endl;
    
  }
  for(size_t idx=0; idx<weight.size(); idx++){
    for (size_t jdx=0; jdx<weight[0].size(); jdx++){
      std::cout << idx << "-th truth , " <<jdx <<"-th univ , "<<jdx<<" evt_weight : " << weight[idx][jdx] << std::endl;
    }
  } 
  return weight;
}

REGISTER_WEIGHTCALC(PhotoNuclearWeightCalc)

}  // namespace evwgh
