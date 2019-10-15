/**
 * \class evwgh::Geant4WeightCalc
 * \brief Updated hadron reinteraction event reweighting, using geant4reweight
 * \author K. Duffy <kduffy@fnal.gov>, 2019/10
 *
 * Reweight events based on hadron reinteraction probabilities.
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
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
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
#include "larcore/Geometry/Geometry.h"
#include "geant4reweight/src/ReweightBase/G4Reweighter.hh"
#include "geant4reweight/src/ReweightBase/G4ReweightTraj.hh"
#include "geant4reweight/src/ReweightBase/G4ReweightStep.hh"
#include "geant4reweight/src/PropBase/G4ReweightParameterMaker.hh"

// local include
#include "BetheBlochForG4ReweightValid.h"

namespace evwgh {

class Geant4WeightCalc : public WeightCalc {
public:
  Geant4WeightCalc() {}

  void Configure(fhicl::ParameterSet const& p, CLHEP::HepRandomEngine& engine);

  std::vector<std::vector<double> > GetWeight(art::Event& e);

  /**
   * \struct ParticleDef
   * \brief A reweightable particle definition
   */
  // class ParticleDef {
  // public:
  //   ParticleDef() {}
  //   ParticleDef(std::string _name, std::string objname,
  //               int _pdg, float _sigma, TFile* probFile)
  //       : name(_name), pdg(_pdg), par_sigma(_sigma) {
  //     pint = dynamic_cast<TH1D*>(probFile->Get(objname.c_str()));
  //     assert(pint);
  //
  //     // Reconstitute the cross section vs. KE
  //     char hname[100];
  //     snprintf(hname, 100, "_xs_%s", name.c_str());
  //     xs = dynamic_cast<TH1D*>(pint->Clone(hname));
  //     assert(xs);
  //
  //     for (int j=1; j<xs->GetNbinsX()+1; j++) {
  //       float p1 = pint->GetBinContent(j);
  //       float p2 = pint->GetBinContent(j-1);
  //       float v = 0;
  //
  //       if (p1 > p2 && p1 < 1) {
  //         v = -1.0 * log((1.0 - p1) / (1.0 - p2));
  //       }
  //
  //       xs->SetBinContent(j, v);
  //     }
  //   }
  //
  //   std::string name;  //!< String name
  //   int pdg;  //!< PDG code
  //   float par_sigma;  //!< Variation sigma set by user
  //   TH1D* pint;  //!< Interaction probability as a function of KE
  //   TH1D* xs;  //!< Derived effective cross section
  //   std::vector<double> sigmas;  //!< Sigmas for universes
  // };

private:
  std::string fMCParticleProducer;  //!< Label for the MCParticle producer
  std::string fMCTruthProducer;  //!< Label for the MCTruth producer
  CLHEP::RandGaussQ* fGaussRandom;  //!< Random number generator
  // std::map<int, ParticleDef> fParticles;  //!< Particles to reweight
  unsigned fNsims;  //!< Number of multisims
  // float fXSUncertainty;  //!< Flat cross section uncertainty
  G4Reweighter *theReweighter;
  G4ReweightParameterMaker *ParMaker;
  std::vector<std::map<std::string, double>> UniverseVals; //!< Vector of maps relating parameter name to value (defines parameter values that will be evaluated in universes). Each map should have one entry per parameter we are considering

  art::ServiceHandle < geo::Geometry > fGeometryService;

  bool fMakeOutputTrees; ///!< Fcl parameter to decide whether to save output tree (useful for validations but not necessary when running in production)
  TTree *fOutTree_MCTruth; //!< Output tree for quick validations: on entry per MCTruth object
  TTree *fOutTree_Particle; //!< Output tree for quick validations: one entry per neutrino-induced pi+, pi-, or proton
  int event_num; //!< Variables for both output trees
  int run_num; //!< Variables for both output trees
  int subrun_num; //!< Variables for both output trees
  double p_track_length; //!< Variables for by-particle output tree
  int p_PDG; //!< Variables for by-particle output tree
  std::string p_final_proc; //!< Variables for by-particle output tree
  double p_init_momentum; //!< Variables for by-particle output tree
  double p_final_momentum; //!< Variables for by-particle output tree
  std::vector< double > p_energies; //!< Variables for by-particle output tree
  std::vector< int > p_sliceInts; //!< Variables for by-particle output tree
  int p_nElasticScatters; //!< Variables for by-particle output tree
  std::vector<double> p_inel_weight; //!< Variables for by-particle output tree
  std::vector<double> p_elastic_weight; //!< Variables for by-particle output
  std::vector<std::vector<double> > e_inel_weight; //!< Variables for by-event output tree
  std::vector<std::vector<double> > e_elastic_weight; //!< Variables for by-event output tree

  DECLARE_WEIGHTCALC(Geant4WeightCalc)
};


void Geant4WeightCalc::Configure(fhicl::ParameterSet const& p,
                                        CLHEP::HepRandomEngine& engine)
{
  std::cout << "Using Geant4WeightCalc for reinteraction weights" << std::endl;

  // Get configuration for this function
  fhicl::ParameterSet const& pset = p.get<fhicl::ParameterSet>(GetName());
  fMCParticleProducer = pset.get<std::string>("MCParticleProducer", "largeant");
  fMCTruthProducer = pset.get<std::string>("MCTruthProducer", "generator");
  fMakeOutputTrees = pset.get< bool >( "makeoutputtree", false );
  std::string mode = pset.get<std::string>("mode");
  std::string FracsFileName = pset.get< std::string >( "fracsfile" );
  std::string XSecFileName = pset.get< std::string >( "xsecfile" );
  std::vector< fhicl::ParameterSet > FitParSets = pset.get< std::vector< fhicl::ParameterSet > >("parameters");
  fNsims = pset.get<int> ("number_of_multisims", 0);

  // Prepare random generator
  fGaussRandom = new CLHEP::RandGaussQ(engine);

  // Get input files
  TFile FracsFile( FracsFileName.c_str(), "OPEN" );
  TFile XSecFile( XSecFileName.c_str(), "OPEN" );

  // Configure G4Reweighter
  ParMaker = new G4ReweightParameterMaker( FitParSets );

  theReweighter = new G4Reweighter( &FracsFile, ParMaker->GetFSHists(), ParMaker->GetElasticHist() );
  theReweighter->SetTotalGraph(&XSecFile);

  // Make output trees to save things for quick and easy validation
  art::ServiceHandle<art::TFileService> tfs;

  if (fMakeOutputTrees){
    fOutTree_Particle = tfs->make<TTree>("ByParticleValidTree","");
    fOutTree_Particle->Branch("event_num",&event_num);
    fOutTree_Particle->Branch("run_num",&run_num);
    fOutTree_Particle->Branch("subrun_num",&subrun_num);
    fOutTree_Particle->Branch("weight",&p_inel_weight);
    fOutTree_Particle->Branch("elastic_weight",&p_elastic_weight);
    fOutTree_Particle->Branch("track_length",&p_track_length);
    fOutTree_Particle->Branch("PDG",&p_PDG);
    fOutTree_Particle->Branch("final_proc",&p_final_proc);
    fOutTree_Particle->Branch("init_momentum",&p_init_momentum);
    fOutTree_Particle->Branch("final_momentum",&p_final_momentum);
    fOutTree_Particle->Branch("energies",&p_energies);
    fOutTree_Particle->Branch("sliceInts",&p_sliceInts);
    fOutTree_Particle->Branch("nElasticScatters",&p_nElasticScatters);

    fOutTree_MCTruth = tfs->make<TTree>("ByMCTruthValidTree","");
    fOutTree_MCTruth->Branch("event_num",&event_num);
    fOutTree_MCTruth->Branch("run_num",&run_num);
    fOutTree_MCTruth->Branch("subrun_num",&subrun_num);
    fOutTree_MCTruth->Branch("inelastic_weight",&e_inel_weight);
    fOutTree_MCTruth->Branch("elastic_weight",&e_elastic_weight);

  }


  // Read input parameter sets and set up universes
  size_t n_parsets = FitParSets.size();
  std::vector<std::string> FitParNames;
  std::vector<double> FitParNominals;
  std::vector<double> FitParSigmas;
  std::map<std::string, double> theNominals;

  for (size_t i_parset=0; i_parset<n_parsets; ++i_parset){
    fhicl::ParameterSet theSet = FitParSets.at(i_parset);
    std::string theName = theSet.get<std::string>("Name");
    double theNominal = theSet.get<double>("Nominal",1.);
    double theSigma = theSet.get<double>("Sigma",0.);

    FitParNames.push_back(theName);
    FitParNominals.push_back(theNominal);
    FitParSigmas.push_back(theSigma);

    theNominals[theName] = theNominal;
  }

  if (mode=="pm1sigma"){
    // pm1sigma mode: 0 = +1sigma, 1 = -1sigma of a single parameter. All other parameters at nominal
    for (size_t i_parset=0; i_parset<n_parsets; ++i_parset){
      // For each parameter, first create a nominal parameter set (one for +1sigma and one for -1sigma)
      std::map<std::string, double> tmp_vals_p1sigma(theNominals);
      std::map<std::string, double> tmp_vals_m1sigma(theNominals);
      // Now reset the +1sigma and -1sigma values for this parameter set only
      tmp_vals_p1sigma[FitParNames.at(i_parset)] = FitParNominals.at(i_parset)+FitParSigmas.at(i_parset);
      tmp_vals_m1sigma[FitParNames.at(i_parset)] = FitParNominals.at(i_parset)-FitParSigmas.at(i_parset);
      // Finally, add these universes into the vector
      UniverseVals.push_back(tmp_vals_p1sigma);
      UniverseVals.push_back(tmp_vals_m1sigma);
    } // end loop over parsets (i)
    fNsims = 2;
  } // pm1sigma
  else if (mode=="multisim"){
    // multisim mode: parameter values sample within the given uncertainty for all parameters simultaneously
    // Loop over universes j
    for (unsigned j=0; j<fNsims; j++){
      // In each multisim universe, loop through all parameters. For each parameter, generate a new random number from Nominal-Sigma to Nominal+Sigma.
      std::map<std::string, double> tmp_vals;
      for (size_t i_parset=0; i_parset<n_parsets; ++i_parset){
        double r = fGaussRandom->fire(-1.0,1.0);
        tmp_vals[FitParNames.at(i_parset)] = FitParNominals.at(i_parset)+(FitParSigmas.at(i_parset)*r);
      } // loop over parameters (i_parset)
      // Now save this universe
      UniverseVals.push_back(tmp_vals);
    } // loop over Nsims (j)
  } // multisim
  else{
    // Anything else mode: Set parameters to user-defined nominal value
    UniverseVals.push_back(theNominals);
    fNsims = 1;
  } // any other mode
}


std::vector<std::vector<double> >
Geant4WeightCalc::GetWeight(art::Event& e) {

  // Get event/run/subrun numbers for output
  run_num = e.run();
  subrun_num = e.subRun();
  event_num = e.id().event();

  // Get MCParticles for each MCTruth in this event
  art::Handle<std::vector<simb::MCTruth> > truthHandle;
  e.getByLabel(fMCTruthProducer, truthHandle);
  const art::FindManyP<simb::MCParticle> truthParticles(truthHandle, e, fMCParticleProducer);
  assert(truthParticles.isValid());

  // Initialize the vector of event weights
  std::vector<std::vector<double> > weight(truthHandle->size());
  // These two are just for saving to the output tree for fast validation
  e_inel_weight.resize(truthHandle->size());
  e_elastic_weight.resize(truthHandle->size());

  // Loop over sets of MCTruth-associated particles
  for (size_t itruth=0; itruth<truthParticles.size(); itruth++) {

    // Initialize weight vector for this MCTruth
    weight[itruth].resize(fNsims, 1.0);
    // These ones are just for saving to the output tree for fast validation
    e_inel_weight[itruth].resize(fNsims,1.0);
    e_elastic_weight[itruth].resize(fNsims,1.0);

    // Loop over MCParticles in the event
    auto const& mcparticles = truthParticles.at(itruth);

    for (size_t i=0; i<mcparticles.size(); i++) {
      // Initialize weight vectors for saving to the output tree for fast validation
      p_inel_weight.resize(fNsims,1.0);
      p_elastic_weight.resize(fNsims,1.0);

      const simb::MCParticle& p = *mcparticles.at(i);
      int pdg = p.PdgCode();
      int mcpID = p.TrackId();
      std::string EndProcess  = p.EndProcess();

      double mass = 0.;
      if( pdg == 211 ) mass = 139.57;
      else if( pdg == 2212 ) mass = 938.28;

      // We only want to record weights for true charged pions and protons, so skip other particles
      // Actually for now, consider pi+ only
      // TODO add pi- and protons!
      //if ( ( TMath::Abs(pdg) == 211 || pdg == 2212 ) ){
      if (pdg == 211){
        // Get GEANT trajectory points: weighting will depend on position and momentum at each trajectory point so calculate those
        std::vector<double> trajpoint_X;
        std::vector<double> trajpoint_Y;
        std::vector<double> trajpoint_Z;
        std::vector<double> trajpoint_PX;
        std::vector<double> trajpoint_PY;
        std::vector<double> trajpoint_PZ;
        std::vector<int> elastic_indices;

        //Get the list of processes from the true trajectory
        const std::vector< std::pair< size_t, unsigned char > > & processes = p.Trajectory().TrajectoryProcesses();
        std::map< size_t, std::string > process_map;
        for( auto it = processes.begin(); it != processes.end(); ++it ){
          process_map[ it->first ] = p.Trajectory().KeyToProcess( it->second );
        }

        for( size_t i = 0; i < p.NumberTrajectoryPoints(); ++i ){
          double X = p.Position(i).X();
          double Y = p.Position(i).Y();
          double Z = p.Position(i).Z();
          geo::Point_t testpoint1 { X, Y, Z };
          const TGeoMaterial* testmaterial1 = fGeometryService->Material( testpoint1 );
          //For now, just going to reweight the points within the LAr of the TPC
          // TODO check if this is right
          if ( !strcmp( testmaterial1->GetName(), "LAr" ) ){
            trajpoint_X.push_back( X );
            trajpoint_Y.push_back( Y );
            trajpoint_Z.push_back( Z );

            trajpoint_PX.push_back( p.Px(i) );
            trajpoint_PY.push_back( p.Py(i) );
            trajpoint_PZ.push_back( p.Pz(i) );

            auto itProc = process_map.find(i);
            if( itProc != process_map.end() && itProc->second == "hadElastic" ){
              //Push back the index relative to the start of the reweightable steps
              elastic_indices.push_back( trajpoint_X.size() - 1 );
            }

          }

        } // end loop over trajectory points

        // Now find daughters of the MCP
        std::vector<int> daughter_PDGs;
        std::vector<int> daughter_IDs;
        for( int i_mcp = 0; i_mcp < p.NumberDaughters(); i_mcp++ ){
          int daughterID = p.Daughter(i_mcp);
          for (auto test_mcp : mcparticles){
            if (test_mcp->TrackId() == daughterID){
              int pid = test_mcp->PdgCode();
              daughter_PDGs.push_back(pid);
              daughter_IDs.push_back( test_mcp->TrackId() );
              break;
            }
          }
        } // end loop over daughters

        // --- Now that we have all the information about the track we need, here comes the reweighting part! --- //

        //Make a G4ReweightTraj -- This is the reweightable object
        G4ReweightTraj theTraj(mcpID, pdg, 0, event_num, std::make_pair(0,0));

        //Create its set of G4ReweightSteps and add them to the Traj (note: this needs to be done once per MCParticle but will be valid for all weight calculations)
        std::vector< G4ReweightStep* > allSteps;

        size_t nSteps = trajpoint_PX.size();
        if( nSteps < 2 ) continue;

        p_nElasticScatters = elastic_indices.size();

        for( size_t istep = 1; istep < nSteps; ++istep ){

          // if( istep == trajpoint_PX.size() - 1 && std::find( elastic_indices.begin(), elastic_indices.end(), j ) != elastic_indices.end() )
          //   std::cout << "Warning: last step an elastic process" << std::endl;

          std::string proc = "default";
          if( istep == trajpoint_PX.size() - 1 )
            proc = EndProcess;
          else if( std::find( elastic_indices.begin(), elastic_indices.end(), istep ) != elastic_indices.end() )
            proc = "hadElastic";

          double deltaX = ( trajpoint_X.at(istep) - trajpoint_X.at(istep-1) );
          double deltaY = ( trajpoint_Y.at(istep) - trajpoint_Y.at(istep-1) );
          double deltaZ = ( trajpoint_Z.at(istep) - trajpoint_Z.at(istep-1) );

          double len = sqrt(
            std::pow( deltaX, 2 )  +
            std::pow( deltaY, 2 )  +
            std::pow( deltaZ, 2 )
          );

          double preStepP[3] = {
            trajpoint_PX.at(istep-1)*1.e3,
            trajpoint_PY.at(istep-1)*1.e3,
            trajpoint_PZ.at(istep-1)*1.e3
          };

          double postStepP[3] = {
            trajpoint_PX.at(istep)*1.e3,
            trajpoint_PY.at(istep)*1.e3,
            trajpoint_PZ.at(istep)*1.e3
          };

          if( istep == 1 ){
            theTraj.SetEnergy( sqrt( preStepP[0]*preStepP[0] + preStepP[1]*preStepP[1] + preStepP[2]*preStepP[2] + mass*mass ) );
          }

          G4ReweightStep * theStep = new G4ReweightStep( mcpID, pdg, 0, event_num, preStepP, postStepP, len, proc );
          theStep->SetDeltaX( deltaX );
          theStep->SetDeltaY( deltaY );
          theStep->SetDeltaZ( deltaZ );

          theTraj.AddStep( theStep );

          for( size_t k = 0; k < daughter_PDGs.size(); ++k ){
            theTraj.AddChild( new G4ReweightTraj(daughter_IDs[k], daughter_PDGs[k], mcpID, event_num, std::make_pair(0,0) ) );
          }
        } // end loop over nSteps (istep)
        p_track_length = theTraj.GetTotalLength();

        p_init_momentum = sqrt( theTraj.GetEnergy()*theTraj.GetEnergy() - mass*mass );
        p_final_momentum = sqrt(
            std::pow( theTraj.GetStep( theTraj.GetNSteps() - 1 )->GetPreStepPx(), 2 ) +
            std::pow( theTraj.GetStep( theTraj.GetNSteps() - 1 )->GetPreStepPy(), 2 ) +
            std::pow( theTraj.GetStep( theTraj.GetNSteps() - 1 )->GetPreStepPz(), 2 )
        );

        std::vector< std::pair< double, int > > thin_slice = ThinSliceBetheBloch( &theTraj, .5 );

        p_energies.clear();
        p_sliceInts.clear();
        for( size_t islice = 0; islice < thin_slice.size(); ++islice ){
          p_energies.push_back( thin_slice[islice].first );
          p_sliceInts.push_back( thin_slice[islice].second );
        }

        // Loop through universes (j)
        for (size_t j=0; j<weight[0].size(); j++) {
          float w, el_w;

          // I think this is the only bit that needs to change for different universes -- all the above is jut about the track, which doesn't change based on universe
          ParMaker->SetNewVals(UniverseVals.at(j));
          theReweighter->SetNewHists(ParMaker->GetFSHists());

          //Get the weight from the G4ReweightTraj
          w = theReweighter->GetWeight( &theTraj );
          // Total weight is the product of track weights in the event
          weight[itruth][j] *= std::max((float)0.0, w);

          // Do the same for elastic weight (should be 1 unless set to non-nominal )
          el_w = theReweighter->GetElasticWeight( &theTraj );
          weight[itruth][j] *= std::max((float)0.0,el_w);

          // just for the output tree
          p_inel_weight[j] = w;
          p_elastic_weight[j] = el_w;
          e_inel_weight[itruth][j] *= std::max((float)0.0,w);
          e_elastic_weight[itruth][j] *= std::max((float)0.0,el_w);

        } // loop through universes (j)
      } // if ( ( TMath::Abs(pdg) == 211 || pdg == 2212 ) )
      if (fMakeOutputTrees) fOutTree_Particle->Fill();
    } // loop over mcparticles (i)
    if (fMakeOutputTrees) fOutTree_MCTruth->Fill();
  } // loop over sets of MCtruth-associated particles (itruth)

return weight;

}

REGISTER_WEIGHTCALC(Geant4WeightCalc)

}  // namespace evwgh
