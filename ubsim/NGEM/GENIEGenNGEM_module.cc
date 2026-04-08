////////////////////////////////////////////////////////////////////////
//
//
// GENIE neutrino event generator
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef EVGEN_GENIEGEN_H
#define EVGEN_GENIEGEN_H

#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>
#include <memory>
#include <unistd.h>

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TDatabasePDG.h"
#include "TSystem.h"
#include "TStopwatch.h"

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Assns.h"
#include "art/Framework/Core/EDProducer.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// LArSoft includes
#include "lardataalg/MCDumpers/MCDumpers.h" // sim::dump namespace
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "nutools/EventGeneratorBase/GENIE/GENIEHelper.h"
#include "lardata/Utilities/AssociationUtil.h"

// dk2nu extensions
#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/NuChoice.h"
#include "dk2nu/genie/GDk2NuFlux.h"

#include "GENIE/Framework/EventGen/EventRecord.h"
#include "nutools/EventGeneratorBase/GENIE/EVGBAssociationUtil.h"
#include "nutools/EventGeneratorBase/evgenbase.h"

// Optional Modifications, if you don't want to modify normal GENIE neutrino generation events
#include "ManuallyDecayPi0sToTwoPhotons.h"
#include "DeleteOneRandomPhoton.h"
#include "GenerateIsotropicSinglePhoton.h"
#include "ModifyParticle.h"

///Event Generation using GENIE, cosmics or single particles
namespace evgen {
  /**
   * @brief A module to check the results from the Monte Carlo generator
   *
   * Note on random number generator
   * --------------------------------
   * 
   * GENIE uses a TRandom generator for its purposes.
   * Since art's RandomNumberGenerator service only provides
   * `CLHEP::HepRandomEngine`, the standard LArSoft/art mechanism for handling
   * the random stream can't be used.
   * GENIEHelper, interface to GENIE provided by nutools, creates a TRandom
   * that GENIE can use. It initializes it with a random seed read from
   * *RandomSeed* configuration parameter. This and all the other parameters
   * are inherited from the art module (that is, `GENIEGen`) configuration.
   * LArSoft meddles with this mechanism to provide support for the standard
   * "Seed" parameter and NuRandomService service.
   * 
   * Configuration parameters
   * -------------------------
   * 
   * - *RandomSeed* (integer, optional): if specified, this value is used as
   *   seed for GENIE random number generator engine
   * - *Seed* (unsigned integer, optional): if specified, this value is used as
   *   seed for GENIE random number generator engine; if *RandomSeed* is also
   *   specified, this value is ignored (but in the future this could turn into
   *   a configuration error)
   * 
   * As custom, if the random seed is not provided by the configuration, one is
   * fetched from `NuRandomService` (if available), with the behaviour in
     * lar::util::FetchRandomSeed().
   */
  class GENIEGenNGEM : public art::EDProducer {
  public:
    explicit GENIEGenNGEM(fhicl::ParameterSet const &pset);
    virtual ~GENIEGenNGEM();                        

    void produce(art::Event& evt);  
    void beginJob();
    void beginRun(art::Run& run);
    void beginSubRun(art::SubRun& sr);
    void endSubRun(art::SubRun& sr);

  private:

    std::string ParticleStatus(int StatusCode);
    std::string ReactionChannel(int ccnc,int mode);
    
    void FillHistograms(simb::MCTruth mc);

    evgb::GENIEHelper  *fGENIEHelp;       ///< GENIEHelper object
    bool fDefinedVtxHistRange;///use defined hist range; it is useful to have for asymmetric ranges like in DP FD.
    std::vector< double > fVtxPosHistRange; 

    int                 fPassEmptySpills; ///< whether or not to kill evnets with no interactions
    TStopwatch          fStopwatch;       ///keep track of how long it takes to run the job

    double fGlobalTimeOffset;             /// The start of a simulated "beam gate".
    double fRandomTimeOffset;             /// The width of a simulated "beam gate".
    ::sim::BeamType_t fBeamType;          /// The type of beam

    double fPrevTotPOT;      ///< Total POT from subruns previous to current subrun
    double fPrevTotGoodPOT;  ///< Total good POT from subruns previous to current subrun

    TH1F* fGenerated[6];  ///< Spectra as generated

    TH1F* fVertexX;    ///< vertex location of generated events in x
    TH1F* fVertexY;    ///< vertex location of generated events in y
    TH1F* fVertexZ;    ///< vertex location of generated events in z

    TH2F* fVertexXY;   ///< vertex location in xy
    TH2F* fVertexXZ;   ///< vertex location in xz
    TH2F* fVertexYZ;   ///< vertex location in yz

    TH1F* fDCosX;      ///< direction cosine in x
    TH1F* fDCosY;      ///< direction cosine in y
    TH1F* fDCosZ;      ///< direction cosine in z

    TH1F* fMuMomentum; ///< momentum of outgoing muons
    TH1F* fMuDCosX;    ///< direction cosine of outgoing mu in x
    TH1F* fMuDCosY;    ///< direction cosine of outgoing mu in y
    TH1F* fMuDCosZ;    ///< direction cosine of outgoing mu in z

    TH1F* fEMomentum;  ///< momentum of outgoing electrons
    TH1F* fEDCosX;     ///< direction cosine of outgoing e in x
    TH1F* fEDCosY;     ///< direction cosine of outgoing e in y
    TH1F* fEDCosZ;     ///< direction cosine of outgoing e in z

    TH1F* fCCMode;      ///< CC interaction mode
    TH1F* fNCMode;      ///< CC interaction mode

    TH1F* fDeltaE;     ///< difference in neutrino energy from MCTruth::Enu() vs TParticle
    TH1F* fECons;      ///< histogram to determine if energy is conserved in the event

    bool ManuallyDecayPi0s;  ///< whether to manually decay pi0s to photons
    bool DeleteRandomGamma;   ///< whether to delete a random gamma
    bool GenerateIsotropicSinglePhoton; ///< replace all particles with single isotropic gamma
    std::vector<double> SinglePhotonEnergyBinEdges;   ///< configurable via FHiCL
    std::vector<double> SinglePhotonEnergyBinProbs; ///< configurable via FHiCL
    std::vector<double> SinglePhotonCosThetaBinEdges; ///< configurable via FHiCL
    std::vector<double> SinglePhotonCosThetaBinProbs; ///< configurable via FHiCL

    bool fModifyParticle; ///< whether to modify particles according to configuration parameters
    int fModifyParticlePdg; ///< pdg code of particle to modify
    std::vector<std::string> fModifyParticleVar; ///< variables to modify 
    //(options: "E", "p", "mass", "x", "y", "z", "T", "theta", "phi"); configurable via FHiCL
    std::vector<std::vector<double>> fModifyParticleVarBinEdges; ///< configurable via FHiCL
    std::vector<std::vector<double>> fModifyParticleVarBinProbs; ///< configurable via FHiCL
    int fModifyParticleRandomSeed; ///< random seed for particle modification; configurable via FHiCL
    bool fManuallyDecayModifiedPi0s; ///< whether to manually decay pi0s after modification. must use if mass is modified.

  };
}

namespace evgen{

  //____________________________________________________________________________
  GENIEGenNGEM::GENIEGenNGEM(fhicl::ParameterSet const& pset)
    : fGENIEHelp(0)
    , fDefinedVtxHistRange (pset.get< bool >("DefinedVtxHistRange"))
    , fVtxPosHistRange (pset.get< std::vector<double> >("VtxPosHistRange"))
    , fPassEmptySpills (pset.get< bool   >("PassEmptySpills"))
    , fGlobalTimeOffset(pset.get< double >("GlobalTimeOffset",0))
    , fRandomTimeOffset(pset.get< double >("RandomTimeOffset",1600.)) // BNB default value
    , fBeamType(::sim::kBNB)
    , ManuallyDecayPi0s(pset.get<bool>("ManuallyDecayPi0s", false))
    , DeleteRandomGamma(pset.get<bool>("DeleteRandomGamma", false))
    , GenerateIsotropicSinglePhoton(pset.get<bool>("GenerateIsotropicSinglePhoton", false))
    , SinglePhotonEnergyBinEdges(pset.get<std::vector<double>>("SinglePhotonEnergyBinEdges", std::vector<double>()))
    , SinglePhotonEnergyBinProbs(pset.get<std::vector<double>>("SinglePhotonEnergyBinProbs", std::vector<double>()))
    , SinglePhotonCosThetaBinEdges(pset.get<std::vector<double>>("SinglePhotonCosThetaBinEdges", std::vector<double>()))
    , SinglePhotonCosThetaBinProbs(pset.get<std::vector<double>>("SinglePhotonCosThetaBinProbs", std::vector<double>()))
    , fModifyParticle(pset.get<bool>("ModifyParticle", false))
    , fModifyParticlePdg(pset.get<int>("ModifyParticlePdg", 0))
    , fModifyParticleVar(pset.get<std::vector<std::string>>("ModifyParticleVar", std::vector<std::string>()))
    , fModifyParticleVarBinEdges(pset.get<std::vector<std::vector<double>>>("ModifyParticleVarBinEdges", std::vector<std::vector<double>>()))
    , fModifyParticleVarBinProbs(pset.get<std::vector<std::vector<double>>>("ModifyParticleVarBinProbs", std::vector<std::vector<double>>()))
    , fModifyParticleRandomSeed(pset.get<int>("ModifyParticleRandomSeed", 0))
    , fManuallyDecayModifiedPi0s(pset.get<bool>("ManuallyDecayModifiedPi0s", false))
  {
    fStopwatch.Start();

    produces< std::vector<simb::MCTruth> >();
    produces< std::vector<simb::MCFlux>  >();
    produces< std::vector<simb::GTruth>  >();
    produces< sumdata::RunData, art::InRun >();
    produces< sumdata::POTSummary, art::InSubRun >();
    produces< art::Assns<simb::MCTruth, simb::MCFlux> >();
    produces< art::Assns<simb::MCTruth, simb::GTruth> >();
    produces< std::vector<sim::BeamGateInfo> >();

    // dk2nu additions
    produces< std::vector<bsim::Dk2Nu>     >();
    produces< std::vector<bsim::NuChoice>  >();
    produces< art::Assns<simb::MCTruth, bsim::Dk2Nu>    >();
    produces< art::Assns<simb::MCTruth, bsim::NuChoice> >();

    std::string beam_type_name = pset.get<std::string>("BeamName");

    if(beam_type_name == "numi") 

      fBeamType = ::sim::kNuMI;

    else if(beam_type_name == "booster") 

      fBeamType = ::sim::kBNB;

    else

      fBeamType = ::sim::kUnknown;

    art::ServiceHandle<geo::Geometry> geo;

    signed int temp_seed; // the seed read by GENIEHelper is a signed integer...
    fhicl::ParameterSet GENIEconfig(pset);
    if (!GENIEconfig.get_if_present("RandomSeed", temp_seed)) { // TODO use has_key() when it becomes available
      // no RandomSeed specified; check for the LArSoft-style "Seed" instead:
      // obtain the random seed from a service,
      // unless overridden in configuration with key "Seed"
      unsigned int seed;
      if (!GENIEconfig.get_if_present("Seed", seed))
        seed = art::ServiceHandle<rndm::NuRandomService>()->getSeed();
      
      // The seed is not passed to RandomNumberGenerator,
      // since GENIE uses a TRandom generator that is owned by the GENIEHelper.
      // Instead, we explicitly configure the random seed for GENIEHelper:
      GENIEconfig.put("RandomSeed", seed);
    } // if no RandomSeed present
    
    fGENIEHelp = new evgb::GENIEHelper(GENIEconfig, 
				       geo->ROOTGeoManager(),
				       geo->ROOTFile(),
				       geo->TotalMass(pset.get< std::string>("TopVolume").c_str()));
    
  }

  //____________________________________________________________________________
  GENIEGenNGEM::~GENIEGenNGEM()
  {  
    if(fGENIEHelp) delete fGENIEHelp;
    fStopwatch.Stop();
    mf::LogInfo("GENIEProductionTime") << "real time to produce file: " << fStopwatch.RealTime();
  }

  //____________________________________________________________________________
  void GENIEGenNGEM::beginJob(){
    fGENIEHelp->Initialize();

    fPrevTotPOT = 0.;
    fPrevTotGoodPOT = 0.;

    // Get access to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;

    fGenerated[0] = tfs->make<TH1F>("fGenerated_necc","",  100, 0.0, 20.0);
    fGenerated[1] = tfs->make<TH1F>("fGenerated_nebcc","", 100, 0.0, 20.0);
    fGenerated[2] = tfs->make<TH1F>("fGenerated_nmcc","",  100, 0.0, 20.0);
    fGenerated[3] = tfs->make<TH1F>("fGenerated_nmbcc","", 100, 0.0, 20.0);
    fGenerated[4] = tfs->make<TH1F>("fGenerated_nnc","",   100, 0.0, 20.0);
    fGenerated[5] = tfs->make<TH1F>("fGenerated_nbnc","",  100, 0.0, 20.0);
    
    fDCosX = tfs->make<TH1F>("fDCosX", ";dx/ds", 200, -1., 1.);
    fDCosY = tfs->make<TH1F>("fDCosY", ";dy/ds", 200, -1., 1.);
    fDCosZ = tfs->make<TH1F>("fDCosZ", ";dz/ds", 200, -1., 1.);

    fMuMomentum = tfs->make<TH1F>("fMuMomentum", ";p_{#mu} (GeV/c)", 500, 0., 50.);
    fMuDCosX    = tfs->make<TH1F>("fMuDCosX", ";dx/ds;", 200, -1., 1.);
    fMuDCosY    = tfs->make<TH1F>("fMuDCosY", ";dy/ds;", 200, -1., 1.);
    fMuDCosZ    = tfs->make<TH1F>("fMuDCosZ", ";dz/ds;", 200, -1., 1.);

    fEMomentum  = tfs->make<TH1F>("fEMomentum", ";p_{e} (GeV/c)", 500, 0., 50.);
    fEDCosX     = tfs->make<TH1F>("fEDCosX", ";dx/ds;", 200, -1., 1.);
    fEDCosY     = tfs->make<TH1F>("fEDCosY", ";dy/ds;", 200, -1., 1.);
    fEDCosZ     = tfs->make<TH1F>("fEDCosZ", ";dz/ds;", 200, -1., 1.);

    fCCMode = tfs->make<TH1F>("fCCMode", ";CC Interaction Mode;", 4, 0., 4.);
    fCCMode->GetXaxis()->SetBinLabel(1, "QE");
    fCCMode->GetXaxis()->SetBinLabel(2, "Res");
    fCCMode->GetXaxis()->SetBinLabel(3, "DIS");
    fCCMode->GetXaxis()->SetBinLabel(4, "Coh");
    fCCMode->GetXaxis()->CenterLabels();

    fNCMode = tfs->make<TH1F>("fNCMode", ";NC Interaction Mode;", 4, 0., 4.);
    fNCMode->GetXaxis()->SetBinLabel(1, "QE");
    fNCMode->GetXaxis()->SetBinLabel(2, "Res");
    fNCMode->GetXaxis()->SetBinLabel(3, "DIS");
    fNCMode->GetXaxis()->SetBinLabel(4, "Coh");
    fNCMode->GetXaxis()->CenterLabels();

    fDeltaE = tfs->make<TH1F>("fDeltaE", ";#Delta E_{#nu} (GeV);", 200, -1., 1.); 
    fECons  = tfs->make<TH1F>("fECons", ";#Delta E(#nu,lepton);", 500, -5., 5.);

    art::ServiceHandle<geo::Geometry> geo;
    double x = 2.1*geo->DetHalfWidth();
    double y = 2.1*geo->DetHalfHeight();
    double z = 2.*geo->DetLength();
    int xdiv = TMath::Nint(2*x/5.);
    int ydiv = TMath::Nint(2*y/5.);
    int zdiv = TMath::Nint(2*z/5.);

    if (fDefinedVtxHistRange == false)
      {
	fVertexX = tfs->make<TH1F>("fVertexX", ";x (cm)", xdiv, -0.1*x, x);
	fVertexY = tfs->make<TH1F>("fVertexY", ";y (cm)", ydiv, -y,     y);
	fVertexZ = tfs->make<TH1F>("fVertexZ", ";z (cm)", zdiv, -0.1*z, z);
    
	fVertexXY = tfs->make<TH2F>("fVertexXY", ";x (cm);y (cm)", xdiv, -0.1*x, x, ydiv,     -y, y);
	fVertexXZ = tfs->make<TH2F>("fVertexXZ", ";z (cm);x (cm)", zdiv, -0.2*z, z, xdiv, -0.1*x, x);
	fVertexYZ = tfs->make<TH2F>("fVertexYZ", ";z (cm);y (cm)", zdiv, -0.2*z, z, ydiv,     -y, y);
      }
    else
      {
        fVertexX = tfs->make<TH1F>("fVertexX", ";x (cm)", xdiv, fVtxPosHistRange[0], fVtxPosHistRange[1]);
	fVertexY = tfs->make<TH1F>("fVertexY", ";y (cm)", ydiv, fVtxPosHistRange[2], fVtxPosHistRange[3]);
	fVertexZ = tfs->make<TH1F>("fVertexZ", ";z (cm)", zdiv, fVtxPosHistRange[4], fVtxPosHistRange[5]);
    
	fVertexXY = tfs->make<TH2F>("fVertexXY", ";x (cm);y (cm)", xdiv, fVtxPosHistRange[0], fVtxPosHistRange[1], ydiv,     fVtxPosHistRange[2], fVtxPosHistRange[3]);
	fVertexXZ = tfs->make<TH2F>("fVertexXZ", ";z (cm);x (cm)", zdiv, fVtxPosHistRange[4], fVtxPosHistRange[5], xdiv,     fVtxPosHistRange[0], fVtxPosHistRange[1]);
	fVertexYZ = tfs->make<TH2F>("fVertexYZ", ";z (cm);y (cm)", zdiv, fVtxPosHistRange[4], fVtxPosHistRange[5], ydiv,     fVtxPosHistRange[2], fVtxPosHistRange[3]);
      }

  }

  //____________________________________________________________________________
  void GENIEGenNGEM::beginRun(art::Run& run)
  {

    // grab the geometry object to see what geometry we are using
    art::ServiceHandle<geo::Geometry> geo;
    std::unique_ptr<sumdata::RunData> runcol(new sumdata::RunData(geo->DetectorName()));

    run.put(std::move(runcol));

    return;
  }

  //____________________________________________________________________________
  void GENIEGenNGEM::beginSubRun(art::SubRun& sr)
  {

    fPrevTotPOT = fGENIEHelp->TotalExposure();
    fPrevTotGoodPOT = fGENIEHelp->TotalExposure();

    return;
  }

  //____________________________________________________________________________
  void GENIEGenNGEM::endSubRun(art::SubRun& sr)
  {

    std::unique_ptr<sumdata::POTSummary> p(new sumdata::POTSummary());
    
    p->totpot = fGENIEHelp->TotalExposure() - fPrevTotPOT;
    p->totgoodpot = fGENIEHelp->TotalExposure() - fPrevTotGoodPOT;

    sr.put(std::move(p));

    return;
  }

  //____________________________________________________________________________
  void GENIEGenNGEM::produce(art::Event& evt)
  {
    std::unique_ptr< std::vector<simb::MCTruth> > truthcol  (new std::vector<simb::MCTruth>);
    std::unique_ptr< std::vector<simb::MCFlux>  > fluxcol   (new std::vector<simb::MCFlux >);
    std::unique_ptr< std::vector<simb::GTruth>  > gtruthcol (new std::vector<simb::GTruth >);
    std::unique_ptr< art::Assns<simb::MCTruth, simb::MCFlux> > tfassn(new art::Assns<simb::MCTruth, simb::MCFlux>);
    std::unique_ptr< art::Assns<simb::MCTruth, simb::GTruth> > tgtassn(new art::Assns<simb::MCTruth, simb::GTruth>);
    std::unique_ptr< std::vector<sim::BeamGateInfo> > gateCollection(new std::vector<sim::BeamGateInfo>);

    std::unique_ptr< std::vector<bsim::Dk2Nu> > 
      dk2nucol(new std::vector<bsim::Dk2Nu>);
    std::unique_ptr< std::vector<bsim::NuChoice> > 
      nuchoicecol(new std::vector<bsim::NuChoice>);    
    std::unique_ptr< art::Assns<simb::MCTruth, bsim::Dk2Nu> > 
      dk2nuassn(new art::Assns<simb::MCTruth, bsim::Dk2Nu>);
    std::unique_ptr< art::Assns<simb::MCTruth, bsim::NuChoice> > 
      nuchoiceassn(new art::Assns<simb::MCTruth, bsim::NuChoice>);

    while(truthcol->size() < 1){
      while(!fGENIEHelp->Stop()){
	
	simb::MCTruth truth;
	simb::MCFlux  flux;
	simb::GTruth  gTruth;

	// GENIEHelper returns a false in the sample method if 
	// either no neutrino was generated, or the interaction
	// occurred beyond the detector's z extent - ie something we
	// would never see anyway.
	if(fGENIEHelp->Sample(truth, flux, gTruth)){

	  std::cout << "truth.NParticles(): " << truth.NParticles() << std::endl;
	  for (int i = 0; i < truth.NParticles(); ++i) {
	    std::cout << "    pdg: " << truth.GetParticle(i).PdgCode();
	    std::cout << ", track_id: " << truth.GetParticle(i).TrackId();
	    std::cout << ", mother: " << truth.GetParticle(i).Mother();
	    std::cout << ", mass: " << truth.GetParticle(i).Mass();
	    std::cout << ", position: (" << truth.GetParticle(i).Vx() << ", " << truth.GetParticle(i).Vy() << ", " << truth.GetParticle(i).Vz() << ")";
	    std::cout << ", momentum: (" << truth.GetParticle(i).Px() << ", " << truth.GetParticle(i).Py() << ", " << truth.GetParticle(i).Pz() << ")";
	    std::cout << ", Gvtx: (" << truth.GetParticle(i).Gvx() << ", " << truth.GetParticle(i).Gvy() << ", " << truth.GetParticle(i).Gvz() << ", " << truth.GetParticle(i).Gvt() << ")";
	    std::cout << ", status code: " << truth.GetParticle(i).StatusCode();
	    std::cout << ", process: " << truth.GetParticle(i).Process();
	    std::cout << std::endl;
	  }

	  if (GenerateIsotropicSinglePhoton) {
	    std::cout << "Generating isotropic single photon" << std::endl;
	    GenerateIsotropicSinglePhotonWithBins(
						  SinglePhotonEnergyBinEdges, SinglePhotonEnergyBinProbs,
						  SinglePhotonCosThetaBinEdges, SinglePhotonCosThetaBinProbs,
        truth
						  );
	    std::cout << "truth.NParticles(): " << truth.NParticles() << std::endl;
	    for (int i = 0; i < truth.NParticles(); ++i) {
	      std::cout << "    pdg: " << truth.GetParticle(i).PdgCode();
	      std::cout << ", track_id: " << truth.GetParticle(i).TrackId();
	      std::cout << ", mother: " << truth.GetParticle(i).Mother();
	      std::cout << ", mass: " << truth.GetParticle(i).Mass();
	      std::cout << ", position: (" << truth.GetParticle(i).Vx() << ", " << truth.GetParticle(i).Vy() << ", " << truth.GetParticle(i).Vz() << ")";
	      std::cout << ", momentum: (" << truth.GetParticle(i).Px() << ", " << truth.GetParticle(i).Py() << ", " << truth.GetParticle(i).Pz() << ")";
	      std::cout << ", Gvtx: (" << truth.GetParticle(i).Gvx() << ", " << truth.GetParticle(i).Gvy() << ", " << truth.GetParticle(i).Gvz() << ", " << truth.GetParticle(i).Gvt() << ")";
	      std::cout << ", status code: " << truth.GetParticle(i).StatusCode();
	      std::cout << ", process: " << truth.GetParticle(i).Process();
	      std::cout << std::endl;
	    }
	  }

	  if (ManuallyDecayPi0s) {
	    std::cout << "Manually decaying pi0s" << std::endl;
	    ManuallyDecayPi0sToTwoPhotons(truth);
	    std::cout << "truth.NParticles(): " << truth.NParticles() << std::endl;
	    for (int i = 0; i < truth.NParticles(); ++i) {
	      std::cout << "    pdg: " << truth.GetParticle(i).PdgCode();
	      std::cout << ", track_id: " << truth.GetParticle(i).TrackId();
	      std::cout << ", mother: " << truth.GetParticle(i).Mother();
	      std::cout << ", mass: " << truth.GetParticle(i).Mass();
	      std::cout << ", position: (" << truth.GetParticle(i).Vx() << ", " << truth.GetParticle(i).Vy() << ", " << truth.GetParticle(i).Vz() << ")";
	      std::cout << ", momentum: (" << truth.GetParticle(i).Px() << ", " << truth.GetParticle(i).Py() << ", " << truth.GetParticle(i).Pz() << ")";
	      std::cout << ", Gvtx: (" << truth.GetParticle(i).Gvx() << ", " << truth.GetParticle(i).Gvy() << ", " << truth.GetParticle(i).Gvz() << ", " << truth.GetParticle(i).Gvt() << ")";
	      std::cout << ", status code: " << truth.GetParticle(i).StatusCode();
	      std::cout << ", process: " << truth.GetParticle(i).Process();
	      std::cout << std::endl;
	    }
	  }

	  if (DeleteRandomGamma) {
	    std::cout << "Deleting random gamma" << std::endl;
	    DeleteOneRandomPhoton(truth);
	    std::cout << "truth.NParticles(): " << truth.NParticles() << std::endl;
	    for (int i = 0; i < truth.NParticles(); ++i) {
	      std::cout << "    pdg: " << truth.GetParticle(i).PdgCode();
	      std::cout << ", track_id: " << truth.GetParticle(i).TrackId();
	      std::cout << ", mother: " << truth.GetParticle(i).Mother();
	      std::cout << ", mass: " << truth.GetParticle(i).Mass();
	      std::cout << ", position: (" << truth.GetParticle(i).Vx() << ", " << truth.GetParticle(i).Vy() << ", " << truth.GetParticle(i).Vz() << ")";
	      std::cout << ", momentum: (" << truth.GetParticle(i).Px() << ", " << truth.GetParticle(i).Py() << ", " << truth.GetParticle(i).Pz() << ")";
	      std::cout << ", Gvtx: (" << truth.GetParticle(i).Gvx() << ", " << truth.GetParticle(i).Gvy() << ", " << truth.GetParticle(i).Gvz() << ", " << truth.GetParticle(i).Gvt() << ")";
	      std::cout << ", status code: " << truth.GetParticle(i).StatusCode();
	      std::cout << ", process: " << truth.GetParticle(i).Process();
	      std::cout << std::endl;
	    }
	  }

    if (fModifyParticle) {
      std::cout << "Modifying particles with pdg code " << fModifyParticlePdg << std::endl;
      if (fModifyParticleVar.size() != fModifyParticleVarBinEdges.size() || fModifyParticleVar.size() != fModifyParticleVarBinProbs.size()) {
        std::cerr << "Error: fModifyParticleVar, fModifyParticleVarBinEdges, and fModifyParticleVarBinProbs must have the same size" << std::endl;
        return;
      }
      for (uint i_v = 0; i_v < fModifyParticleVar.size(); ++i_v) {
        std::cout << "    Modifying variable " << fModifyParticleVar.at(i_v) << std::endl;
        ModifyParticle(fModifyParticlePdg, fModifyParticleVar.at(i_v), fModifyParticleVarBinEdges.at(i_v), fModifyParticleVarBinProbs.at(i_v), fModifyParticleRandomSeed, truth);
        if (fManuallyDecayModifiedPi0s){ //fModifyParticleVar.at(i_v) == "mass" && fModifyParticlePdg == 111) {
          std::cout << "Manually decaying pi0s" << std::endl;
          ManuallyDecayPi0sToTwoPhotons(truth);
	      }
      }
      std::cout << "truth.NParticles(): " << truth.NParticles() << std::endl;
      for (int i = 0; i < truth.NParticles(); ++i) {
        std::cout << "    pdg: " << truth.GetParticle(i).PdgCode();
        std::cout << ", track_id: " << truth.GetParticle(i).TrackId();
        std::cout << ", mother: " << truth.GetParticle(i).Mother();
        std::cout << ", mass: " << truth.GetParticle(i).Mass();
        std::cout << ", position: (" << truth.GetParticle(i).Vx() << ", " << truth.GetParticle(i).Vy() << ", " << truth.GetParticle(i).Vz() << ")";
        std::cout << ", momentum: (" << truth.GetParticle(i).Px() << ", " << truth.GetParticle(i).Py() << ", " << truth.GetParticle(i).Pz() << ")";
        std::cout << ", Gvtx: (" << truth.GetParticle(i).Gvx() << ", " << truth.GetParticle(i).Gvy() << ", " << truth.GetParticle(i).Gvz() << ", " << truth.GetParticle(i).Gvt() << ")";
        std::cout << ", status code: " << truth.GetParticle(i).StatusCode();
        std::cout << ", process: " << truth.GetParticle(i).Process();
        std::cout << std::endl;
      }
    }

	  truthcol ->push_back(truth);
	  fluxcol  ->push_back(flux);
	  gtruthcol->push_back(gTruth);
	  util::CreateAssn(*this, evt, *truthcol, *fluxcol, *tfassn, fluxcol->size()-1, fluxcol->size());
	  util::CreateAssn(*this, evt, *truthcol, *gtruthcol, *tgtassn, gtruthcol->size()-1, gtruthcol->size());
	    
	  FillHistograms(truth);
	    
	  genie::GFluxI* fdriver = fGENIEHelp->GetFluxDriver(true);
	  genie::flux::GDk2NuFlux* dk2nuDriver = 
	    dynamic_cast<genie::flux::GDk2NuFlux*>(fdriver);
	  if ( dk2nuDriver ) {
	    const bsim::Dk2Nu& dk2nuObj = dk2nuDriver->GetDk2Nu();
	    dk2nucol   ->push_back(dk2nuObj);
	    const bsim::NuChoice& nuchoiceObj = dk2nuDriver->GetNuChoice();
	    nuchoicecol->push_back(nuchoiceObj);
	    util::CreateAssn(*this, evt, *truthcol, *dk2nucol, *dk2nuassn,
			     dk2nucol->size()-1, dk2nucol->size());
	    util::CreateAssn(*this, evt, *truthcol, *nuchoicecol, *nuchoiceassn,
			     nuchoicecol->size()-1, nuchoicecol->size());
	  }

	  // check that the process code is not unsupported by GENIE
	  // (see issue #18025 for reference);
	  // if it is, print all the information we can about this truth record
	  if (truth.NeutrinoSet() && (truth.GetNeutrino().InteractionType() == simb::kNuanceOffset)) {
	    mf::LogWarning log("GENIEmissingProcessMapping");
	    log << "Found an interaction that is not represented by the interaction type code in GENIE:"
	            "\nMCTruth record:"
	            "\n"
	      ;
	    sim::dump::DumpMCTruth(log, truth, 2U); // 2 trajectory points per line
	    log <<
	            "\nGENIE truth record:"
	            "\n"
	      ;
	    sim::dump::DumpGTruth(log, gTruth);
	  } // if 
	    
	}// end if genie was able to make an event

      }// end event generation loop
      
      // check to see if we are to pass empty spills
      if(truthcol->size() < 1 && fPassEmptySpills){
	MF_LOG_DEBUG("GENIEGenNGEM") << "no events made for this spill but "
				 << "passing it on and ending the event anyway";
	break;
      }

    }// end loop while no interactions are made

    // Create a simulated "beam gate" for these neutrino events.
    // We're creating a vector of these because, in a
    // distant-but-possible future, we may be generating more than one
    // beam gate within a simulated time window.
    gateCollection->push_back(sim::BeamGateInfo( fGlobalTimeOffset, fRandomTimeOffset, fBeamType ));

    // put the collections in the event 
    evt.put(std::move(truthcol));
    evt.put(std::move(fluxcol));
    evt.put(std::move(gtruthcol));
    evt.put(std::move(tfassn));
    evt.put(std::move(tgtassn));
    evt.put(std::move(gateCollection));

    evt.put(std::move(dk2nucol));
    evt.put(std::move(nuchoicecol));
    evt.put(std::move(dk2nuassn));
    evt.put(std::move(nuchoiceassn));

    return;
  }

  //......................................................................
  std::string GENIEGenNGEM::ParticleStatus(int StatusCode)
  {
    int code = StatusCode;
    std::string ParticleStatusName;

    switch(code)
      {
      case -1:
	ParticleStatusName = "kIStUndefined";
	break;
      case 0:
	ParticleStatusName = "kIStInitialState";
	break;
      case 1:
	ParticleStatusName = "kIStStableFinalState";
	break;
      case 2:
	ParticleStatusName = "kIStIntermediateState";
	break;
      case 3:
	ParticleStatusName = "kIStDecayedState";
	break;
      case 11:
	ParticleStatusName = "kIStNucleonTarget";
	break;
      case 12:
	ParticleStatusName = "kIStDISPreFragmHadronicState";
	break;
      case 13:
	ParticleStatusName = "kIStPreDecayResonantState";
	break;
      case 14:
	ParticleStatusName = "kIStHadronInTheNucleus";
	break;
      case 15:
	ParticleStatusName = "kIStFinalStateNuclearRemnant";
	break;
      case 16:
	ParticleStatusName = "kIStNucleonClusterTarget";
	break;
      default:
	ParticleStatusName = "Status Unknown";
      }
    return ParticleStatusName;
  }

  //......................................................................
  std::string GENIEGenNGEM::ReactionChannel(int ccnc,int mode)
  {
    std::string ReactionChannelName=" ";

    if(ccnc==0)
      ReactionChannelName = "kCC";
    else if(ccnc==1)
      ReactionChannelName = "kNC";
    else std::cout<<"Current mode unknown!! "<<std::endl;

    if(mode==0)
      ReactionChannelName += "_kQE";
    else if(mode==1)
      ReactionChannelName += "_kRes";
    else if(mode==2)
      ReactionChannelName += "_kDIS";
    else if(mode==3)
      ReactionChannelName += "_kCoh";
    else std::cout<<"interaction mode unknown!! "<<std::endl;

    return ReactionChannelName;
  }

  //......................................................................
  void GENIEGenNGEM::FillHistograms(simb::MCTruth mc)
  {
    // Decide which histograms to put the spectrum in
    int id = -1;
    if (mc.GetNeutrino().CCNC()==simb::kCC) {
      fCCMode->Fill(mc.GetNeutrino().Mode());
      if      (mc.GetNeutrino().Nu().PdgCode() ==  12) id = 0;
      else if (mc.GetNeutrino().Nu().PdgCode() == -12) id = 1;
      else if (mc.GetNeutrino().Nu().PdgCode() ==  14) id = 2;
      else if (mc.GetNeutrino().Nu().PdgCode() == -14) id = 3;
      else return;
    }
    else {
      fNCMode->Fill(mc.GetNeutrino().Mode());
      if (mc.GetNeutrino().Nu().PdgCode() > 0) id = 4;
      else                                     id = 5;
    }
    if (id==-1) abort();
  
    // Fill the specta histograms
    fGenerated[id]->Fill(mc.GetNeutrino().Nu().E() );
      
    ///< fill the vertex histograms from the neutrino - that is always 
    ///< particle 0 in the list
    simb::MCNeutrino       mcnu = mc.GetNeutrino();
    const simb::MCParticle nu   = mcnu.Nu();

    fVertexX->Fill(nu.Vx());
    fVertexY->Fill(nu.Vy());
    fVertexZ->Fill(nu.Vz());

    fVertexXY->Fill(nu.Vx(), nu.Vy());
    fVertexXZ->Fill(nu.Vz(), nu.Vx());
    fVertexYZ->Fill(nu.Vz(), nu.Vy());

    double mom = nu.P();
    if(std::abs(mom) > 0.){
      fDCosX->Fill(nu.Px()/mom);
      fDCosY->Fill(nu.Py()/mom);
      fDCosZ->Fill(nu.Pz()/mom);
    }


    MF_LOG_DEBUG("GENIEInteractionInformation") 
      << std::endl
      << "REACTION:  " << ReactionChannel(mc.GetNeutrino().CCNC(),mc.GetNeutrino().Mode()) 
      << std::endl
      << "-----------> Particles in the Stack = " << mc.NParticles() << std::endl
      << std::setiosflags(std::ios::left) 
      << std::setw(20) << "PARTICLE"
      << std::setiosflags(std::ios::left) 
      << std::setw(32) << "STATUS"
      << std::setw(18) << "E (GeV)"
      << std::setw(18) << "m (GeV/c2)"
      << std::setw(18) << "Ek (GeV)"
      << std::endl << std::endl;

    const TDatabasePDG* databasePDG = TDatabasePDG::Instance();

    // Loop over the particle stack for this event 
    for(int i = 0; i < mc.NParticles(); ++i){
      simb::MCParticle part(mc.GetParticle(i));
      std::string name = databasePDG->GetParticle(part.PdgCode())->GetName();
      int code = part.StatusCode();
      std::string status = ParticleStatus(code);
      double mass = part.Mass();
      double energy = part.E(); 
      double Ek = (energy-mass); // Kinetic Energy (GeV)
      if(status=="kIStStableFinalState"||status=="kIStHadronInTheNucleus")
	MF_LOG_DEBUG("GENIEFinalState") 
	  << std::setiosflags(std::ios::left) << std::setw(20) << name
	  << std::setiosflags(std::ios::left) << std::setw(32) <<status
	  << std::setw(18)<< energy
	  << std::setw(18)<< mass
	  << std::setw(18)<< Ek <<std::endl;
      else 
	MF_LOG_DEBUG("GENIEFinalState")
	  << std::setiosflags(std::ios::left) << std::setw(20) << name
	  << std::setiosflags(std::ios::left) << std::setw(32) << status
	  << std::setw(18) << energy
	  << std::setw(18) << mass <<std::endl; 
    }


    if(mc.GetNeutrino().CCNC() == simb::kCC){
  
      ///look for the outgoing lepton in the particle stack
      ///just interested in the first one
      for(int i = 0; i < mc.NParticles(); ++i){
	simb::MCParticle part(mc.GetParticle(i));
	if(abs(part.PdgCode()) == 11){
	  fEMomentum->Fill(part.P());
	  fEDCosX->Fill(part.Px()/part.P());
	  fEDCosY->Fill(part.Py()/part.P());
	  fEDCosZ->Fill(part.Pz()/part.P());
	  fECons->Fill(nu.E() - part.E());
	  break;
	}
	else if(abs(part.PdgCode()) == 13){
	  fMuMomentum->Fill(part.P());
	  fMuDCosX->Fill(part.Px()/part.P());
	  fMuDCosY->Fill(part.Py()/part.P());
	  fMuDCosZ->Fill(part.Pz()/part.P());
	  fECons->Fill(nu.E() - part.E());
	  break;
	}
      }// end loop over particles
    }//end if CC interaction

    return;
  }

}

namespace evgen{

  DEFINE_ART_MODULE(GENIEGenNGEM)

}

#endif 
