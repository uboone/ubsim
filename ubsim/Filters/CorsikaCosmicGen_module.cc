////////////////////////////////////////////////////////////////////////
// Class:       CorsikaCosmicGen
// Plugin Type: producer (art v3_01_02)
// File:        CorsikaCosmicGen_module.cc
//
// Generated at Fri Feb  5 10:36:29 2021 by David Caratelli using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TDatabasePDG.h"
#include <memory>

// art extensions
#include "nurandom/RandomUtils/NuRandomService.h"

// larsoft includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/EventGeneratorBase/evgenbase.h"
#include "larcorealg/Geometry/geo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcoreobj/SummaryData/RunData.h"


// CRY related include files
#include "CRYSetup.h"
#include "CRYParticle.h"
#include "CRYGenerator.h"
#include "CLHEP/Random/RandEngine.h"
#include "nutools/EventGeneratorBase/CRY/CRYHelper.h"


class CorsikaCosmicGen;


class CorsikaCosmicGen : public art::EDProducer {
public:
  explicit CorsikaCosmicGen(fhicl::ParameterSet const& p);
  virtual ~CorsikaCosmicGen();
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CorsikaCosmicGen(CorsikaCosmicGen const&) = delete;
  CorsikaCosmicGen(CorsikaCosmicGen&&) = delete;
  CorsikaCosmicGen& operator=(CorsikaCosmicGen const&) = delete;
  CorsikaCosmicGen& operator=(CorsikaCosmicGen&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;
  void beginRun(art::Run& run) override;

private:

  // Declare member data here.

  CRYSetup*      fSetup;
  CRYGenerator*  fGen;
  std::vector<double> fTopLayerDims; //<Dimensions of each layer (linearized multi-dim array layers vs dims (x1,x2,y1,y2,z1,z2))
  std::vector<double> fBottomLayerDims; //<Dimensions of each layer (linearized multi-dim array layers vs dims (x1,x2,y1,y2,z1,z2))
  //double fSampleTime;//<t0 for the particle (at the top box position)
  double fEnergyThresholdLow; //<Only keep particles above this energy
  double fEnergyThresholdHigh; //<Only keep particles below this energy
  double fTimeLow, fTimeHigh; // time interval for simulation [in seconds]
  std::string fCRYConfigStr; //<Configuration string that gets passed to CRY
  int fParticlesPerEvent; //<How many particles to keep per event
  CLHEP::HepRandomEngine& fEngine;

  TRandom3 rand;

};


CorsikaCosmicGen::CorsikaCosmicGen(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fTopLayerDims{p.get< std::vector<double> >( "TopLayerDims" )},
  fBottomLayerDims{p.get< std::vector<double> >( "BottomLayerDims" )},
  //fSampleTime{p.get<double>("SampleTime",0)},
  fEnergyThresholdLow{p.get<double>("EnergyThresholdLow",0.1)},
  fEnergyThresholdHigh{p.get<double>("EnergyThresholdHigh",2.0)},
  fTimeLow{p.get<double>("TimeLow",-0.001)}, // in seconds
  fTimeHigh{p.get<double>("TimeHigh",0.002)}, // in seconds
  fCRYConfigStr{p.get<std::string>("CRYConfigStr")},
  fParticlesPerEvent{p.get<int>("ParticlesPerEvent",1)},
  // create a default random engine; obtain the random seed from NuRandomService,
  // unless overridden in configuration with key "Seed"
  fEngine(art::ServiceHandle<rndm::NuRandomService>()->registerAndSeedEngine(createEngine(0), p, "Seed"))

  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  produces< std::vector<simb::MCTruth> >();
  produces< sumdata::RunData, art::InRun >();

  // Find the pointer to the CRY data tables
  std::string crydatadir;
  const char* datapath = getenv("CRYDATAPATH");
  if( datapath != 0) crydatadir = datapath;
  else{
    mf::LogError("CRYHelper") << "no variable CRYDATAPATH set for cry data location, bail";
    exit(0);
  }

  // Construct the event generator object
  fSetup = new CRYSetup(fCRYConfigStr, crydatadir);
  evgb::RNGWrapper<CLHEP::HepRandomEngine>::set(&fEngine, &CLHEP::HepRandomEngine::flat);
  fSetup->setRandomFunction(evgb::RNGWrapper<CLHEP::HepRandomEngine>::rng);
  fGen = new CRYGenerator(fSetup);


}

CorsikaCosmicGen::~CorsikaCosmicGen(){
  delete fGen;
  delete fSetup;
}

void CorsikaCosmicGen::beginRun(art::Run& run)
{
  // grab the geometry object to see what geometry we are using
  art::ServiceHandle<geo::Geometry> geo;
  run.put(std::make_unique<sumdata::RunData>(geo->DetectorName()), art::fullRun());
}



void CorsikaCosmicGen::produce(art::Event& e)
{
  // Implementation of required member function here.

  std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

  simb::MCTruth mctruth;

   // Generator time at start of sample
  int idctr = 0;

  while (idctr<fParticlesPerEvent) {
    std::vector<CRYParticle*> parts;
    fGen->genEvent(&parts);
    for (unsigned int i=0; i<parts.size(); ++i) {
      // Take ownership of the particle from the vector
      std::unique_ptr<CRYParticle> cryp(parts[i]);

      // Get the energies of the particles
      double ke = cryp->ke()*1.0E-3; // MeV to GeV conversion
      if (ke<fEnergyThresholdLow) continue;
      if (ke>fEnergyThresholdHigh) continue;

      double m = 0.; // in GeV
      static TDatabasePDG*  pdgt = TDatabasePDG::Instance();
      TParticlePDG* pdgp = pdgt->GetParticle(cryp->PDGid());
      if (pdgp) m = pdgp->Mass();

      double etot = ke + m;
      double ptot = etot*etot-m*m;
      if (ptot>0.0) ptot = sqrt(ptot);
      else          ptot = 0.0;

      // Sort out the momentum components. Remember that the NOvA
      // frame has y up and z along the beam. So uvw -> zxy
      double px = ptot * cryp->v();
      double py = ptot * cryp->w();
      double pz = ptot * cryp->u();

      // Particle start position. CRY distributes uniformly in x-y
      // plane at fixed z, where z is the vertical direction. This
      // requires some offsets and rotations to put the particles at
      // the surface in the geometry as well as some rotations
      // since the coordinate frame has y up and z along the
      // beam.
      double vx = cryp->y()*100.0 + 0.5*(fTopLayerDims[1]+fTopLayerDims[0]);
      double vy = fTopLayerDims[3];
      double vz = cryp->x()*100.0 + 0.5*(fTopLayerDims[5]+fTopLayerDims[4]);
      double t  = rand.Uniform(fTimeLow,fTimeHigh); // seconds

      //project to bottom boxes y position and keep if new position is in bottom box
      if(py==0.) continue; //ignore horizontal particles
      double dt=(fBottomLayerDims[3]-vy)/py;
      double nv[]={vx + dt*px, vy + dt*py, vz + dt*pz};

      if( nv[0]<fBottomLayerDims[0] || nv[0]>fBottomLayerDims[1]
          || nv[2]<fBottomLayerDims[4] || nv[2]>fBottomLayerDims[5])
          continue;

      std::cout << "CorsikaCosmicGen...creating new particle with position : "
		<< vx << ", " << vy << ", " << vz 
		<< " and time : " << t*1e9 
		<< std::endl;

      // Boiler plate...
      int istatus    =  1;
      int imother1   = evgb::kCosmicRayGenerator;
      std::string primary("primary");
      simb::MCParticle p(idctr,cryp->PDGid(),primary,imother1,m,istatus);
      TLorentzVector pos(vx,vy,vz,t*1e9);// time needs to be in ns to match GENIE, etc
      TLorentzVector mom(px,py,pz,etot);
      p.AddTrajectoryPoint(pos,mom);

      mctruth.Add(p);
      ++idctr;
      break;
    } // Loop on particles in event
  } // Loop on events simulated

  mctruth.SetOrigin(simb::kCosmicRay);

  truthcol->push_back(mctruth);
  e.put(std::move(truthcol));

  return;


}

DEFINE_ART_MODULE(CorsikaCosmicGen)
