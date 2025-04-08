////////////////////////////////////////////////////////////////////////
// Class:       UBPhotonLibraryPropagation
// Plugin Type: producer (art v2_05_00)
// File:        UBPhotonLibraryPropagation_module.cc
//
// Generated at Tue Mar 21 07:45:42 2017 by Wesley Ketchum using cetskelgen
// from cetlib version v1_21_00.
// Modified by Vincent Basque January 2025 to include the transport time of the photons calling the propagationtime class from legacylarg4
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
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "nurandom/RandomUtils/NuRandomService.h"

#include <memory>
#include <iostream>

#include "larcore/Geometry/Geometry.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/Simulation/PhotonVoxels.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"

#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larsim/IonizationScintillation/ISCalcSeparate.h"

#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larsim/PhotonPropagation/PropagationTimeModel.h" //ns timing new!


namespace phot {
  class UBPhotonLibraryPropagation;
}


class phot::UBPhotonLibraryPropagation : public art::EDProducer {
public:
  explicit UBPhotonLibraryPropagation(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  UBPhotonLibraryPropagation(UBPhotonLibraryPropagation const &) = delete;
  UBPhotonLibraryPropagation(UBPhotonLibraryPropagation &&) = delete;
  UBPhotonLibraryPropagation & operator = (UBPhotonLibraryPropagation const &) = delete;
  UBPhotonLibraryPropagation & operator = (UBPhotonLibraryPropagation &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  double fRiseTimeFast;
  double fRiseTimeSlow;
  bool   fDoSlowComponent;
  bool   fIncludePhotPropTimeUBSim;
  bool   fUsingScaleFactor;

  std::vector<art::InputTag> fEDepTags;
  std::vector<double> fPhotonScale;

  larg4::ISCalcSeparate fISAlg;

  CLHEP::HepRandomEngine& fPhotonEngine;
  CLHEP::HepRandomEngine& fScintEngine;

  double GetScintYield(sim::SimEnergyDeposit const&, detinfo::LArProperties const&);
  // Old Method
  //double GetScintTime(double scint_time, double rise_time, double, double);
  // New Method
  double GetScintTime(double rise_time, double scint_time, CLHEP::RandFlat& randflatscinttime);
  double bi_exp(double t, double tau1, double tau2);
  double single_exp(double t, double tau2);
  // propagation time model
  std::unique_ptr<PropagationTimeModel> fPropTimeModel;

};


phot::UBPhotonLibraryPropagation::UBPhotonLibraryPropagation(fhicl::ParameterSet const & p) :
  art::EDProducer(p),
  fRiseTimeFast(p.get<double>("RiseTimeFast",-1.0)),
  fRiseTimeSlow(p.get<double>("RiseTimeSlow",-1.0)),
  fDoSlowComponent(p.get<bool>("DoSlowComponent")),
  fIncludePhotPropTimeUBSim(p.get<bool>("IncludePhotPropTimeUBSim")),
  fUsingScaleFactor(p.get<bool>("UsingScaleFactor")),
  fEDepTags(p.get< std::vector<art::InputTag> >("EDepModuleLabels")),
  fPhotonScale(p.get< std::vector<double> >("PhotonScale", std::vector<double>())),
  fPhotonEngine(art::ServiceHandle<rndm::NuRandomService>()
                ->registerAndSeedEngine(createEngine(0, "HepJamesRandom", "photon"),
                                        "HepJamesRandom", "photon", p, "SeedPhoton")),
  fScintEngine(art::ServiceHandle<rndm::NuRandomService>()
               ->registerAndSeedEngine(createEngine(0, "HepJamesRandom", "scinttime"),
                                       "HepJamesRandom", "scinttime", p, "SeedScintTime"))
{
  while(fPhotonScale.size() < fEDepTags.size())
    fPhotonScale.push_back(1.);
  produces< std::vector<sim::SimPhotons> >();

  // Parameterized Simulation
  fhicl::ParameterSet VUVTimingParams;
  fhicl::ParameterSet VISTimingParams; //this is not used in MicroBooNE because it is for the reflected light but it needs to be set here to work

  // validate configuration    
  if(fIncludePhotPropTimeUBSim && !p.get_if_present<fhicl::ParameterSet>("VUVTiming", VUVTimingParams)) 
    {
      throw art::Exception(art::errors::Configuration)
	<< "Propagation time simulation requested, but VUVTiming not specified." << "\n";
    }

  //setting reflected light which is last parameter to false makes it such that it doesn't go took for VISTimingParams
  if (fIncludePhotPropTimeUBSim)
    fPropTimeModel = std::make_unique<PropagationTimeModel>(
							    VUVTimingParams, VISTimingParams, fScintEngine, false, false);
}

double phot::UBPhotonLibraryPropagation::GetScintYield(sim::SimEnergyDeposit const& edep,
						       detinfo::LArProperties const& larp)
{
  double yieldRatio = larp.ScintYieldRatio();
  if(larp.ScintByParticleType()){
    switch(edep.PdgCode()) {
    case 2212:
      yieldRatio = larp.ProtonScintYieldRatio();
      break;
    case 13:
    case -13:
      yieldRatio = larp.MuonScintYieldRatio();
      break;
    case 211:
    case -211:
      yieldRatio = larp.PionScintYieldRatio();
      break;
    case 321:
    case -321:
      yieldRatio = larp.KaonScintYieldRatio();
      break;
    case 1000020040:
      yieldRatio = larp.AlphaScintYieldRatio();
      break;
    case 11:
    case -11:
    case 22:
      yieldRatio = larp.ElectronScintYieldRatio();
      break;
    default:
      yieldRatio = larp.ElectronScintYieldRatio();
    }
  }
  return yieldRatio;
}

void phot::UBPhotonLibraryPropagation::produce(art::Event & e)
{
  art::ServiceHandle<PhotonVisibilityService> pvs;
  art::ServiceHandle<sim::LArG4Parameters> lgpHandle;
  const detinfo::LArProperties* larp = lar::providerFrom<detinfo::LArPropertiesService>();
  auto const detprop = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob();
  
  art::ServiceHandle<art::RandomNumberGenerator> rng;  

  //auto const& module_label = moduleDescription().moduleLabel();

  CLHEP::RandPoissonQ randpoisphot{fPhotonEngine};
  CLHEP::RandFlat randflatscinttime{fScintEngine};
  
  const size_t NOpChannels = pvs->NOpChannels();
  double yieldRatio;
  double nphot,nphot_fast,nphot_slow;

  //auto fSCE = lar::providerFrom<spacecharge::SpaceChargeService>();

  sim::OnePhoton photon;
  photon.Energy = 9.7e-6;
  photon.SetInSD = false;
  
  std::unique_ptr< std::vector<sim::SimPhotons> > photCol ( new std::vector<sim::SimPhotons>);
  auto & photonCollection(*photCol);

  //size_t edep_reserve_size=0;
  std::vector< std::vector<sim::SimEnergyDeposit> const*> edep_vecs;
  for(auto label : fEDepTags){
    auto const& edep_handle = e.getValidHandle< std::vector<sim::SimEnergyDeposit> >(label);
    edep_vecs.push_back(edep_handle);
    //edep_reserve_size += edep_handle->size();
  }

  for(size_t i_op=0; i_op<NOpChannels; ++i_op){
    photonCollection.emplace_back(i_op);
  }

  size_t nvec = 0;
  for(auto const& edeps : edep_vecs){
    double scale_factor = fPhotonScale[nvec++];

    for(auto const& edep : *edeps){
      /*
	std::cout << "Processing edep with trackID=" 
	<< edep.TrackID()
	<< " pdgCode="
	<< edep.PdgCode() 
	<< " energy="
	<< edep.Energy()
	<< "(x,y,z)=("
	<< edep.X() << "," << edep.Y() << "," << edep.Z() << ")"
	<< std::endl;
      */
      double const xyz[3] = { edep.X(), edep.Y(), edep.Z() };
      
      photon.InitialPosition = TVector3(xyz[0],xyz[1],xyz[2]);
      
      geo::Point_t pos(xyz[0], xyz[1], xyz[2]);
      geo::Point_t const ScintPoint = {xyz[0], xyz[1], xyz[2]};
      phot::MappedCounts_t Visibilities = pvs->GetAllVisibilities(pos);

      if(!Visibilities)
	continue;
      
      yieldRatio = GetScintYield(edep,*larp);
      larg4::ISCalcData isdata = fISAlg.CalcIonAndScint(detprop, edep);
      nphot = isdata.numPhotons;
      nphot_fast = yieldRatio*nphot;
      nphot_slow = nphot - nphot_fast;
      
      for(size_t i_op=0; i_op<NOpChannels; ++i_op){
	if (Visibilities[i_op] < 1e-9) continue; // voxel is not visible at this optical channel
	auto nph_fast = 0.0;
	auto nph_slow = 0.0;
	if(fUsingScaleFactor) //scaling factor for increasing the light outside the TPC. This is off by default
	  {
	    nph_fast = randpoisphot.fire(nphot_fast*Visibilities[i_op]*scale_factor);
	    nph_slow = randpoisphot.fire(nphot_slow*Visibilities[i_op]*scale_factor);
	  }
	else
	  {
	    nph_fast = randpoisphot.fire(nphot_fast*Visibilities[i_op]);
	    nph_slow = randpoisphot.fire(nphot_slow*Visibilities[i_op]);
	  }
	
	//skip if there is nophotons that are visible to the OpChannels
	if(nph_fast+nph_slow == 0) continue;

	// initialize transport time distribution vector
	std::vector<double> transport_time;
	// calculate propagation times if included, does not matter whether fast or slow photon
	if (fIncludePhotPropTimeUBSim) {
	  transport_time.resize(nph_fast + nph_slow);
	  fPropTimeModel->propagationTime(transport_time, ScintPoint, i_op, false);
	}

	if(nph_fast > 0)
	  {

	    for (int i = 0; i < nph_fast; i++)
	      {
		double dTime = edep.T() + GetScintTime(fRiseTimeFast,larp->ScintFastTimeConst(),randflatscinttime);
		if(fIncludePhotPropTimeUBSim) dTime += transport_time[i];
		photon.Time = dTime;
		photonCollection[i_op].insert(photonCollection[i_op].end(),1,photon);
	      }
	  }
	if(fDoSlowComponent && nph_slow>0){

          for (int i = 0; i < nph_slow; i++)
	    {
	      double dTime = edep.T() + GetScintTime(fRiseTimeSlow,larp->ScintSlowTimeConst(),randflatscinttime);
	      if(fIncludePhotPropTimeUBSim) dTime += transport_time[i+nph_fast];
	      photon.Time = dTime;
	      photonCollection[i_op].insert(photonCollection[i_op].end(),1,photon);
	    }
	}//end doing slow component
      }//end loop over OpChannels
      
    }//end loop over edeps
  }//end loop over edep vectors
  
  e.put(std::move(photCol));
  
}

// Old Method
/*double phot::UBPhotonLibraryPropagation::GetScintTime(double scint_time, double rise_time,
      double r1, double r2)
{
  //no rise time
  if(rise_time<0.0)
    return -1 * scint_time * std::log(r1);

  while(1){
    double t = -1.0*scint_time*std::log(1-r1);
    //double g = (scint_time+rise_time)/scint_time * std::exp(-1.0*t/scint_time)/scint_time;
    if ( r2 <= (std::exp(-1.0*t/rise_time)*(1-std::exp(-1.0*t/rise_time))/scint_time/scint_time*(scint_time+rise_time)) )
      return -1 * t;
  }
  
  return 1;
  
  }*/
// New Method
// Returns the time within the time distribution of the scintillation process, when the photon was created.
// Scintillation light has an exponential decay which here is given by the decay time, tau2,
// and an exponential increase, which here is given by the rise time, tau1.
// randflatscinttime is passed to use the saved seed from the RandomNumberSaver in order to be able to reproduce the same results.
double phot::UBPhotonLibraryPropagation::single_exp(double t, double tau2)
{ 
  return exp((-1.0 * t) / tau2) / tau2;
}

double phot::UBPhotonLibraryPropagation::bi_exp(double t, double tau1, double tau2)
{
  return (((exp((-1.0 * t) / tau2) * (1.0 - exp((-1.0 * t) / tau1))) / tau2) / tau2) *
    (tau1 + tau2);
}

double phot::UBPhotonLibraryPropagation::GetScintTime(double tau1, double tau2, CLHEP::RandFlat& randflatscinttime)
{
  // tau1: rise time (originally defaulted to -1) and tau2: decay time
  //ran1, ran2 = random numbers for the algorithm
  if ((tau1 == 0.0) || (tau1 == -1.0)) { return -tau2 * log(randflatscinttime()); }
  while (1) {
    auto ran1 = randflatscinttime();
    auto ran2 = randflatscinttime();
    auto d = (tau1 + tau2) / tau2;
    auto t = -tau2 * log(1 - ran1);
    auto g = d * single_exp(t, tau2);
    if (ran2 <= bi_exp(t, tau1, tau2) / g) { return t; }
  }
}


void phot::UBPhotonLibraryPropagation::beginJob()
{

}

DEFINE_ART_MODULE(phot::UBPhotonLibraryPropagation)
