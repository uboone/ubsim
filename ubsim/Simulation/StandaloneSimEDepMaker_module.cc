////////////////////////////////////////////////////////////////////////
// Class:       StandaloneSimEDepMaker
// Plugin Type: producer (art v3_01_02)
// File:        StandaloneSimEDepMaker_module.cc
//
// Generated at Sat Jun 29 16:43:37 2019 by Wesley Ketchum using cetskelgen
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

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"

#include <memory>

#include "lardataobj/Simulation/SimEnergyDeposit.h"

namespace sim {
  class StandaloneSimEDepMaker;
}


class sim::StandaloneSimEDepMaker : public art::EDProducer {
public:
  explicit StandaloneSimEDepMaker(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  StandaloneSimEDepMaker(StandaloneSimEDepMaker const&) = delete;
  StandaloneSimEDepMaker(StandaloneSimEDepMaker&&) = delete;
  StandaloneSimEDepMaker& operator=(StandaloneSimEDepMaker const&) = delete;
  StandaloneSimEDepMaker& operator=(StandaloneSimEDepMaker&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  typedef enum Mode{
    kSingle = 0,
    kLine = 1
  } Mode_t;

private:

  Mode_t fMode;
  
  int fNumPhotons;
  int fNumElectrons;
  double fEnergy;
  geo::Point_t fXYZ;
  double fTime;

  bool fDoGrid;
  unsigned int fCryostat,fTPC;
  float fXFiducialCut;
  float fYFiducialCut;
  float fZFiducialCut;
  unsigned int fNGridX;
  unsigned int fNGridY;
  unsigned int fNGridZ;
  bool fVerbose;

  sim::SimEnergyDeposit MakeSingle(geo::Point_t);
  void FillGridPoints();

  std::vector<geo::Point_t> fGridPoints;
  size_t fGridIndex;

};


sim::StandaloneSimEDepMaker::StandaloneSimEDepMaker(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fMode(Mode_t(p.get<int>("Mode"))),
  fNumPhotons(p.get<int>("NumPhotons")),
  fNumElectrons(p.get<int>("NumElectrons")),
  fEnergy(p.get<double>("Energy")),
  fXYZ({p.get<float>("X",0.0),p.get<float>("Y",0.0),p.get<float>("Z",0.0)}),
  fTime(p.get<double>("Time")),
  fDoGrid(p.get<bool>("DoGrid")),
  fCryostat(p.get<unsigned int>("Cryostat",0)),
  fTPC(p.get<unsigned int>("TPC",0)),
  fXFiducialCut(p.get<float>("XFiducialCut",0.0)),
  fYFiducialCut(p.get<float>("YFiducialCut",0.0)),
  fZFiducialCut(p.get<float>("ZFiducialCut",0.0)),
  fNGridX(p.get<unsigned int>("NGridX",1)),
  fNGridY(p.get<unsigned int>("NGridY",1)),
  fNGridZ(p.get<unsigned int>("NGridZ",1)),
  fVerbose(p.get<bool>("Verbose",false))

  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  produces< std::vector<sim::SimEnergyDeposit> >();

  fGridIndex=0;
  if(fDoGrid){
    FillGridPoints();
  }
}

void sim::StandaloneSimEDepMaker::FillGridPoints()
{

  fGridPoints.clear();

  art::ServiceHandle<geo::Geometry> geo;
  auto const& bbox = geo->TPC(fCryostat,fTPC).ActiveBoundingBox();
  
  geo::BoxBoundedGeo bbox_fid(bbox.MinX()+fXFiducialCut,
			      bbox.MaxX()-fXFiducialCut,
			      bbox.MinY()+fYFiducialCut,
			      bbox.MaxY()-fYFiducialCut,
			      bbox.MinZ()+fZFiducialCut,
			      bbox.MaxZ()-fZFiducialCut);
  if(fVerbose)
    std::cout << bbox_fid.MinX() << " " << bbox_fid.MaxX() << " " 
	      << bbox_fid.MinY() << " " << bbox_fid.MaxY() << " "
	      << bbox_fid.MinZ() << " " << bbox_fid.MaxZ() << " "
	      << std::endl;

  float x_incr=bbox_fid.SizeX()/(fNGridX+1);
  float y_incr=bbox_fid.SizeY()/(fNGridY+1);
  float z_incr=bbox_fid.SizeZ()/(fNGridZ+1);

  for(size_t i_x=1; i_x<fNGridX+1; ++i_x){
    float xval = bbox_fid.MinX()+i_x*x_incr;
    for(size_t i_y=1; i_y<fNGridY+1; ++i_y){
      float yval = bbox_fid.MinY()+i_y*y_incr;
      for(size_t i_z=1; i_z<fNGridZ+1; ++i_z){
	float zval = bbox_fid.MinZ()+i_z*z_incr;
	if(fVerbose){
	  std::cout << "Making point (x,y,z)=(" 
		    << xval << ","
		    << yval << ","
		    << zval << ")"
		    << std::endl;
	}
	fGridPoints.emplace_back(xval,yval,zval);
      }
    }
  }


  for(float x = bbox_fid.MinX(); bbox_fid.ContainsX(x); x+=bbox_fid.SizeX()/(fNGridX+1)){
    for(float y = bbox_fid.MinY(); bbox_fid.ContainsY(y); y+=bbox_fid.SizeY()/(fNGridY+1)){
      for(float z = bbox_fid.MinZ(); bbox_fid.ContainsZ(z); z+=bbox_fid.SizeZ()/(fNGridZ+1)){

	fGridPoints.emplace_back(x,y,z);
      }
    }
  }
  
}

sim::SimEnergyDeposit sim::StandaloneSimEDepMaker::MakeSingle(geo::Point_t pt)
{
  if(fVerbose)
    std::cout << "Making edep at point (x,y,z)=(" 
	      << pt.X() << ","
	      << pt.Y() << ","
	      << pt.Z() << ")"
	      << std::endl;

  return sim::SimEnergyDeposit(fNumPhotons,
			       fNumElectrons,
			       fEnergy,
			       pt,
			       pt,
			       fTime,
			       fTime,
			       -1, //invalid track id
			       0); //empty pdg code
}

void sim::StandaloneSimEDepMaker::produce(art::Event& e)
{
  std::unique_ptr< std::vector<sim::SimEnergyDeposit> > edepColPtr(new std::vector<sim::SimEnergyDeposit>);
  
  //if doing a grid
  if(fDoGrid){
    if(fMode==Mode_t::kSingle){
      
      if(fVerbose)
	std::cout << "On grid point " << fGridIndex+1 << " of " << fGridPoints.size() << std::endl;

      edepColPtr->push_back(MakeSingle(fGridPoints[fGridIndex]));

      //increment index
      ++fGridIndex;
      if(fGridIndex>=fGridPoints.size()) fGridIndex=0;
    }
  }
  else{
    if(fMode==Mode_t::kSingle){
      edepColPtr->push_back(MakeSingle(fXYZ));
    }
  }


  e.put(std::move(edepColPtr));

}

DEFINE_ART_MODULE(sim::StandaloneSimEDepMaker)
