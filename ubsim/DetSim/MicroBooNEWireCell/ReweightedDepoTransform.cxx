#include "ReweightedDepoTransform.h"

#include "larevt/CalibrationDBI/Interface/ElectronicsCalibService.h"
#include "larevt/CalibrationDBI/Interface/ElectronicsCalibProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "ubevt/Database/TPCEnergyCalib/TPCEnergyCalibService.h"
#include "ubevt/Database/TPCEnergyCalib/TPCEnergyCalibProvider.h"
#include "larcore/Geometry/Geometry.h"

#include "TFile.h"
#include "TH2F.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Persist.h"
#include "WireCellAux/SimpleDepo.h"

#include <memory>
#include <math.h>

WIRECELL_FACTORY(wclsReweightedDepoTransform, wcls::ReweightedDepoTransform,
		 wcls::IArtEventVisitor, WireCell::IDepoFramer)


using namespace WireCell;

wcls::ReweightedDepoTransform::ReweightedDepoTransform()
    : DepoTransform()
    , m_scaleDATA_perplane{1.0, 1.0, 1.0}
    , m_scaleMC_perplane{1.0, 1.0, 1.0}
{
}

wcls::ReweightedDepoTransform::~ReweightedDepoTransform()
{
}

void wcls::ReweightedDepoTransform::visit(art::Event & event)
{
    /// post configuration
    /// how to access energy calibration database 
    std::cout << "wclsReweightedDepoTransform: DATA YZ correction map read in!\n"; 
    const lariov::TPCEnergyCalibProvider& energyCalibProvider
        = art::ServiceHandle<lariov::TPCEnergyCalibService>()->GetProvider();

    /// energyCalibProvider.YZdqdxCorrection(view, y, z); 
    //view: plane 0, 1, 2
    //y, z: yz positions in units of cm
    for(Int_t xbin=1; xbin <=m_hists_orig[0]->GetNbinsX(); xbin++){
        for(Int_t ybin=1; ybin <=m_hists_orig[0]->GetNbinsY(); ybin++){
            for(size_t iplane=0; iplane<m_hists_orig.size(); iplane++){
                double zcenter = m_hists_orig[iplane]->GetXaxis()->GetBinCenter(xbin); // cm
                double ycenter = m_hists_orig[iplane]->GetYaxis()->GetBinCenter(ybin); // cm
                double scaleMC = m_hists_orig[iplane]->GetBinContent(xbin, ybin);
                double scaleDATA = energyCalibProvider.YZdqdxCorrection(iplane, ycenter, zcenter);


		//we want a scale to apply to compensate for differnces in relative scale of data and MC adc
		//this is a scale to make the mc look more like data
		//would multiply this by an MC scale
		//these are adc/electron, from stopping muon calibration
		//(from fcl, CalAreaConstants per plane)
		double scale_gain_mc2data = m_scaleMC_perplane[iplane]/m_scaleDATA_perplane[iplane];
		if(std::isnan(scale_gain_mc2data) || std::isinf(scale_gain_mc2data) || 
		   scale_gain_mc2data<0.0001 || scale_gain_mc2data>1000.){
		    std::cout << "WARNING! mc2data GAIN SCALE OUT OF RANGE! " << scale_gain_mc2data 
			      << " ... setting to 1.0" << std::endl;
		    scale_gain_mc2data = 1.0;
		}
		//std::cout << "Using ... scale_gain_mc2data = " << scale_gain_mc2data << std::endl;

		//Wes, 11 March 2019
		//adding some protection in the case that scaleMC or scaleDATA are unreasonable ranges
		//Saying reasonable is <0.0001 or >1000
		//In either case, set the ratio to 1.

		double scale_corr = 1.0;
		if(scaleMC>0.0001 && scaleMC<1000. && scaleDATA>0.0001 && scaleDATA<1000.){
		    scale_corr = scaleMC/scaleDATA;
		    //scale_corr = scaleDATA/scaleMC;
		}
		else {
		    std::cout << "WARNING! for xbin " << xbin 
			      << " and ybin " << ybin
			      << " we had scaleMC=" << scaleMC 
			      << " and scaleDATA=" << scaleDATA
			      << " so we're setting correction scale to " << scale_corr << std::endl;
		}
		
		if(std::isnan(scale_corr) || std::isinf(scale_corr)){
		    std::cout << "WARNING! for xbin " << xbin 
			      << " and ybin " << ybin
			      << " we had scaleMC=" << scaleMC 
			      << " and scaleDATA=" << scaleDATA
			      << " and somehow scale ratio was still a nan/inf so setting to 1.0" << std::endl;
		    scale_corr = 1.0;
		}
		    
		//scale_gain_mc2data = 1.0;
		//scale_corr = 1.0;
  
		//divide by the mc2data scale, since it would be multiplied in denominator of scale_corr
                m_hists[iplane]->SetBinContent(xbin, ybin, scale_corr/scale_gain_mc2data);
		//std::cout << "final scale appled was " << scale_corr << "/" << scale_gain_mc2data << "="
		//	  << scale_corr/scale_gain_mc2data << std::endl;


		/*
		std::cout << iplane << " " << xbin << " " << ybin << " " << zcenter << " " << ycenter << " "
			  << m_scaleMC_perplane[iplane] << " " << m_scaleDATA_perplane[iplane] << " "
			  << scale_gain_mc2data << " "
			  << scaleMC << " " << scaleDATA << " " << scale_corr << " " << 1./scale_corr << " "
			  << scale_corr/scale_gain_mc2data << " "
			  << m_hists[iplane]->GetBinContent(xbin,ybin) << std::endl;
		*/

            }
        }
    }

    /// m_scaleDATA_perplane[iplane] can be used
    /// m_scaleMC_perplane[iplane] can be used
}


void wcls::ReweightedDepoTransform::configure(const WireCell::Configuration& cfg)
{
    DepoTransform::configure(cfg);
    std::cout << "wclsReweightedDepoTransform: MC YZ correction histograms read in!\n"; 
    m_fileName_MC = cfg["filenameMC"].asString();
    if (m_fileName_MC.empty()){
        THROW(ValueError() << errmsg{"Must configure a filenameMC!"});
    }
    auto histnames = cfg["histnames"]; // plane 0, 1, 2
    for (auto histname : histnames) {
        m_histnames.push_back(histname.asString());
    }
    if (m_histnames.size() !=3){
        THROW(ValueError() << errmsg{"Size of MC correction maps/hists should be 3!"});
    }

    std::string filepath = Persist::resolve(m_fileName_MC);
    if (filepath.empty()){
        THROW(IOError() << errmsg{"Cannot find the file: " + m_fileName_MC});
    }
    TFile* f = new TFile(filepath.c_str());
    for(size_t i=0; i<m_histnames.size(); i++)
    {
        TH2F* h = (TH2F*)f->Get(m_histnames.at(i).c_str());
	TH2F* h1 = (TH2F*)h->Clone("h1");
        if (!h){
            THROW(IOError() << errmsg{"Cannot find hist: " + m_histnames.at(i)});
        }
        m_hists_orig.push_back(h);
        m_hists.push_back(h1);
    }
   
    auto scaleDATA_planes = cfg["scaleDATA_perplane"]; 
    if (scaleDATA_planes.empty()){
        std::cout << "wclsReweightedDepoTransform: DATA per plane scaling factors not found and default (1.0) will be used!\n";
    }
    else{
        m_scaleDATA_perplane.clear();
        for (auto scale_perplane : scaleDATA_planes) {
            m_scaleDATA_perplane.push_back(scale_perplane.asDouble());
        }
    }
    std::cout << "DATA scaling per plane: "<< m_scaleDATA_perplane[0] << " " << m_scaleDATA_perplane[1] << " " << m_scaleDATA_perplane[2] << std::endl;

    auto scaleMC_planes = cfg["scaleMC_perplane"]; 
    if (scaleMC_planes.empty()){
        std::cout << "wclsReweightedDepoTransform: MC per plane scaling factors not found and default (1.0) will be used!\n";
    }
    else{
        m_scaleMC_perplane.clear();
        for (auto scale_perplane : scaleMC_planes) {
            m_scaleMC_perplane.push_back(scale_perplane.asDouble());
        }
    }
    std::cout << "MC scaling per plane: "<< m_scaleMC_perplane[0] << " " << m_scaleMC_perplane[1] << " " << m_scaleMC_perplane[2] << std::endl;

}

IDepo::pointer wcls::ReweightedDepoTransform::modify_depo(WirePlaneId wpid, IDepo::pointer depo){
    Point pos = depo->pos();
    Int_t ybin = m_hists[wpid.index()]->GetYaxis()->FindBin(pos.y()/units::cm);
    Int_t zbin = m_hists[wpid.index()]->GetXaxis()->FindBin(pos.z()/units::cm);
    double scale = m_hists[wpid.index()]->GetBinContent(zbin, ybin);

    /*
    std::cout << "Edep scaling: " << depo->charge() << " at (" 
	      << pos.x()/units::cm << "," << pos.y()/units::cm << "," << pos.z()/units::cm << ") "
	      << "ybin=" << ybin << " zbin=" << zbin 
	      << "wpid=" << wpid.index() << " scale=" << scale
	      << " final_charge=" << depo->charge()*scale << std::endl;
    */

    //added by Wes to just make sure we don't do insane scaling...
    if(scale<0.0001 || scale>1000. || std::isnan(scale) || std::isinf(scale) ){

	std::cout << "WARNING! Found crazy scale " << scale 
		  << " in modify_depo for pos (y,z) " 
		  << pos.y()/units::cm << "," << pos.z()/units::cm << ")"
		  << " --> Setting it to 1.0" << std::endl;
	scale=1.;
    }



    using WireCell::Aux::SimpleDepo;
    auto newdepo = std::make_shared<SimpleDepo>(depo->time(), pos, depo->charge()*scale, 
            depo, depo->extent_long(), depo->extent_tran());
    return newdepo;
}


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
