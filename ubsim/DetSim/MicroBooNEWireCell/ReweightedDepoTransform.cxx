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
#include "WireCellIface/SimpleDepo.h"

#include <memory>

WIRECELL_FACTORY(wclsReweightedDepoTransform, wcls::ReweightedDepoTransform,
		 wcls::IArtEventVisitor, WireCell::IDepoFramer)


using namespace WireCell;

wcls::ReweightedDepoTransform::ReweightedDepoTransform()
    : DepoTransform()
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
    for(Int_t xbin=1; xbin <=m_hists[0]->GetNbinsX(); xbin++){
        for(Int_t ybin=1; ybin <=m_hists[0]->GetNbinsY(); ybin++){
            for(size_t iplane=0; iplane<m_hists.size(); iplane++){
                double zcenter = m_hists[iplane]->GetXaxis()->GetBinCenter(xbin); // cm
                double ycenter = m_hists[iplane]->GetYaxis()->GetBinCenter(ybin); // cm
                double scaleMC = m_hists[iplane]->GetBinContent(xbin, ybin);
                double scaleDATA = energyCalibProvider.YZdqdxCorrection(iplane, ycenter, zcenter);
                m_hists[iplane]->SetBinContent(xbin, ybin, scaleMC/scaleDATA);
            }
        }
    }
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
        if (!h){
            THROW(IOError() << errmsg{"Cannot find hist: " + m_histnames.at(i)});
        }
        m_hists.push_back(h);
    }
}

IDepo::pointer wcls::ReweightedDepoTransform::modify_depo(WirePlaneId wpid, IDepo::pointer depo){
    Point pos = depo->pos();
    Int_t ybin = m_hists[wpid.index()]->GetYaxis()->FindBin(pos.y()/units::cm);
    Int_t zbin = m_hists[wpid.index()]->GetXaxis()->FindBin(pos.z()/units::cm);
    double scale = m_hists[wpid.index()]->GetBinContent(zbin, ybin);

    auto newdepo = std::make_shared<SimpleDepo>(depo->time(), pos, depo->charge()*scale, 
            depo, depo->extent_long(), depo->extent_tran());
    return newdepo;
}


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
