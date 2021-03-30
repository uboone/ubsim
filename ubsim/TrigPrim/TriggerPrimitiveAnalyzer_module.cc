////////////////////////////////////////////////////////////////////////
// Class:       TriggerPrimitiveAnalyzer
// Plugin Type: analyzer (art v2_11_02)
// File:        TriggerPrimitiveAnalyzer_module.cc
//
// Generated at Mon Jun  3 13:29:25 2019 by Claire Hinrichs using cetskelgen and modified by Mark Ross-Lonergan and Harrison Seigal
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

//#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "TH1D.h"
#include <string>
#include "TCanvas.h"
#include <iostream>
#include <memory>
#include <cmath>
#include <algorithm>
#include "TF1.h"
#include "TTree.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "TH2D.h"
#include "TF2.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "larreco/HitFinder/HitFinderTools/ICandidateHitFinder.h"
#include "larreco/HitFinder/HitFinderTools/IPeakFitter.h"
#include "TGraphErrors.h"
#include "TDecompSVD.h"
#include "TMath.h"
#include "TStopwatch.h"

class TriggerPrimitiveAnalyzer : public art::EDAnalyzer {
    public:

        explicit TriggerPrimitiveAnalyzer(fhicl::ParameterSet const & pset);

        // The compiler-generated destructor is fine for non-base
        // classes without bare pointers or other resource use.

        // Plugins should not be copied or assigned.
        TriggerPrimitiveAnalyzer(TriggerPrimitiveAnalyzer const &) = delete;
        TriggerPrimitiveAnalyzer(TriggerPrimitiveAnalyzer &&) = delete;
        TriggerPrimitiveAnalyzer & operator = (TriggerPrimitiveAnalyzer const &) = delete;
        TriggerPrimitiveAnalyzer & operator = (TriggerPrimitiveAnalyzer &&) = delete;

        // Required functions.
        void analyze(art::Event const & e) override;

        // Selected optional functions.
        void beginJob() override;

        void endJob() override;

    private:




        TTree *event_tree;

        int fPrimMode;
        std::string fAllHitsInstanceName;
        std::string fWireModuleLabel;
        void reconfigure(fhicl::ParameterSet const& p) ;

        std::vector<int> m_channel;
        std::vector<int> m_view;
        std::vector<float> m_max_ADC;
        std::vector<float> m_max_ADC_tick;
        std::vector<float> m_integral_sum;
        std::vector<float> m_tot;
        std::vector<float> m_first_tick; 
        std::vector<float> m_integral_over_n;


};


TriggerPrimitiveAnalyzer::TriggerPrimitiveAnalyzer(fhicl::ParameterSet const & pset): EDAnalyzer(pset)
                                                                                      // :
                                                                                      // Initialize member data here.
{
    this->reconfigure(pset);

    return;
}

void TriggerPrimitiveAnalyzer::reconfigure(fhicl::ParameterSet const& p){
    fPrimMode = p.get<int>("PrimMode");
    fAllHitsInstanceName = p.get<std::string>("AllHitsInstanceName","");
    fWireModuleLabel = p.get<std::string>("WireModuleLabel","compress");
    return;
}


void TriggerPrimitiveAnalyzer::beginJob()
{

    art::ServiceHandle<art::TFileService const> tfs;

    event_tree = tfs->make<TTree>("event_tree","event_tree");
    event_tree->Branch("channel",&m_channel);
    event_tree->Branch("view",&m_view);
    event_tree->Branch("max_ADC",&m_max_ADC);
    event_tree->Branch("max_ADC_tick",&m_max_ADC_tick);
    event_tree->Branch("integral_sum",&m_integral_sum);
    event_tree->Branch("tot",&m_tot);;
    event_tree->Branch("first_tick",&m_first_tick);
    event_tree->Branch("integral_over_n",&m_integral_over_n);


}


void TriggerPrimitiveAnalyzer::analyze(art::Event const & e){

    //std::cout<<"Starting: "<<std::endl;
    art::ValidHandle<std::vector<recob::Wire>> const & wiredata = e.getValidHandle<std::vector<recob::Wire>>(fWireModuleLabel);

    //std::cout<<"Iterating over WireData: "<<wiredata->size()<<std::endl;

    for(size_t rdIter = 0; rdIter < wiredata->size(); ++rdIter){
        //std::cout<<rdIter<<" / "<<wiredata->size()<<std::endl;

        // get the reference to the current non-deconvolved recob::Wire                                                                                               
        art::Ptr<recob::Wire> wireVec(wiredata, rdIter);
        //art::Ptr<raw::RawDigit> rawdigits = RawDigits.at(rdIter);
        auto channel = wireVec->Channel();
        auto zsROIs = wireVec->SignalROI();
        art::ServiceHandle<geo::Geometry const> geom;
        //std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
        //geo::WireID wid  = wids[0];
        //geo::PlaneID::PlaneID_t plane = wid.Plane;

        std::vector<recob::Hit>  filteredHitVec;


        for (auto iROI = zsROIs.begin_range(); iROI != zsROIs.end_range(); ++iROI) {
            auto ROI = *iROI;
            const size_t firstTick = ROI.begin_index();
            const size_t endTick = ROI.end_index();

            //std::string titlestring = "roi_original"+ std::to_string(channel);
            //TH1D horig(titlestring.c_str(), "roi_original;Tick;ADC", endTick + 1 - firstTick, firstTick, endTick + 1);
            //horig.SetLineColor(kBlack);

            float integralsum = 0; 
            float maxpeak = 0;
            //float  tot = 0;
            float peaktime = 0;
            float integralN=0;

            for (size_t iTick = ROI.begin_index(); iTick < ROI.end_index(); iTick++ ){
                //std::cout<<"iTick = "<<iTick<<std::endl;
                //std::cout<<"ROI[iTick] = "<<ROI[iTick]<<std::endl;
                integralsum +=  std::abs(ROI[iTick]);
                if (std::abs(ROI[iTick]) > maxpeak){
                    maxpeak = std::abs(ROI[iTick]);
                    peaktime = iTick;
                    //std::cout<<"maxpeak"<<maxpeak<<std::endl;
                    //std::cout<<"peaktime"<<peaktime<<std::endl;
                }
            
                if(iTick-firstTick<=12) integralN += std::abs(ROI[iTick]);
            }


            /*
            std::cout<<wid<<std::endl;
            std::cout<<wireVec->NSignal()<<" "<<wireVec->Signal().size()<<std::endl;
            for(auto &s: wireVec->Signal()) if(s!=0)std::cout<<s<<" ";
            std::cout<<std::endl;
            */

             m_channel.push_back((int)channel);
             m_view.push_back((int)wireVec->View());
             m_max_ADC.push_back((float)maxpeak);
             m_max_ADC_tick.push_back((float)peaktime);

             m_integral_sum.push_back((float)integralsum);
             m_tot.push_back((float)(endTick-firstTick));
             m_first_tick.push_back((float)firstTick); 
             m_integral_over_n.push_back((float)integralN);


        } // end of roi loop

    } // end of channel loop

    event_tree->Fill();
}


void TriggerPrimitiveAnalyzer::endJob()
{
    // Implementation of optional member function here.
}

DEFINE_ART_MODULE(TriggerPrimitiveAnalyzer)
