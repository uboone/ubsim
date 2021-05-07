////////////////////////////////////////////////////////////////////////
// Class:       TriggerPrimitiveAnalyzer
// Plugin Type: analyzer (art v2_11_02)
// File:        TriggerPrimitiveAnalyzer_module.cc
//
// Generated at Mon Jun  3 13:29:25 2019 by Claire Hinrichs using cetskelgen and modified by Mark Ross-Lonergan and Harrison Seigal
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

//#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
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
        std::string m_geantModuleLabel;
        void reconfigure(fhicl::ParameterSet const& p) ;

        std::vector<int> m_channel;
        std::vector<int> m_view;
        std::vector<float> m_max_ADC;
        std::vector<float> m_max_ADC_tick;
        std::vector<float> m_integral_sum;
        std::vector<float> m_tot;
        std::vector<float> m_first_tick; 
        std::vector<float> m_integral_over_n;
        std::vector<float> m_MCmuon_init_px;
        std::vector<float> m_MCmuon_init_py;
        std::vector<float> m_MCmuon_init_pz;
        std::vector<float> m_MCmuon_init_x;
        std::vector<float> m_MCmuon_init_y;
        std::vector<float> m_MCmuon_init_z;
        std::vector<float> m_MCmuon_init_E;
        std::vector<float> m_run;
        std::vector<float> m_subrun;


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
    m_geantModuleLabel = p.get<std::string>("MCParticleModuleLabel","largeant");
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
    event_tree->Branch("muon_init_px",&m_MCmuon_init_px);
    event_tree->Branch("muon_init_py",&m_MCmuon_init_py);
    event_tree->Branch("muon_init_pz",&m_MCmuon_init_pz);
    event_tree->Branch("muon_init_x",&m_MCmuon_init_x);
    event_tree->Branch("muon_init_y",&m_MCmuon_init_y);
    event_tree->Branch("muon_init_z",&m_MCmuon_init_z);
    event_tree->Branch("muon_init_E",&m_MCmuon_init_E);
    event_tree->Branch("run",&m_run);
    event_tree->Branch("subrun",&m_subrun);

}


void TriggerPrimitiveAnalyzer::analyze(art::Event const & e){

    //std::cout<<"Starting: "<<std::endl;
    art::ValidHandle<std::vector<recob::Wire>> const & wiredata = e.getValidHandle<std::vector<recob::Wire>>(fWireModuleLabel);

    std::vector<art::Ptr<simb::MCParticle>> mcParticleVector;//loading MCParticles
    art::ValidHandle<std::vector<simb::MCParticle>> const & mcParticleHandle= e.getValidHandle<std::vector<simb::MCParticle>>(m_geantModuleLabel);
    art::fill_ptr_vector(mcParticleVector,mcParticleHandle);

    //auto run = e.getRun();
    //auto subrun = e.getSubRun();
    //std::cout<<"run = "<<run<<std:endl;
    //std::cout<<"subrun = "<<subrun<<std:endl;

    //std::cout<<"mcParticleVector.size() = "<<mcParticleVector.size()<<std::endl;
    m_run.push_back((int)e.run());
    m_subrun.push_back((int)e.subRun());
    for(size_t j=0;j< mcParticleVector.size();j++){ 
        const art::Ptr<simb::MCParticle> mcp = mcParticleVector[j]; 
 	if (mcp->PdgCode()==13 or mcp->PdgCode()==-13){
		m_MCmuon_init_E.push_back((float)mcp->E()); 
    		m_MCmuon_init_px.push_back((float)mcp->Px()); 
    		m_MCmuon_init_py.push_back((float)mcp->Py()); 
    		m_MCmuon_init_pz.push_back((float)mcp->Pz()); 
    		m_MCmuon_init_x.push_back((float)mcp->Position()[0]); 
    		m_MCmuon_init_y.push_back((float)mcp->Position()[1]); 
    		m_MCmuon_init_z.push_back((float)mcp->Position()[2]); 
    		//std::cout<<"PARG: "<<j<<" PDG "<<mcp->PdgCode()<<" EndProcess: "<<mcp->EndProcess()<<std::endl; 
    		//std::cout<<"Energy"<<mcp->E()<<std::endl; 
		//std::cout<<"X"<<mcp->Position()[0]<<std::endl;
		//std::cout<<"Y"<<mcp->Position()[1]<<std::endl;
		//std::cout<<"Z"<<mcp->Position()[2]<<std::endl;
		//std::cout<<"Px"<<mcp->Px()<<std::endl;
		//std::cout<<"Py"<<mcp->Py()<<std::endl;
		//std::cout<<"Pz"<<mcp->Pz()<<std::endl;
    	}
    }
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

    m_channel.clear();
    m_view.clear();
    m_max_ADC.clear();
    m_max_ADC_tick.clear();
    m_integral_sum.clear();
    m_tot.clear();
    m_first_tick.clear(); 
    m_integral_over_n.clear();

    m_MCmuon_init_E.clear(); 
    m_MCmuon_init_px.clear(); 
    m_MCmuon_init_py.clear(); 
    m_MCmuon_init_pz.clear(); 
    m_MCmuon_init_x.clear(); 
    m_MCmuon_init_y.clear(); 
    m_MCmuon_init_z.clear(); 
    m_run.clear();
    m_subrun.clear();

}


void TriggerPrimitiveAnalyzer::endJob()
{
    // Implementation of optional member function here.
}

DEFINE_ART_MODULE(TriggerPrimitiveAnalyzer)
