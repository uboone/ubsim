#ifndef WFALGODIGITIZEDSPE_CXX
#define WFALGODIGITIZEDSPE_CXX

#include "WFAlgoDigitizedSPE.h"
#include "WFAlgoUtilities.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

namespace {
  auto normal_response()
  {
    std::vector<float> wf;
    opdet::SetResponseNormal_BNLv1(wf);
    return wf;
  }
    
  auto opch28_response()
  {
    std::vector<float> wf;
    opdet::SetResponseNormal_BNLv1(wf);
    return wf;
  }
    
  auto default_clock(const std::vector<float> &wf)
  {
    detinfo::ElecClock const result{0., 1600., 1000.}; // 1.6ms frame period, 1GHz frequency
    if(result.Ticks() >= static_cast<int>(wf.size()))
      throw opdet::UBOpticalException(Form("Invalid WF index (%d) for the WF of size %zu",
                                           result.Ticks(),
                                           wf.size()
                                           )
                                      );
    return result;
  }
  }

namespace opdet {

  //--------------------------------------------------------
  WFAlgoDigitizedSPE::WFAlgoDigitizedSPE()
    : fSPE_Normal{normal_response()}
    , fSPE_OpCh28{opch28_response()}
    , fSPETime_Normal{default_clock(fSPE_Normal)}
    , fSPETime_OpCh28{default_clock(fSPE_OpCh28)}
  {}

  //------------------------------
  void WFAlgoDigitizedSPE::Reset()
  //------------------------------
  {
    WFAlgoSPEBase::Reset();
  }

  //--------------------------------------------------------------
  void WFAlgoDigitizedSPE::Process(std::vector<float> &wf,
                                   const detinfo::DetectorClocksData& clockData,
				   const ::detinfo::ElecClock &start_time)
  //--------------------------------------------------------------
  {
    // Predefine variables to save time later
    ::detinfo::ElecClock rel_spe_start = start_time.WithTime(0);

    auto const& fSPETime = GetClock(OpChannel());
    auto const& fSPE = GetSPE(OpChannel());

    double unit_time = fSPETime.TickPeriod();

    for(auto const &t : fPhotonTime) {

      //
      // Check if this photon should be added or not
      //

      // Time in electronics clock frame (with T0)
      //double time = ::detinfo::DetectorClocksService::GetME().G4ToElecTime(t);
      double time = clockData.G4ToElecTime(t);

      if(fEnableSpread)  time +=  RandomServer::GetME().Gaus(fT0,fT0Sigma) * 1.e-3 ;
      else time += fT0 * 1.e-3;

      // If before waveform vector, ignore
      if(time < start_time.Time()) continue;

      // If after waveform vector, ignore
      if(time > (start_time.Time() + start_time.Time((int)(wf.size())))) continue;

      //
      // Add signal
      //

      // Figure out time stamp of the beginning of SPE
      rel_spe_start = rel_spe_start.WithTime(time - fSPETime.Time() - start_time.Time());
      //float maxval = 0;
      //int   argmax = 0;
      //int start = rel_spe_start.Ticks();
      //std::cout<<fSPE[51]<<std::endl;

      auto thisgain = RandomServer::GetME().Gaus(fGain,fGainSigma*fGain);

      for(size_t i=0; i < fSPE.size(); ++i ) {

	if(rel_spe_start.Ticks() >= (int)(wf.size())) break;

	if(rel_spe_start.Ticks() >= 0) {

	  float val = 0.;

	  if(!fEnableSpread) 

	    val = (fGain * fSPE.at(i));

	  else

	    val = (  thisgain * fSPE.at(i) );
	    
	  wf.at(rel_spe_start.Ticks()) += val;
	  /*
	  if(val > maxval) {
	    maxval = val;
	    argmax = rel_spe_start.Ticks() - start;
	  }
	  */
	}

        rel_spe_start = rel_spe_start.AdvanceTimeBy(unit_time);
	if( (unit_time * i) > 300 && fabs(fGain * fSPE.at(i)) < 1)
	  break;
      }
      /*
      std::cout << "OpCh: " << OpChannel() << std::endl;
      std::cout << "  Max Value: " << maxval << " (gain) " << fGain << " added @ " << argmax << " (local vector)" << std::endl;
      maxval = 0;
      argmax = 0;
      for(int i=start; i<(start+20); ++i) {
	if(wf[i] <= maxval) continue;
	maxval = wf[i];
	argmax = i;
      }
      std::cout << "  Current local max: " << maxval << " @ " << argmax << " (local vector)" << std::endl;
      */
    }

  }
  
  //-------------------------------------------------------------------------
  const std::vector<float>& WFAlgoDigitizedSPE::GetSPE(const int opch) const
  //-------------------------------------------------------------------------
  { return ((opch%100) == 28 ? fSPE_OpCh28 : fSPE_Normal); }

  //----------------------------------------------------------------------------
  const ::detinfo::ElecClock& WFAlgoDigitizedSPE::GetClock(const int opch) const
  //----------------------------------------------------------------------------
  { return ((opch%100) == 28 ? fSPETime_OpCh28 : fSPETime_Normal); }

}

#endif
