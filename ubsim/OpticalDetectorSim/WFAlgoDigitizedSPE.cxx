#ifndef WFALGODIGITIZEDSPE_CXX
#define WFALGODIGITIZEDSPE_CXX

#include "WFAlgoDigitizedSPE.h"
#include "WFAlgoUtilities.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

namespace opdet {
  
  //--------------------------------------------------------
  WFAlgoDigitizedSPE::WFAlgoDigitizedSPE() : WFAlgoSPEBase()
  //--------------------------------------------------------
  {
    Reset();
    // We initialize abnormal ch here b/c we do not want to reset every time
    // SPE is reset
    fAbnormCh = 0;
 
    //fSPETime = detinfo::DetectorClocksService::GetME().OpticalClock();
    //auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();
    //auto spe_time = ts->OpticalClock();
    //spe_time.SetTime(0);
    ::detinfo::ElecClock spe_time(0., 1600., 1000.); // 1.6ms frame period, 1GHz frequency

    std::vector<float> wf;
    SetResponseNormal_BNLv1(wf);
    SetSPE(wf,spe_time,0);
    
    SetResponseOpCh28_BNLv1(wf);
    SetSPE(wf,spe_time,28);

  }

  //------------------------------
  void WFAlgoDigitizedSPE::Reset()
  //------------------------------
  {
    WFAlgoSPEBase::Reset();
  }

  //--------------------------------------------------------------
  void WFAlgoDigitizedSPE::Process(std::vector<float> &wf,
				   const ::detinfo::ElecClock &start_time)
  //--------------------------------------------------------------
  {
    // Predefine variables to save time later
    ::detinfo::ElecClock rel_spe_start = start_time;

    rel_spe_start.SetTime(0);

    auto const& fSPETime = GetClock(OpChannel());
    auto const& fSPE = GetSPE(OpChannel());

    double unit_time = fSPETime.TickPeriod();

    auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();
    for(auto const &t : fPhotonTime) {

      //
      // Check if this photon should be added or not
      //

      // Time in electronics clock frame (with T0)
      //double time = ::detinfo::DetectorClocksService::GetME().G4ToElecTime(t);
      double time = ts->G4ToElecTime(t);

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
      rel_spe_start.SetTime(time - fSPETime.Time() - start_time.Time());
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

	rel_spe_start += unit_time;
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
  
  //----------------------------------------------------------------
  void WFAlgoDigitizedSPE::SetSPE( const std::vector<float> &wf,
				   const detinfo::ElecClock &time_info,
				   const int opch)
  //----------------------------------------------------------------
  {
    if(time_info.Time() < 0 || time_info.Ticks() >= (int)(wf.size()))
      
      throw UBOpticalException(Form("Invalid WF index (%d) for the WF of size %zu",
				    time_info.Ticks(),
				    wf.size()
				    )
			       );

    auto& fSPE = GetSPEWriteable(opch);
    auto& fSPETime = GetClockWriteable(opch);
    fSPE.clear();
    fSPE.reserve(wf.size());

    //double area = 0;
    for(auto const &v : wf) {
      fSPE.push_back(v);
      //area += v;
    }

    // Normalize area
    //for(auto &v : fSPE) v /= area;

    fSPETime = time_info;

  }

  //----------------------------------------------------------------------
  std::vector<float>& WFAlgoDigitizedSPE::GetSPEWriteable(const int opch)
  //----------------------------------------------------------------------
  { return ((opch%100) == fAbnormCh ? fSPE_OpCh28 : fSPE_Normal); }

  //-------------------------------------------------------------------------
  const std::vector<float>& WFAlgoDigitizedSPE::GetSPE(const int opch) const
  //-------------------------------------------------------------------------
  { return ((opch%100) == fAbnormCh ? fSPE_OpCh28 : fSPE_Normal); }

  //-------------------------------------------------------------------------
  ::detinfo::ElecClock& WFAlgoDigitizedSPE::GetClockWriteable(const int opch)
  //-------------------------------------------------------------------------
  { return ((opch%100) == fAbnormCh ? fSPETime_OpCh28 : fSPETime_Normal); }

  //----------------------------------------------------------------------------
  const ::detinfo::ElecClock& WFAlgoDigitizedSPE::GetClock(const int opch) const
  //----------------------------------------------------------------------------
  { return ((opch%100) == fAbnormCh ? fSPETime_OpCh28 : fSPETime_Normal); }

}

#endif
