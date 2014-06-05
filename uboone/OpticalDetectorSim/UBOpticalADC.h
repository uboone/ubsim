/**
 * \file UBOpticalADC.h
 *
 * \ingroup OpticalDetectorSim
 * 
 * \brief Class def header for a class UBOpticalADC
 *
 * @author kazuhiro
 */

/** \addtogroup OpticalDetectorSim

    @{*/
#ifndef UBOPTICALADC_H
#define UBOPTICALADC_H

#include "UBADCBase.h"

namespace opdet {
  /**
     \class UBOpticalADC
     A generator class that handles WF generation algorithm to create an waveform.
  */
  class UBOpticalADC : public UBADCBase {
    
  public:
    
    /// Default constructor
    UBOpticalADC();
    
    /// Default destructor
    virtual ~UBOpticalADC(){}

    /// Function to reset algorithm configuration
    virtual void Reset();

    /// Function to enable gain/T0 spread
    void EnableSpread(bool doit=true) { fSPE.EnableSpread(doit); }

    /// Function to add G4 photon in G4 time (ns as that is G4 natural unit)
    void AddPhoton(double g4time){ fInputPhotonTime.push_back(g4time); }

    /// Function to set G4 photons in G4 time (ns as that is G4 natural unit)
    void SetPhotons(const std::vector<double>& g4time);

    /// Method to generate waveform for a specific channel
    void GenWaveform(const unsigned int ch, 
		     std::vector<unsigned short> &high_wf,
		     std::vector<unsigned short> &low_wf  );

    /// Method to retrieve signal G4 photon time vector
    const std::vector<double>& SignalPhotonTime() const { return fPhotonTime; }

    /// Method to retrieve dark-noise photon time vector
    const std::vector<double>& DarkPhotonTime() const { return fDarkPhotonTime; }
    
  protected:

    /// Function to generate dark noise photon timings (called in GenWaveform)
    void GenDarkNoise(double dark_rate,double period);

  protected:

    /// G4 photon time for signal in G4 clock
    std::vector<double> fInputPhotonTime;

    /// Dark noise photon time in G4 clock
    std::vector<double> fDarkPhotonTime;

    /// Photon time that is injected to the waveform (after QE applied)
    std::vector<double> fPhotonTime;

    /// Algorithm to generate SPE waveform
    WFAlgoAnalyticalSPE fSPE;

    /// Algorithm to generate Pedestal
    WFAlgoPedestal      fPED;

  };
}

#endif
/** @} */ // end of doxygen group 

