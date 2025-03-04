/**
 * \file WFAlgoBase.h
 *
 * \ingroup OpticalDetector
 * 
 * \brief Class def header for a class WFAlgoBase
 *
 * @author kazuhiro
 */

/** \addtogroup OpticalDetectorSim

    @{*/
#ifndef WFALGOBASE_H
#define WFALGOBASE_H

#include <vector>
#include "lardataalg/DetectorInfo/ElecClock.h"
#include "UBOpticalException.h"
namespace detinfo {
  class DetectorClocksData;
}

namespace opdet {
  /**
     \class WFAlgoBase
     Base class for WF generation algorithms for optical detector
  */
  class WFAlgoBase{
    
  public:
    
    /// Default constructor
    WFAlgoBase() : fOpChannel(-1) {}

    /// Default destructor
    virtual ~WFAlgoBase(){}

    /**
       Core function: the algorithm is supposed to add relevant singal/noise ADC
       values to the input vector std::vector<float> wf. The tick=0 timing is 
       provided by detinfo::ElecClock start_time input variable which can also used
       to retrieve G4 time offset to electronics clock counting.
     */
    virtual void Process(std::vector<float> &wf,
                         const detinfo::DetectorClocksData& clockData,
			 const ::detinfo::ElecClock &start_time) = 0;

    inline void SetOpChannel(int ch) { fOpChannel = ch; }

    inline int OpChannel() const { return fOpChannel; }

    virtual void Reset() {}

  protected:

    int fOpChannel;
    /**
       A utility function to combine 2 vectors. out_wf is the output vector
       to which in_wf is added. The function adds in_wf from index=start_index
       of out_wf. Function only adds relevant portion of in_wf, which means
       extra part of in_wf that does not fit in out_wf is ommitted.
    */
    void CombineVector(std::vector<float> &out_wf,
		       const size_t start_index,
		       const std::vector<float> &in_wf) const;
    
  };
}

#endif
/** @} */ // end of doxygen group 
