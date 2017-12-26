/**
 * \file LocalLinearityCalculatorAlgo.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class LocalLinearityCalculatorAlgo
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Base

    @{*/
#ifndef HITCOSMICTAG_LOCALLINEARITYCALCULATORALGO_H
#define HITCOSMICTAG_LOCALLINEARITYCALCULATORALGO_H

#include "BaseAlgorithm.h"

namespace cosmictag {
  /**
     \class LocalLinearityCalculatorAlgo
     Algorithm base class for prohibiting the match
     between a charge cluster and a flash \n
  */
  class LocalLinearityCalculatorAlgo : public BaseAlgorithm{
    
  public:
    
    /// Default constructor
    LocalLinearityCalculatorAlgo(const std::string name="noname") 
       : BaseAlgorithm(kHitOrderer,name)
    {}
 
    /// Default destructor
    virtual ~LocalLinearityCalculatorAlgo(){}
    
    /**
     * @brief CORE FUNCTION: order the hits
     */
    virtual int CalculateLocalLinearity(SimpleCluster&) const = 0;
  };
}

#endif
/** @} */ // end of doxygen group 

