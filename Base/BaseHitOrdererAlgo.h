/**
 * \file BaseHitOrdererAlgo.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class BaseOrdererAlgo
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Base

    @{*/
#ifndef HITCOSMICTAG_BASEHITORDERERALGO_H
#define HITCOSMICTAG_BASEHITORDERERALGO_H

#include "BaseAlgorithm.h"

namespace cosmictag {
  /**
     \class BaseOrdererAlgo
     Algorithm base class for prohibiting the match
     between a charge cluster and a flash \n
  */
  class BaseHitOrdererAlgo : public cosmictag::BaseAlgorithm{
    
  public:
    
    /// Default constructor
    BaseHitOrdererAlgo(const std::string name="noname") 
       : cosmictag::BaseAlgorithm(cosmictag::kHitOrderer,name)
    {}
 
    /// Default destructor
    virtual ~BaseHitOrdererAlgo(){}
    
    /**
     * @brief CORE FUNCTION: order the hits
     */
    virtual bool OrderHits(const SimpleCluster&) const = 0;
  };
}

#endif
/** @} */ // end of doxygen group 

