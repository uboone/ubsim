/**
 * \file BaseOrdererAlgo.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class BaseOrdererAlgo
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Base

    @{*/
#ifndef HITCOSMICTAG_BASEORDERERALGO_H
#define HITCOSMICTAG_BASEORDERERALGO_H

#include "BaseAlgorithm.h"

namespace cosmictag {
  /**
     \class BaseOrdererAlgo
     Algorithm base class for prohibiting the match
     between a charge cluster and a flash \n
  */
  class BaseOrdererAlgo : public BaseAlgorithm{
    
  public:
    
    /// Default constructor
    BaseOrdererAlgo(const std::string name="noname") : BaseAlgorithm(kHitOrderer,name)
    {}
 
    /// Default destructor
    virtual ~BaseOrdererAlgo(){}
    
    /**
     * @brief CORE FUNCTION: order the hits
     */
    virtual bool OrderHits(const SimpleCluster&) = 0;
  };
}

#endif
/** @} */ // end of doxygen group 

