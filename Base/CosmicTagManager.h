/**
 * \file CosmicTagManager.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class CosmicTagManager
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Base

    @{*/
#ifndef HITCOSMICTAG_COSMICTAGMANAGER_H
#define HITCOSMICTAG_COSMICTAGMANAGER_H

#include "LoggerFeature.h"
#include "HitCosmicTagFMWKInterface.h"
#include "BaseAlgorithm.h"
#include "BaseHitOrderer.h"

namespace cosmictag {
  /**
     \class FlashMatchManager
  */
  class CosmicTagManager : public LoggerFeature {

  public:
    
    /// Default constructor
    CosmicTagManager(const std::string name="CosmicTagManager");
    
    /// Default destructor
    ~CosmicTagManager(){}

    /// Name getter
    const std::string& Name() const;

    /// Configuration
    void Configure(const Config_t& cfg);

    /// Algorithm getter
    cosmictag::BaseAlgorithm* GetAlgo(cosmictag::AlgoType type);

    /// Custom algorithm getter
    cosmictag::BaseAlgorithm* GetCustomAlgo(std::string name);
		 
    /// Emplacer of a Cluster
    void Emplace(cosmictag::SimpleCluster && obj);

    /**
       CORE FUNCTION: executes algorithms to find a match of TPC object and flash provided by users. \n
       The execution takes following steps:             \n
       0) TPC filter algorithm if provided (optional)   \n
       1) Flash filter algorithm if provided (optional) \n
       3) Flash matching algorithm (required)           \n
       4) Returns match information for created TPC object & flash pair which respects the outcome of 3)
     */
    bool Match();

    /// Clears locally kept TPC object (QClusterArray_t) and flash (FlashArray_t), both provided by a user
    void Reset()
    { _tpc_object_v.clear(); _flash_v.clear(); }

    /// Configuration option: true => allows an assignment of the same flash to multiple TPC objects
    void CanReuseFlash(bool ok=true)
    { _allow_reuse_flash = ok; }

    void PrintConfig();


  private:

    void AddCustomAlgo(BaseAlgorithm* alg);

    BaseHitOrderer*     _alg_hit_orderer;     ///< Order Hits Algorithm
   
    /**
       A set of custom algorithms (not to be executed but to be configured)
    */
    std::map<std::string,flashana::BaseAlgorithm*> _custom_alg_m;

    /// TPC object information collection (provided by a user)
    SimpleCluster _cluster;


    /// Configuration readiness flag
    bool _configured;
    /// Configuration file
    std::string _config_file;
    /// Name
    std::string _name;
    /// Request boolean to store full matching result (per Match function call)
    bool _store_full;   
  };
}

#endif
/** @} */ // end of doxygen group 

