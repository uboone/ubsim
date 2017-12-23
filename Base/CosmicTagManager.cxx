#ifndef HITCOSMICTAG_COSMICTAGMANAGER_CXX
#define HITCOSMICTAG_COSMICTAGMANAGER_CXX

#include <sstream>
#include <map>
#include <set>
#include "CosmicTagManager.h"

#include "HitCosmicTagException.h"
#include "HitOrdererFactory.h"
#include "CustomAlgoFactory.h"

//#include "TPCFilterFactory.h"
//#include "FlashMatchFactory.h"
//#include "FlashHypothesisFactory.h"
//#include "FlashProhibitFactory.h"


namespace cosmictag {

  CosmicTagManager::CosmicTagManager(const std::string name)
    : _alg_hit_orderer(nullptr)
    , _name(name)
  {
    _configured = false;
  }

  const std::string& CosmicTagManager::Name() const
  { return _name; }


  void CosmicTagManager::AddCustomAlgo(BaseAlgorithm* alg)
  {
    if(_custom_alg_m.find(alg->AlgorithmName()) != _custom_alg_m.end()) {
      std::stringstream ss;
      ss << "Duplicate name: " << alg->AlgorithmName() << std::endl;
      throw HitCosmicTagException(ss.str());
      }
    _custom_alg_m[alg->AlgorithmName()] = alg;
  }
  
  void CosmicTagManager::Configure(const Config_t& main_cfg) 
  {
   
    auto const& mgr_cfg = main_cfg.get<cosmictag::Config_t>(Name());

    this->set_verbosity((msg::Level_t)(mgr_cfg.get<unsigned int>("Verbosity")));

    auto const hit_orderer_name = mgr_cfg.get<std::string>("HitOrdererAlgo","");
    
    std::vector<std::string> custom_algo_v;
    custom_algo_v = mgr_cfg.get<std::vector<std::string> >("CustomAlgo",custom_algo_v);

    if(!hit_orderer_name.empty()) _alg_hit_orderer = HitOrdererFactory::get().create(hit_orderer_name,hit_orderer_name);
    
    for(auto const& name : custom_algo_v)
      if(!name.empty()) AddCustomAlgo(CustomAlgoFactory::get().create(name,name));
    
    // checks
  
    if (_alg_hit_orderer) {

      _alg_hit_orderer->Configure(main_cfg.get<cosmictag::Config_t>(_alg_hit_orderer->AlgorithmName()));
      
    }

    _configured = true;
  }

  BaseAlgorithm* CosmicTagManager::GetAlgo(cosmictag::AlgoType type)
  {
    if (!_configured)
      CT_WARNING() << "Algorithm may be not configured yet!" << std::endl;

    // Figure out the type of a provided algorithm
    switch (type) {

    // TPC filter
    case kHitOrderer:
      return _alg_hit_orderer;
/*
    // Flash filter
    case kFlashFilter:
      return _alg_flash_filter;

    // Match prohibit algo
    case kMatchProhibit:
      return _alg_match_prohibit;

    // Flash matching
    case kFlashMatch:
      return _alg_flash_match;

    // Flash hypothesis
    case kFlashHypothesis:
      return _alg_flash_hypothesis;
*/
    // Fuck it
    default:
      std::stringstream ss;
      ss << "Unsupported algorithm type: " << type;
      throw HitCosmicTagException(ss.str());
    }
    return nullptr;
  }

  cosmictag::BaseAlgorithm* CosmicTagManager::GetCustomAlgo(std::string name)
  {
    if(_custom_alg_m.find(name) == _custom_alg_m.end()) {
      CT_ERROR() << Form("Algorithm name %s not found!",name.c_str()) << std::endl;
      throw HitCosmicTagException();
    }
    return _custom_alg_m[name];
  }

  

  void CosmicTagManager::Emplace(cosmictag::SimpleCluster && obj)
  { _cluster = obj; }

  
  // CORE FUNCTION
  bool CosmicTagManager::Match()
  {
    // Clear some history variables
    
    
    // Create also a result container
    
    if (!_alg_hit_orderer)
      throw HitCosmicTagException("Hit orderer algorithm is reuqired! (not attached)");

    if (!_configured)
      throw HitCosmicTagException("Have not configured yet!");

    //if(_tpc_object_v.empty() || _flash_v.empty()) return result;

    //
    // Filter stage: for both TPC and Flash
    //

    return false;

  }

  void CosmicTagManager::PrintConfig() {
    
    std::cout << "---- FLASH MATCH MANAGER PRINTING CONFIG     ----" << std::endl
	      << "_name = " << _name << std::endl;
        /*
	      << "_alg_flash_filter?" << std::endl;
    if (_alg_flash_filter)
      std::cout << "\t" << _alg_flash_filter->AlgorithmName() << std::endl;
    std::cout << "_alg_tpc_filter?" << std::endl;
    if (_alg_tpc_filter)
      std::cout << "\t" << _alg_tpc_filter->AlgorithmName() << std::endl;
    std::cout << "_alg_match_prohibit?" << std::endl;
    if (_alg_match_prohibit)
      std::cout << "\t" << _alg_match_prohibit->AlgorithmName() << std::endl;
    std::cout << "_alg_flash_hypothesis?" << std::endl;
    if (_alg_flash_hypothesis)
      std::cout << "\t" << _alg_flash_hypothesis->AlgorithmName() << std::endl;
    std::cout << "_alg_flash_match?" << std::endl;
    if (_alg_flash_match)
      std::cout << "\t" << _alg_flash_match->AlgorithmName() << std::endl;
    std::cout << "_custom_alg_m?" << std::endl;
    for (auto& name_ptr : _custom_alg_m)
      std::cout << "\t" << name_ptr.first << std::endl;
    std::cout << "---- END FLASH MATCH MANAGER PRINTING CONFIG ----" << std::endl;
    */
  }
}

#endif
