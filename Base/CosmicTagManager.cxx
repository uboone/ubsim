#ifndef HITCOSMICTAG_COSMICTAGMANAGER_CXX
#define HITCOSMICTAG_COSMICTAGMANAGER_CXX

#include <sstream>
#include <map>
#include <set>
#include "CosmicTagManager.h"

#include "HitCosmicTagException.h"
#include "StartHitFinderFactory.h"
#include "HitOrdererFactory.h"
#include "DqDsCalculatorFactory.h"
#include "DqDsSmootherFactory.h"
#include "LocalLinearityCalculatorFactory.h"
#include "CustomAlgoFactory.h"


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

    auto const start_hit_finder_name      = mgr_cfg.get<std::string>("StartHitFinderAlgo","");
    auto const hit_orderer_name           = mgr_cfg.get<std::string>("HitOrdererAlgo","");
    auto const dqds_calculator_name       = mgr_cfg.get<std::string>("DqDsCalculatorAlgo","");
    auto const dqds_smoother_name         = mgr_cfg.get<std::string>("DqDsSmootherAlgo","");
    auto const linearity_calculator_name  = mgr_cfg.get<std::string>("LocalLinearityCalculatorAlgo","");
    
    std::vector<std::string> custom_algo_v;
    custom_algo_v = mgr_cfg.get<std::vector<std::string> >("CustomAlgo",custom_algo_v);

    if(!start_hit_finder_name.empty()) _alg_start_hit_finder = StartHitFinderFactory::get().create(start_hit_finder_name,start_hit_finder_name);
    if(!hit_orderer_name.empty()) _alg_hit_orderer = HitOrdererFactory::get().create(hit_orderer_name,hit_orderer_name);
    if(!dqds_calculator_name.empty()) _alg_dqds_calculator = DqDsCalculatorFactory::get().create(dqds_calculator_name,dqds_calculator_name);
    if(!dqds_smoother_name.empty()) _alg_dqds_smoother = DqDsSmootherFactory::get().create(dqds_smoother_name,dqds_smoother_name);
    if(!linearity_calculator_name.empty()) _alg_linearity_calculator = LocalLinearityCalculatorFactory::get().create(linearity_calculator_name,linearity_calculator_name);

    for(auto const& name : custom_algo_v)
      if(!name.empty()) AddCustomAlgo(CustomAlgoFactory::get().create(name,name));
    

    // Checks
  
    if (_alg_start_hit_finder) {

      _alg_start_hit_finder->Configure(main_cfg.get<cosmictag::Config_t>(_alg_start_hit_finder->AlgorithmName()));
      
    }
    if (_alg_hit_orderer) {

      _alg_hit_orderer->Configure(main_cfg.get<cosmictag::Config_t>(_alg_hit_orderer->AlgorithmName()));
      
    }
    if (_alg_dqds_calculator) {

      _alg_dqds_calculator->Configure(main_cfg.get<cosmictag::Config_t>(_alg_dqds_calculator->AlgorithmName()));
      
    }
    if (_alg_dqds_smoother) {

      _alg_dqds_smoother->Configure(main_cfg.get<cosmictag::Config_t>(_alg_dqds_smoother->AlgorithmName()));
      
    }
    if (_alg_linearity_calculator) {

      _alg_linearity_calculator->Configure(main_cfg.get<cosmictag::Config_t>(_alg_linearity_calculator->AlgorithmName()));
      
    }

    _configured = true;
  }

  BaseAlgorithm* CosmicTagManager::GetAlgo(cosmictag::AlgoType type)
  {
    if (!_configured)
      CT_WARNING() << "Algorithm may be not configured yet!" << std::endl;

    // Figure out the type of a provided algorithm
    switch (type) {

    // Start Hit Finder
    case kStartHitFinder:
      return _alg_start_hit_finder;

    // Hit Orderer
    case kHitOrderer:
      return _alg_hit_orderer;

    // dq/ds Calculator
    case kDqDsCalculator:
      return _alg_dqds_calculator;

    // dq/ds Smoother
    case kDqDsSmoother:
      return _alg_dqds_smoother;

    // Local Linearity 
    case kLocalLinearity:
      return _alg_linearity_calculator;

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
      //CT_ERROR() << Form("Algorithm name %s not found!",name.c_str()) << std::endl;
      throw HitCosmicTagException();
    }
    return _custom_alg_m[name];
  }

  

  void CosmicTagManager::Emplace(cosmictag::SimpleCluster && obj)
  { _cluster = obj; }

  void CosmicTagManager::SetStartHit(cosmictag::SimpleHit && hit)
  { _start_hit = hit; }

  
  // CORE FUNCTION
  bool CosmicTagManager::Match()
  {
    // Clear some history variables
    
    
    // Create also a result container
    
    if (!_alg_start_hit_finder)
      throw HitCosmicTagException("Start hit finder algorithm is reuqired! (not attached)");
    if (!_alg_hit_orderer)
      throw HitCosmicTagException("Hit orderer algorithm is reuqired! (not attached)");
    if (!_alg_dqds_calculator)
      throw HitCosmicTagException("dq/ds calculator algorithm is reuqired! (not attached)");
    if (!_alg_dqds_smoother)
      throw HitCosmicTagException("dq/ds smoother algorithm is reuqired! (not attached)");
    if (!_alg_linearity_calculator)
      throw HitCosmicTagException("local linearity algorithm is reuqired! (not attached)");

    if (!_configured)
      throw HitCosmicTagException("Have not configured yet!");

    //if(_tpc_object_v.empty() || _flash_v.empty()) return result;

    //
    // Filter stage: for both TPC and Flash
    //

    return false;

  }

  void CosmicTagManager::PrintConfig() {
    
    std::cout << "---- COSMIC TAG MANAGER PRINTING CONFIG     ----" << std::endl
	      << "_name = " << _name << std::endl 
	      << "_alg_start_hit_finder?" << std::endl;
    if (_alg_start_hit_finder)
      std::cout << "\t" << _alg_start_hit_finder->AlgorithmName() << std::endl;
    std::cout << "_alg_hit_orderer?" << std::endl;
    if (_alg_hit_orderer)
      std::cout << "\t" << _alg_hit_orderer->AlgorithmName() << std::endl;
    std::cout << "_alg_dqds_calculator?" << std::endl;
    if (_alg_dqds_calculator)
      std::cout << "\t" << _alg_dqds_calculator->AlgorithmName() << std::endl;
    std::cout << "_alg_dqds_smoother?" << std::endl;
    if (_alg_dqds_smoother)
      std::cout << "\t" << _alg_dqds_smoother->AlgorithmName() << std::endl;
    std::cout << "_alg_linearity_calculator?" << std::endl;
    if (_alg_linearity_calculator)
      std::cout << "\t" << _alg_linearity_calculator->AlgorithmName() << std::endl;
    std::cout << "_custom_alg_m?" << std::endl;
    for (auto& name_ptr : _custom_alg_m)
      std::cout << "\t" << name_ptr.first << std::endl;
    std::cout << "---- END FLASH MATCH MANAGER PRINTING CONFIG ----" << std::endl;
    
  }
}

#endif
