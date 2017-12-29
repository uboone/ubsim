#ifndef CLASSICSTARTHITFINDER_CXX
#define CLASSICSTARTHITFINDER_CXX

#include "ClassicStartHitFinder.h"


namespace cosmictag {

  static ClassicStartHitFinderFactory __global_ClassicStartHitFinderFactory__;

  ClassicStartHitFinder::ClassicStartHitFinder(const std::string name)
    : BaseStartHitFinderAlgo(name)
  {}

  void ClassicStartHitFinder::_Configure_(const Config_t &pset)
  {
    _max_allowed_hit_distance = pset.get<double>("MaxAllowedHitDistance");
  }
  
  int ClassicStartHitFinder::FindStartHit(SimpleCluster& cluster, SimpleHit& start_hit) const
  {

    int                    & _start_index      = cluster._start_index;
    std::vector<SimpleHit> & _s_hit_v          = cluster._s_hit_v;
    bool                   & _start_hit_is_set = cluster._start_hit_is_set;
    //std::vector<double>    & _ds_v        = cluster._ds_v;

    _start_hit_is_set = false;

    //if (plane_no != 2) {
    //  std::cout << "Plane not supported." << std::endl;
    //  return;
    //}

    int wire_no = start_hit.wire;
    int time = start_hit.time;

    TVector3 pt1(time, wire_no, 0);
    
    double min_dist = 1e9;
    int best_hit_id = -1;

    if (_debug) std::cout << "Simple hit vector size " << _s_hit_v.size() << std::endl;

    for (size_t i = 0; i < _s_hit_v.size(); i++) {

      auto sh = _s_hit_v.at(i);
      TVector3 pt2(sh.time, sh.wire, 0);
      double dist = (pt1-pt2).Mag();

      if (dist< min_dist) {
        best_hit_id = i;
        min_dist = dist;
      }
    }

    if (best_hit_id == -1) {
      std::cout << "Could not find start hit." << std::endl;
      return 0;
    }
   
 
    // The best hit may not be the start one, but one very close,
    // So let's find a border hit

    auto almost_best_hit = _s_hit_v.at(best_hit_id);

    CT_DEBUG() << "Almost best hit has wire " << almost_best_hit.wire
               << ", and time " << almost_best_hit.time*4 << std::endl;

    TVector3 pt0(almost_best_hit.time, almost_best_hit.wire, 0);

    // First create a map from wire to hit
    std::map<int,SimpleHit> wire_to_hit;

    for (size_t i = 0; i < _s_hit_v.size(); i++) {

      auto sh = _s_hit_v.at(i); 

      // If we never encountered this wire, just save it
      auto iter = wire_to_hit.find(sh.wire);
      if (iter == wire_to_hit.end()) {
        wire_to_hit[sh.wire] = sh;
        continue;
      }

      // Otherwise pick the one that is closer to the almost_best_hit
      auto previous_hit = iter->second;
      TVector3 pt1(previous_hit.time, previous_hit.wire, 0);
      double dist_previous = (pt0-pt1).Mag();
      TVector3 pt2(sh.time, sh.wire, 0);
      double dist_current = (pt0-pt2).Mag();

      if (dist_current < dist_previous)
        wire_to_hit[sh.wire] = sh;
      else 
        continue;
    }

    // Then try going to the left first
    int n_step_left = 0;
    int best_index_left = -1;
    for (int w = almost_best_hit.wire - 1; w > 0; w--) {

      auto iter = wire_to_hit.find(w);
      if (iter == wire_to_hit.end()) {
        break;
      }

      //SimpleHit temp = iter->second;
      auto it = std::find(_s_hit_v.begin(), _s_hit_v.end(), iter->second);
      best_index_left = it - _s_hit_v.begin();//w;//std::distance(_s_hit_v.begin(), it);
      TVector3 pt1(iter->second.time, iter->second.wire, 0);
      double dist = (pt0-pt1).Mag();
  
      if (dist > _max_allowed_hit_distance)  
        break;

      n_step_left ++;

      CT_DEBUG() << "Found hit on the left, wire " << iter->second.wire
                 << ", time " << iter->second.time << std::endl;
    }

    // Then go rigth
    int n_step_right = 0;
    int best_index_right = -1;
    for (int w = almost_best_hit.wire + 1; w < 3456; w++) {

      auto iter = wire_to_hit.find(w);
      if (iter == wire_to_hit.end()) {
        break;
      }

      auto it = std::find(_s_hit_v.begin(), _s_hit_v.end(), iter->second);
      best_index_right = it - _s_hit_v.begin();//w;//std::distance(_s_hit_v.begin(), it);
      TVector3 pt1(iter->second.time, iter->second.wire, 0);
      double dist = (pt0-pt1).Mag(); 

      if (dist > _max_allowed_hit_distance)   
        break; 

      n_step_right ++;

      CT_DEBUG() << "Found hit on the right, wire " << iter->second.wire 
                 << ", time " << iter->second.time*4 << std::endl; 
    }

    _start_index = 0;

    if (n_step_left == 0 || n_step_right == 0) {
      _start_index = best_hit_id;
    } else if (n_step_left > n_step_right) {
      _start_index = best_index_right;
    } else if (n_step_right >= n_step_left) {
      _start_index = best_index_left;
    }


    CT_DEBUG() << "[StoppingMuonTaggerHelper] Start hit set to w: " 
               << _s_hit_v.at(_start_index).wire << ", and t: " 
               << _s_hit_v.at(_start_index).time*4 << std::endl;

    _start_hit_is_set = true;


    return _s_hit_v.size();
  }
}
#endif
