#ifndef HITCOSMICTAG_DATATYPES_H
#define HITCOSMICTAG_DATATYPES_H

#include <vector>
#include <numeric>
#include <string>

namespace cosmictag {

  /// Enumerator for different types of algorithm
  enum AlgoType{
    kAlgoUnknown = -1,  ///< Algorithm type to f?
    kStartHitFinder = 0,
    kHitOrderer = 1,
    kDqDsCalculator = 2,
    kDqDsSmoother = 3,
    kCustomAlgo,
    kAlgorithmTypeMax   ///< enum flag for algorithm type count & invalid type
  };

  /// Struct to represent a hit
   struct SimpleHit {

    public: 

    double time;     ///< in cm, basically x
    int wire;        ///< in cm, given the wire pitch
    int plane;       ///< plane number
    double integral; ///< integral (~charge)

    double t;        ///< in time ticks
    int w;           ///< in wire number

    SimpleHit () {
      time = t = -1;
      wire = w = -1;
      plane = -1;
    }

    bool operator==(const SimpleHit& x) const {
      return (time == x.time) && 
             (wire == x.wire) &&
             (plane == x.plane);
    }
  };

  /// Struct to represent a cluster of SimpleHits
   struct SimpleCluster {

    public: 

    std::vector<SimpleHit> _s_hit_v;
    int _start_index = -1;
    std::vector<double> _dqds_v;
    std::vector<double> _ds_v;
    std::vector<double> _dqds_slider;
    std::vector<double> _linearity_v;

    bool _linearity_is_set = false;
    bool _hits_ordered = false;
    bool _start_hit_is_set = false;

    SimpleCluster () {
      _s_hit_v.clear();
      _dqds_v.clear();
      _ds_v.clear();
      _dqds_slider.clear();
      _linearity_v.clear();

      _linearity_is_set = false;
      _hits_ordered = false;
      _start_hit_is_set = false;
    }
  };

 
  

  
  
  namespace msg {
    /// Verbosity message level
    enum Level_t {
      kDEBUG,
      kINFO,
      kNORMAL,
      kWARNING,
      kERROR,
      kEXCEPTION,
      kMSG_TYPE_MAX
    };
    
    const std::string kStringPrefix[kMSG_TYPE_MAX] =
      {
	"\033[94m     [DEBUG]  \033[00m", ///< DEBUG message prefix
	"\033[92m      [INFO]  \033[00m", ///< INFO message prefix
	"\033[95m    [NORMAL]  \033[00m", ///< NORMAL message prefix
	"\033[93m   [WARNING]  \033[00m", ///< WARNING message prefix
	"\033[91m     [ERROR]  \033[00m", ///< ERROR message prefix
	"\033[5;1;33;41m [EXCEPTION]  \033[00m"  ///< CRITICAL message prefix
      };
    ///< Prefix of message
  }
}
#endif
