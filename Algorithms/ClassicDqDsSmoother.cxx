#ifndef CLASSICSDQDSSMOOTHER_CXX
#define CLASSICSDQDSSMOOTHER_CXX

#include "ClassicDqDsSmoother.h"

namespace cosmictag {

  static ClassicDqDsSmootherFactory __global_ClassicDqDsSmootherFactory__;

  ClassicDqDsSmoother::ClassicDqDsSmoother(const std::string name)
    : BaseDqDsSmootherAlgo(name)
  {}

  void ClassicDqDsSmoother::_Configure_(const Config_t &pset)
  {
    _slider_window = pset.get<size_t>("SliderWindow");
  }
  
  int ClassicDqDsSmoother::SmoothDqDs(SimpleCluster& cluster) const
  {

    //int                    & _start_index      = cluster._start_index;
    std::vector<SimpleHit> & _s_hit_v          = cluster._s_hit_v;
    std::vector<double>    & _dqds_v           = cluster._dqds_v;
    std::vector<double>    & _dqds_slider      = cluster._dqds_slider;
    //bool                   & _hits_ordered     = cluster._hits_ordered;
    //std::vector<double>    & _ds_v        = cluster._ds_v;

     if (_dqds_v.size() != _s_hit_v.size()) {
      CT_ERROR() << " dqds size is different than hit vector size" << std::endl;
      throw HitCosmicTagException();
    }

    if (_dqds_v.size() < _slider_window * 2) {
      CT_NORMAL() <<  " not enough hits" << std::endl;
      return 0;
    }

    //size_t window = _slider_window;

    if (_slider_window % 2 != 0) {
      CT_ERROR() << ": _slider_window has to be even." << std::endl; 
      throw HitCosmicTagException();
    }

    _dqds_slider.clear();
    _dqds_slider.reserve(_dqds_v.size());


    for(const auto& window : get_windows(_dqds_v, _slider_window)) {

      double median_dqds = get_smooth_trunc_median(window);
      _dqds_slider.push_back(median_dqds);
      //if (_debug) std::cout << "dqds_slider value " << median_dqds << std::endl;


    }

    return _s_hit_v.size();
  }

  template<typename T>
  std::vector<std::vector<T>> get_windows(const std::vector<T>& the_thing,
                                          const size_t window_size)
  {

    // given a vector of values return a vector of the same length
    // with each element being a vector of the values of the local neighbors
    // of the element at position i in the original vector
    // input  : [0,1,2,3,4,5,6,...,...,N-3,N-2,N-1] (input vector of size N)
    // output  (assuming a value of 'w' below == 3):
    // 0th element: [0]
    // 1st element: [0,1,2]
    // 2nd element: [0,1,2,3,4]
    // jth element: [j-w,j-w+1,..,j+w-2,j+w-1]
    
    std::vector<std::vector<T>> data;
    
    auto w = window_size + 2;
    w = (unsigned int)((w - 1)/2);
    auto num = the_thing.size();
    
    data.reserve(num);
    
    for(size_t i = 1; i <= num; ++i) {
      std::vector<T> inner;
      inner.reserve(20);
      // if we are at the beginning of the vector (and risk accessing -1 elements)
      if(i < w)
        {
          for(size_t j = 0; j < 2 * (i%w) - 1; ++j)
            inner.push_back(the_thing[j]);
        }
      // if we are at the end of the vector (and risk going past it)
      else if (i > num - w + 1)
        {
          for(size_t j = num - 2*((num - i)%w)-1 ; j < num; ++j)
            inner.push_back(the_thing[j]);
        }
      // if we are in the middle of the waveform
      else
        {
          for(size_t j = i - w; j < i + w - 1; ++j)
            inner.push_back(the_thing[j]);
        }
      data.emplace_back(inner);
    }

    return data;
  
  }
}
#endif
