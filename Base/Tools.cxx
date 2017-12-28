#ifndef HITCOSMICTAG_TOOLS_CXX
#define HITCOSMICTAG_TOOLS_CXX

#include "Tools.h"

#include <vector>
#include <numeric>
#include <string>

namespace cosmictag {

  /// Enumerator for different types of algorithm
  std::vector<std::vector<double>> get_windows(const std::vector<double>& the_thing,
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
    
    std::vector<std::vector<double>> data;
    
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







    double get_smooth_trunc_median(std::vector<double> v) {

    if (v.size() > 2) {
      // Find and erase max element
      auto it_max = std::max_element(v.begin(), v.end());
      v.erase(it_max);

      // Find and erase min element
      auto it_min = std::min_element(v.begin(), v.end());
      v.erase(it_min);
    }

    double median = -1;

    size_t size = v.size();
    std::sort(v.begin(), v.end());
    if (size % 2 == 0){
      median = (v[size/2 - 1] + v[size/2]) / 2;
    }
    else{
      median = v[size/2];
    }

    return median;
  }



  double mean(const std::vector<double>& data)
  {
    if(data.size() == 0) std::cout << __PRETTY_FUNCTION__ << "You have me nill to mean" << std::endl;
        
    double result = 0.0;

    for(const auto& d : data) 
      result += d;
        
    return (result / ((double)data.size()));
  }





  double cov (const std::vector<double>& data1,
              const std::vector<double>& data2)
  {
    if(data1.size() == 0) std::cout << __PRETTY_FUNCTION__ << "You have me nill to cov" << std::endl;
    if(data2.size() == 0) std::cout << __PRETTY_FUNCTION__ << "You have me nill to cov" << std::endl;

    double result = 0.0;
    auto   mean1  = mean(data1);
    auto   mean2  = mean(data2);
    
    for(size_t i = 0; i < data1.size(); ++i)
      result += (data1[i] - mean1)*(data2[i] - mean2);
    
    return result/((double)data1.size());
      
  }




  double stdev(const std::vector<double>& data)
  {
    if(data.size() == 0) std::cout << __PRETTY_FUNCTION__ << "You have me nill to stdev" << std::endl;

    double result = 0.0;
    auto    avg   = mean(data);
    for(const auto& d: data)
      result += (d - avg)*(d - avg);
    
    return std::sqrt(result/((double)data.size()));
  }




  

}
#endif
