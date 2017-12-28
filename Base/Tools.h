#ifndef HITCOSMICTAG_TOOLS_H
#define HITCOSMICTAG_TOOLS_H

#include <iostream>
#include <vector>
#include <numeric>
#include <string>
#include <algorithm>

namespace cosmictag {

  /// ?
  std::vector<std::vector<double>> get_windows(const std::vector<double>& the_thing,
                                               const size_t window_size);
  

  /// ?
  double get_smooth_trunc_median(std::vector<double> v);

  
  /// ?
  double mean(const std::vector<double>& data);
  

  /// ?
  double cov (const std::vector<double>& data1,
              const std::vector<double>& data2);
  

  /// ?
  double stdev(const std::vector<double>& data);
  
}
#endif
