#ifndef _VQ_H_
#define _VQ_H_

#include <array>
#include <vector>

namespace vq
{
  template <int C, int K>
  using codebook_array = typename std::array<std::array<double, K>, C>; 

  void dummy_test();

  template<int C, int K>
  std::array<std::array<double, K>, C> lbg(const std::vector<std::array<double, K>>& acoustic_vectors);  

  template<int K>
  double dis_eu(const std::array<double, K>& v1, const std::array<double, K>& v2);
};

#endif
