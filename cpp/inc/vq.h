#ifndef _VQ_H_
#define _VQ_H_

#include <array>
#include <vector>
#include <algorithm>
#include <cmath>

namespace vq
{
  template <int C, int K>
  using codebook_array = typename std::array<std::array<double, K>, C>; 

  void dummy_test();

  template<int K>
  double dis_eu(const std::array<double, K>& v1, const std::array<double, K>& v2)
  {
    double dist =0;

    auto i1 = v1.cbegin();
    auto i2 = v2.cbegin();
    for(;i1 != v1.cend(); ++i1, ++i2)
    {
      dist += std::pow(*i1-*i2, 2); 
    } 

    return std::sqrt(dist);
  }

  template<int C, int K>
  std::array<std::array<double, K>, C> lbg(const std::vector<std::array<double, K>>& acoustic_vectors)  
  {
    codebook_array<C, K> code;
    std::fill(code.at(0).begin(), code.at(0).end(), 0);
    for(auto&& a_vec: acoustic_vectors)
    {
      std::transform(a_vec.cbegin(), a_vec.cend(), code.at(0).cbegin(), code.at(0).begin(),
                     [](double x, double acc){ return x+acc;}); 
    }
    for(auto&& x: code.at(0))
    {
      x /= acoustic_vectors.size();
    }

    const double eps = 0.01;

    std::vector<int> owner_ind(acoustic_vectors.size(), 0);
    std::vector<double> owner_dist(acoustic_vectors.size(), 0);
    std::vector<int> owner_ind_count(C, 0);
    std::vector<double> owner_distortion(C, std::numeric_limits<double>::max());
    for(auto code_size= 1; code_size < C;)
    {
      //double centroids
      for(auto j = 0; j < code_size; ++j)
      {
        std::transform(code.at(j).cbegin(), code.at(j).cend(), code.at(j+code_size).begin(),
            [eps](double x){
             return x*(1.0+eps); 
            });
        for(auto && x: code.at(j))
        {
          x *= 1.0-eps;
        }
      }
      code_size *= 2;

      auto distortion_diff_small_enough = false;
      while(!distortion_diff_small_enough)
      {
        std::fill(owner_ind_count.begin(), owner_ind_count.end(), 0);
        //find owner indexes
        auto a_vec_index =0;
        for(auto a_vec_iter = acoustic_vectors.cbegin(); a_vec_iter != acoustic_vectors.cend(); ++a_vec_iter, ++a_vec_index)
        {
          std::vector<double> distances;
          for(auto iter = code.cbegin(); iter != code.cbegin() +code_size; ++iter)
          {
            distances.push_back(vq::dis_eu<K>(*a_vec_iter, *iter));
          } 
          auto min_dist = std::min_element(distances.begin(), distances.end());
          auto cen_index = std::distance(std::begin(distances), min_dist);
          owner_ind.at(a_vec_index) = cen_index;
          owner_dist.at(a_vec_index) = distances.at(cen_index);
          ++owner_ind_count.at(cen_index);
        }

        //actualize centroids
        for(auto cen_iter = code.begin(); cen_iter !=  code.begin() +code_size; ++cen_iter)
        {
          std::fill(cen_iter->begin(), cen_iter->end(), 0);
        }
        a_vec_index =0;
        for(auto a_vec_iter = acoustic_vectors.cbegin(); a_vec_iter != acoustic_vectors.cend(); ++a_vec_iter, ++a_vec_index)
        {
          auto owner_centroid_begin = code.at(owner_ind.at(a_vec_index)).begin();
          std::transform(a_vec_iter->cbegin(), a_vec_iter->cend(), owner_centroid_begin, owner_centroid_begin,
                        [](double x, double y){
                          return x+y; 
                        });
        }
        auto cen_index = 0;
        for(auto cen_iter = code.begin(); cen_iter !=  code.begin() +code_size; ++cen_iter, ++cen_index)
        {
          for(auto && c: *cen_iter)
          {
            c /= owner_ind_count.at(cen_index);
          }
        }

        //compute distortion
        std::vector<double> distortion(owner_distortion.size(),0);
        cen_index = 0;
        for(auto i=0ul; i < acoustic_vectors.size(); ++i)
        {
          auto owner_index = owner_ind.at(i);
          distortion.at(owner_index) += owner_dist.at(i);
        }
        for(auto i=0; i < code_size; ++i)
        {
          distortion_diff_small_enough = true;
          auto d_coef = (owner_distortion.at(i)-distortion.at(i))/distortion.at(i); 
          if(d_coef > eps)
          {
            distortion_diff_small_enough = false;
            owner_distortion = std::move(distortion);
            break;
          }
        }
      }
    }

    return code;
  }

};

#endif
