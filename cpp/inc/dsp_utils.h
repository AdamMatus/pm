#ifndef _DSP_UTILS_H_
#define _DSP_UTILS_H_

#include <cstddef>
#include <array>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_sf_trig.h>

namespace dsp_utils{

const constexpr int BUFFER_LEN = 1024;
const constexpr int FFT_SIZE = 256; 
const constexpr int N = FFT_SIZE;
const constexpr int FRAME_STEP = 100; 
const constexpr int M = FRAME_STEP;
const constexpr int MELL_FILTER_BANKS = 30;
const constexpr int K = MELL_FILTER_BANKS;
const constexpr int MFCC_NUM = 13;

  template <int N>
  using Frame = typename std::array<double, N>;

  enum class Window_type
  {
    hann_generator,
    hamming_generator,
  };

  struct window_generator
  {
    window_generator(int N);
    const int size;
    int index;
    double operator()() {return 1.0;}; //rectngle window
  };

  struct hann_generator: public window_generator
  {
    hann_generator(int N);
    double operator()(double x);
  };

  struct hamming_generator: public window_generator
  {
    hamming_generator(int N);
    double operator()(double x);
  };

  struct cos_dct_gen
  {
    cos_dct_gen(int _K): K{_K}, i{0}
    {}

    double operator()(){
      double c = gsl_sf_cos(2*M_PI*static_cast<double>(i)/K);
      ++i;
      return c;
    }
    int K, i;
    
  };

  template <int K>
  struct mfcc_gen
  {
    mfcc_gen(const std::array<double, K>& mel_frame, const std::array<double, 4*K>& cos_dct_table)
      : mel{mel_frame}, cos{cos_dct_table}, n{0} {}
    
    double operator()(){
      double mfcc;
      auto k = 0;
      {// for k=0
        mfcc = (1.0/K)*mel[k];
      }
      for(k=1; k<K; ++k)//exclude mean val
      {
        mfcc+=(2.0/K)*mel[k]*cos[(k*(2*n+1))%(4*K)]; 
      }
      ++n;
      return mfcc; 
    };

    private:
    const std::array<double, K>& mel;
    const std::array<double, 4*K>& cos;
    int n; 
  };

  void window_frame(std::array<double, 256ul>& fr, const Window_type& win_type);
  void power_fft_frame(std::array<double, 256ul>& fr);
  std::array<double, K> dct_frame(const std::array<double, 30ul>& mel_frame, const std::array<double, 30ul*4>& cos_table);
};
#endif

