#ifndef _MEL_FRAME_GENERATOR_
#define _MEL_FRAME_GENERATOR_

#include <array>
#include <algorithm>

#include "dsp_utils.h"

namespace mel_utils
{

  inline double mel(double f);
  inline double mel2freq(double mel);

  struct triangle_windowed_sum
  {
    triangle_windowed_sum(int f_center, int f_end)
      : cen{f_center}, end{f_end}, i{0}, acc{0}, win{0}
    {
      del = 1/static_cast<float>(cen);
    }

    void operator()(const double& val){
      if(cen == i) del = -1.0/static_cast<float>((end-cen));

      acc += win*val;
      win += del;
      i++;
    }

    const int cen, end;
    int i;
    double acc, win, del;
  };

  template <int K, int N>
  struct mel_frame_generator
  {
    typedef std::array<double,K> arrK;
    typedef std::array<double,N> arrN;

    mel_frame_generator(const arrN& frame, double samplerate)
      : fs{samplerate}, mel_max{mel(fs/2)}, delta_mel{mel_max/(K+1)},
        fram{frame},
        k{0}
    {
      update_filter_samples();
    }

    double operator()(){
      auto mel_coef = std::for_each(fram.cbegin() + f_begin,
                      fram.cbegin() + f_end,
                      triangle_windowed_sum(f_center-f_begin, f_end-f_begin));
      update_filter_samples();
      return mel_coef.acc;
    };

  private:
    const double fs, mel_max, delta_mel;
    const arrN &fram;
    int k;
    int f_begin, f_center, f_end;

    inline void update_filter_samples(){
      f_begin = sample(0);
      f_center = sample(1);
      f_end = sample(2); 

      k++;
    }

    inline int sample(int inc){
      return static_cast<int>(((mel2freq((k+inc)*delta_mel))/fs)*N);
    }
  };

  std::array<double, 30> mel_frame(const std::array<double, 256ul>& fr, int samplerate);

  template <int N, int K>
  std::vector<std::array<double, K>> mfcc_extraction(std::vector<std::array<double, N>>&& speech_frames, int samplerate)
  {
    std::vector<std::array<double, K>> mel_coefs_speech_frames;
    for(auto &frame: speech_frames)
    {
      dsp_utils::window_frame(frame, dsp_utils::Window_type::hamming_generator);
      dsp_utils::power_fft_frame(frame);
      mel_coefs_speech_frames.push_back(mel_frame(frame, samplerate));
    }

    for(auto &mel_frame: mel_coefs_speech_frames)
    {
      std::for_each(mel_frame.begin(), mel_frame.end(), [](double &val){val = std::log10(val); });
    } 
    
    std::array<double, 4*K> cos_table;
    std::generate(cos_table.begin(), cos_table.end(), dsp_utils::cos_dct_gen(4*K));

    std::vector<std::array<double, K>> mfcc;
    for(const auto &mel_frame: mel_coefs_speech_frames)
    {
      mfcc.push_back(dsp_utils::dct_frame(mel_frame, cos_table));
    }
    return mfcc;
  }
};

#endif
