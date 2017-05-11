#ifndef _MEL_FRAME_GENERATOR_
#define _MEL_FRAME_GENERATOR_

#include <array>
#include <algorithm>

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

};

#endif
