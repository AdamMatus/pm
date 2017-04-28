#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <numeric>

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>

#include <sndfile.h>
#include <cstring>

#include <cassert>

#include "matplotlibcpp.h"

namespace mpl = matplotlibcpp;

inline double mel(double f);
inline double mel2freq(double f);

template <typename It>
void test_mpl(It start, It end, size_t start_index = 0);

std::vector<double> load_wav(const char* file, SF_INFO &info);

struct hann_generator
{
  hann_generator(size_t _N): size{_N}, index{0} {}
  double operator()() {
      return 0.5*(1.0 - gsl_sf_cos(2.0*M_PI*(index++)/(size -1)));
    }
  const size_t size;
  size_t index;
};

struct hamming_generator
{
  hamming_generator(size_t _N): size{_N}, index{0} {}
  double operator()() {
      return 0.54 - 0.46*gsl_sf_cos(2.0*M_PI*(index++)/(size -1));
    }
  const size_t size;
  size_t index;
};

struct triangle_windowed_sum
{
  triangle_windowed_sum(size_t f_center, size_t f_end)
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

  const size_t cen, end;
  size_t i;
  double acc, win, del;
};

template <size_t K, size_t N>
struct mel_frame_generator
{
  typedef typename std::array<double, K> arrK;
  typedef typename std::array<double, N> arrN;

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
  size_t f_begin, f_center, f_end;

  inline void update_filter_samples(){
    f_begin = sample(0);
    f_center = sample(1);
    f_end = sample(2); 

    k++;
  }

  inline size_t sample(int inc){
    return static_cast<size_t>(((mel2freq((k+inc)*delta_mel))/fs)*N);
  }
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

constexpr auto BUFFER_LEN = 1024;
constexpr auto FFT_SIZE = 256; 
constexpr auto N = FFT_SIZE;
constexpr auto FRAME_STEP = 100; 
constexpr auto M = FRAME_STEP;
constexpr auto MELL_FILTER_BANKS = 30;
constexpr auto K = MELL_FILTER_BANKS;
constexpr auto MFCC_NUM = 13;

int main (void)
{
  //load single test wav
  const char *file = "s1.wav";
  SF_INFO wavinfo;
  auto wav_probes = load_wav(file, wavinfo);
  sf_format_check(&wavinfo); 
  const auto samplerate = wavinfo.samplerate;
  const auto delta_samplerate = samplerate/FFT_SIZE;

  //test
  test_mpl(wav_probes.cbegin(), wav_probes.cend());
  //test

  //construct matrix of N-size frames
  constexpr auto test_frame = 50;
  std::vector<std::array<double, N>> speech_frames;
  for(size_t signal_offset = 0; signal_offset + N < wav_probes.size(); signal_offset += M)
  {
    auto frame = std::array<double, N>();
    auto start_pos = wav_probes.cbegin() + signal_offset;
    std::copy(start_pos, start_pos + N, frame.begin());
    speech_frames.push_back(std::move(frame));
  }

  //applying window to frames
  auto hamming_window= std::array<double, N>();
  std::generate(hamming_window.begin(), hamming_window.end(), hann_generator(N));
  for(auto &frame: speech_frames)
  {
    std::transform(hamming_window.cbegin(), hamming_window.cend(),
                   frame.cbegin(), frame.begin(), std::multiplies<double>());
  }

  //abs(fft)^2 on matrix of frames
  //
  //magic '+1' comes from definition of gsl_fft_real_radix2_transform:
  //frame[0] - fft(t=0) and frame[N/2] - fft(t=N/2) have no stored representation for imgainary part
  for(auto &frame: speech_frames)
  {
    gsl_fft_real_radix2_transform(frame.data() ,1 ,N);
    std::for_each(frame.begin(), frame.end(), [](double &x) {x = gsl_pow_2(x); });
    std::transform(frame.cbegin() +1, frame.cbegin() + (N/2),
                   frame.crbegin(),
                   frame.begin() +1,
                   std::plus<double>());
  }

  //test
  //test_mpl(speech_frames.at(test_frame).cbegin(), speech_frames.at(test_frame).cend(), 1);
  //test
  
  //Mel-frequency Wrapping - triangle filters bank
  std::vector<std::array<double, K>> mel_coefs_speech_frames;
  for(const auto &speech_frame: speech_frames)
  {
    std::array<double, K> mel_frame;
    std::generate(mel_frame.begin(), mel_frame.end(), mel_frame_generator<K, N>(speech_frame, samplerate));
    mel_coefs_speech_frames.push_back(mel_frame);
  }

  //test
  //test_mpl(mel_coefs_speech_frames.at(test_frame).cbegin(), mel_coefs_speech_frames.at(test_frame).cend(), 1);
  //test

  //Log lCm (loudness)
  for(auto &mel_frame: mel_coefs_speech_frames)
  {
    std::for_each(mel_frame.begin(), mel_frame.end(), [](double &val){
          val = std::log10(val);
        });
  } 

  //test
  test_mpl(mel_coefs_speech_frames.at(test_frame).cbegin(), mel_coefs_speech_frames.at(test_frame).cend(), 1);
  //test
  
  //DCT - final MFCC
  //cos table
  std::array<double, 4*K> cos_table;
  std::generate(cos_table.begin(), cos_table.end(), cos_dct_gen(4*K));

  //test
  test_mpl(cos_table.begin(), cos_table.begin() + K, 0); 
  //test
  
  //computing dct
  std::vector<std::array<double, MFCC_NUM>> mfcc;
  for(const auto &mel_frame: mel_coefs_speech_frames)
  {
    std::array<double, MFCC_NUM> mfcc_frame;
    std::generate(mfcc_frame.begin(), mfcc_frame.end(), mfcc_gen<K>(mel_frame, cos_table));
    mfcc.push_back(std::move(mfcc_frame));
  }

  //test
  /* 
  std::array<double, K> test;
  std::fill(test.begin(), test.end(), 0);
  test.at(1) = 1;
  std::array<double, MFCC_NUM> ot;
  std::generate(ot.begin(), ot.end(), mfcc_gen<K>(test, cos_table));
  test_mpl(ot.begin(), ot.end(), 1);

  std::array<double, 4*2> c;
  std::generate(c.begin(), c.end(), cos_dct_gen(4*2));
  std::array<double, 2> t;
  std::fill(t.begin(), t.end(), 1);
  std::array<double, 2> o;
  std::generate(o.begin(), o.end(), mfcc_gen<2>(t, c));
  t;
  */
  //
  //test

  //test
  test_mpl(mfcc.at(test_frame).cbegin(), mfcc.at(test_frame).cend(), 1);
  //test

  return 0;
}

inline double mel(double f)
{
  return 2595*std::log10(1+f/700);
}

inline double mel2freq(double mel)
{
  return 700*(std::pow(10, mel/2595)-1);
}

template <typename It>
void test_mpl(It start, It end, size_t start_index)
{
  typedef typename std::iterator_traits<It>::value_type T;

  auto t = std::vector<T>(end-start);
  std::iota(t.begin(), t.end(), start_index);

  auto y = std::vector<T>(start, end);

  mpl::plot(t, y);
  mpl::show();
}

std::vector<double> load_wav(const char* file, SF_INFO&  sfinfo)
{
  double data[BUFFER_LEN];

  SNDFILE *infile;
  memset(&sfinfo, 0, sizeof(sfinfo));

  if(!(infile = sf_open(file, SFM_READ, &sfinfo)))
  {
    std::cerr << "Not able to load a file: " << file << std::endl;
    std::cerr << sf_strerror(NULL);
    return std::vector<double>();
  }

  int readcount;
  std::vector<double> wav_probes;
  while((readcount = sf_read_double(infile, data, BUFFER_LEN)))
  {
    wav_probes.insert(wav_probes.end(), data, data + readcount);      
  }

  return wav_probes;
}

