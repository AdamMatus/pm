#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <numeric>

#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>

#include <sndfile.h>
#include <cstring>

#include <cassert>


#include "main.h"
#include "dsp_utils.h"
#include "mel_frame_generator.h"
#include "vq.h"


inline double mel(double f);
inline double mel2freq(double f);

std::vector<double> load_wav(const char* file, SF_INFO &info);

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

const constexpr int BUFFER_LEN = 1024;
const constexpr int FFT_SIZE = 256; 
const constexpr int N = FFT_SIZE;
const constexpr int FRAME_STEP = 100; 
const constexpr int M = FRAME_STEP;
const constexpr int MELL_FILTER_BANKS = 30;
const constexpr int K = MELL_FILTER_BANKS;
const constexpr int MFCC_NUM = 13;

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
  //test_mpl(wav_probes.cbegin(), wav_probes.cend());
  //test

  //construct matrix of N-size frames
  constexpr auto test_frame = 50;
  std::vector<std::array<double, N>> speech_frames;
  for(int signal_offset = 0; signal_offset + N < static_cast<int>(wav_probes.size()); signal_offset += M)
  {
    auto frame = std::array<double, N>();
    auto start_pos = wav_probes.cbegin() + signal_offset;
    std::copy(start_pos, start_pos + N, frame.begin());
    speech_frames.push_back(std::move(frame));
  }

  std::vector<std::array<double, K>> mel_coefs_speech_frames;
  for(auto &frame: speech_frames)
  {
    dsp_utils::window_frame(frame, dsp_utils::Window_type::hamming_generator);
    dsp_utils::power_fft_frame(frame);
    mel_coefs_speech_frames.push_back(mel_utils::mel_frame(frame, samplerate));
  }

  //Log lCm (loudness)
  for(auto &mel_frame: mel_coefs_speech_frames)
  {
    std::for_each(mel_frame.begin(), mel_frame.end(), [](double &val){val = std::log10(val); });
  } 

  //test
  //test_mpl(mel_coefs_speech_frames.at(test_frame).cbegin(), mel_coefs_speech_frames.at(test_frame).cend(), 1);
  //test
  
  //DCT - final MFCC
  //cos table
  std::array<double, 4*K> cos_table;
  std::generate(cos_table.begin(), cos_table.end(), cos_dct_gen(4*K));

  //computing dct
  std::vector<std::array<double, MFCC_NUM>> mfcc;
  for(const auto &mel_frame: mel_coefs_speech_frames)
  {
    mfcc.push_back(dsp_utils::dct_frame(mel_frame, cos_table));
  }

  //test
  //test_mpl(mfcc.at(test_frame).cbegin(), mfcc.at(test_frame).cend(), 1);
  //test
  //
  vq::dummy_test();

  return 0;
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

