#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <numeric>
#include <string>

#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>

#include <sndfile.h>
#include <cstring>

#include <cassert>

#include "main.h"
#include "dsp_utils.h"
#include "mel_frame_generator.h"
#include "vq.h"
#include "speaker.h"

#include "pm_test.h"

std::vector<double> load_wav(const char* file, SF_INFO &info);

const constexpr int BUFFER_LEN = 1024;
const constexpr int FFT_SIZE = 256; 
const constexpr int N = FFT_SIZE;
const constexpr int FRAME_STEP = 100; 
const constexpr int M = FRAME_STEP;
const constexpr int MELL_FILTER_BANKS = 30;
const constexpr int K = MELL_FILTER_BANKS;
const constexpr int MFCC_NUM = 13;

constexpr const char* wav_file {"train/sx.wav"};
std::string s_wav_files_gen(int i)
{
  std::string file(wav_file);
  if(i<10)
    file.at(7) = '0' + i;
  else if(i<20){
    file.at(7) = '0' + i/10;
    file.insert(file.begin()+8,'0'+i%10);
  }
  return file;
}

const constexpr int S_FILES_NUM = 16;
int main (int argc, char *argv[])
{
  std::vector<speaker<16, K>> speakers; 
  for(int actual_speaker=0; actual_speaker<S_FILES_NUM; ++actual_speaker)
  {
    //load single test wav
    SF_INFO wavinfo;
    auto wav_probes = load_wav(s_wav_files_gen(actual_speaker+1).c_str(), wavinfo);
    std::cout << "training " <<s_wav_files_gen(actual_speaker+1).c_str() << "...\n";
    sf_format_check(&wavinfo); 
    const auto samplerate = wavinfo.samplerate;

    auto speech_frames {dsp_utils::frame_signal<N>(wav_probes, M)};
    auto mfcc {mel_utils::mfcc_extraction<256, 30>(std::move(speech_frames), samplerate)};
    auto s1_zero_code{ vq::lbg<16, K>(mfcc)};

    std::string name("sx");
    name.at(1) = '0' + actual_speaker;
    speakers.push_back(speaker<16, K>{name});
    speakers.back().add_code(speaker<16, K>::Code{"zero", std::move(s1_zero_code)});
  }

  //load single test wav
  SF_INFO wavinfo;
  auto fi = argv[argc-1];
  auto wav_probes = load_wav(fi, wavinfo);
  std::cout << fi  << "...\t";
  sf_format_check(&wavinfo); 
  const auto samplerate = wavinfo.samplerate;

  auto speech_frames {dsp_utils::frame_signal<N>(wav_probes, M)};
  auto mfcc {mel_utils::mfcc_extraction<256, 30>(std::move(speech_frames), samplerate)};

  std::vector<double> euclidean_distances;
  for(const auto& speaker: speakers)
  {
    double distortion = vq::compute_distortion<16, 30>(speaker.codebook.at(0).centroids, mfcc);
    euclidean_distances.push_back(distortion); 
  }
  auto result = std::min_element(euclidean_distances.begin(), euclidean_distances.end());
  std::cout << "recognized speaker is: s" << std::distance(euclidean_distances.begin(), result) +1 << std::endl;

  dsp_utils::numeric_verifiaction();

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

