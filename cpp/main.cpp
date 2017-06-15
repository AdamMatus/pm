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

  SNDFILE *infile;

const constexpr int BUFFER_LEN = 1024;
const constexpr int FFT_SIZE = 256; 
const constexpr int N = FFT_SIZE;
const constexpr int FRAME_STEP = 100; 
const constexpr int M = FRAME_STEP;
const constexpr int MELL_FILTER_BANKS = 30;
const constexpr int K = MELL_FILTER_BANKS;
const constexpr int MFCC_NUM = 13;

void static_test_pm(const std::vector<speaker<16, K>>& speakers);

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

    sf_close(infile);
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
    double distortion = vq::compute_distortion<16, 30>(speaker.codebook.at(0).centroids, mfcc)/mfcc.size();
    euclidean_distances.push_back(distortion); 
  }
  auto result = std::min_element(euclidean_distances.begin(), euclidean_distances.end());
  std::cout << "recognized speaker is: s" << std::distance(euclidean_distances.begin(), result) +1 << std::endl;
  std::cout << *result << std::endl;

  dsp_utils::numeric_verifiaction();

  static_test_pm(speakers);

  return 0;
}

void static_test_pm(const std::vector<speaker<16, 30>>& speakers)
{
  const char* files[] = { "test/s1.wav",
                          "test/s2.wav",
                          "test/s3.wav",
                          "test/s4.wav",
                          "test/s5.wav",
                          "test/s6.wav",
                          "test/s7.wav",
                          "test/s8.wav",
                          "test/s9.wav",
                          "test/s10.wav",
                          "test/s11.wav"
                       //   "test/s12.wav",
                       //   "test/s13.wav",
                       //   "test/s14.wav",
                       //   "test/s15.wav",
                       //   "test/s16.wav"
  };
  const char* adamfiles[] = { "adam/s1.wav",
                          "adam/s2.wav",
                          "adam/s3.wav",
                          "adam/s4.wav",
                          "adam/s5.wav",
                          "adam/s6.wav",
                          "adam/s7.wav",
                          "adam/s8.wav",
                          "adam/s9.wav",
                          "adam/s10.wav",
                          "adam/s11.wav",
                          "adam/s12.wav",
                          "adam/s13.wav",
                          "adam/s14.wav",
                          "adam/s15.wav",
                          "adam/s16.wav"
  };

std::cout << "p=[";
for(double d = 2; d < 5; d+=0.05)
{
  int na=0;
  int nr=0;
  for(int i=0; i<11; ++i)
  {
    SF_INFO wavinfo;
    auto fi = files[i];
    auto wav_probes = load_wav(fi, wavinfo);
    sf_format_check(&wavinfo); 
    const auto samplerate = wavinfo.samplerate;

    auto speech_frames {dsp_utils::frame_signal<N>(wav_probes, M)};
    auto mfcc {mel_utils::mfcc_extraction<256, 30>(std::move(speech_frames), samplerate)};

    int j = 0;
    for(const auto& speaker: speakers)
    {
      double distortion = vq::compute_distortion<16, 30>(speaker.codebook.at(0).centroids, mfcc)/mfcc.size();
      auto decision = distortion < d;
      if(decision)
      {
        if(i!=j) ++na;
      }
      j++;
    }
    sf_close(infile);
  }

  for(int i=0; i<16; ++i)
  {
    SF_INFO wavinfo;
    auto fi = adamfiles[i];
    auto wav_probes = load_wav(fi, wavinfo);
    sf_format_check(&wavinfo); 
    const auto samplerate = wavinfo.samplerate;

    auto speech_frames {dsp_utils::frame_signal<N>(wav_probes, M)};
    auto mfcc {mel_utils::mfcc_extraction<256, 30>(std::move(speech_frames), samplerate)};

      double distortion = vq::compute_distortion<16, 30>(speakers.at(11).codebook.at(0).centroids, mfcc)/mfcc.size();
      auto decision = distortion < d;
      if(!decision)
      {
        ++nr;
      }

    sf_close(infile);
  }
  std::cout << 100.0*na/(16.0*15.0) <<','<<nr/0.16 << ';';
}
std::cout << ']';
} 

std::vector<double> load_wav(const char* file, SF_INFO&  sfinfo)
{
  double data[BUFFER_LEN];

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

