#include <mel_frame_generator.h>
#include <algorithm>

const constexpr int BUFFER_LEN = 1024;
const constexpr int FFT_SIZE = 256; 
const constexpr int N = FFT_SIZE;
const constexpr int FRAME_STEP = 100; 
const constexpr int M = FRAME_STEP;
const constexpr int MELL_FILTER_BANKS = 30;
const constexpr int K = MELL_FILTER_BANKS;
const constexpr int MFCC_NUM = 13;

std::array<double, K> mel_utils::mel_frame(const std::array<double, 256ul>& fr, int samplerate)
{
    std::array<double, K> mel_frame;
    std::generate(mel_frame.begin(), mel_frame.end(), mel_frame_generator<K, N>(fr, samplerate));
    return mel_frame;
}

inline double mel_utils::mel(double f)
{
  return 2595*std::log10(1+f/700);
}

inline double mel_utils::mel2freq(double mel)
{
  return 700*(std::pow(10, mel/2595)-1);
}
