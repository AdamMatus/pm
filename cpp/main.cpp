#include <iostream>
#include <iomanip>
#include <vector>

#include <gsl/gsl_sf_bessel.h>
#include <sndfile.h>
#include <cstring>

#include "matplotlibcpp.h"

namespace mpl = matplotlibcpp;

const constexpr auto BUFFER_LEN = 1024;

int main (void)
{
  double x = 5.0;
  double y = gsl_sf_bessel_J0(x);
  std::cout << "J0(" << static_cast<int>(x) << ")"  << " = " << std::setprecision(18) << y << std::endl;

  const char *file = "s1.wav";
  float data[BUFFER_LEN];

  SNDFILE *infile;
  SF_INFO  sfinfo;
  memset(&sfinfo, 0, sizeof(sfinfo));

  if(!(infile = sf_open(file, SFM_READ, &sfinfo)))
  {
    std::cerr << "Not able to load a file: " << file << std::endl;
    std::cerr << sf_strerror(NULL);
    return 1;
  }

  int readcount;
  std::vector<float> wav_probes;
  while((readcount = sf_read_float(infile, data, BUFFER_LEN)))
  {
    wav_probes.insert(wav_probes.end(), data, data + readcount);      
  }

  mpl::plot(wav_probes);
  mpl::show();

  return 0;
}
