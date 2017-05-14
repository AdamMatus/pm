#include "vq.h"

#include <algorithm>
#include <cmath>
#include <tuple>

#include "test_inc/mfcc_test.h"

void vq::dummy_test()
{
  codebook_array<16,29> pies { vq::lbg<16,29>(mfcc_test_frames)};
}

