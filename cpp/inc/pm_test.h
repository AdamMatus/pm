#ifndef _PM_TEST_
#define _PM_TEST_

#include "matplotlibcpp.h"
#include <vector>

namespace mpl = matplotlibcpp;

template <typename It>
void test_mpl(It start, It end, int start_index =0)
{
  typedef typename std::iterator_traits<It>::value_type T;

  auto t = std::vector<T>(end-start);
  std::iota(t.begin(), t.end(), start_index);

  auto y = std::vector<T>(start, end);

  mpl::plot(t, y);
  mpl::show();
}

template <typename It>
void two_dim(It tstart, It tend, It start, It end)
{
  typedef typename std::iterator_traits<It>::value_type T;

  auto y = std::vector<T>(start, end);
  auto t = std::vector<T>(tstart, tend);

  mpl::plot(t, y, ".");
  mpl::show();
}

#endif

