/**
 * @file SmoothingTest.cpp
 * @author F. Gratl
 * @date 23/11/2020
 */

#include "SmoothingTest.h"

#include "autopas/selectors/Smoothing.h"

TEST(SmoothingTest, lowessLastPoint) {
  std::vector<size_t> v_xval{1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 8, 10, 12};
  std::vector<size_t> v_yval{18, 2, 15, 6, 10, 4, 16, 11, 7, 3, 14, 17, 20, 12, 9, 13, 1, 8};

  std::vector<std::pair<size_t, long>> obs;
  obs.reserve(v_xval.size());
  for (size_t i = 0; i < v_xval.size(); ++i) {
    obs.emplace_back(v_xval[i], v_yval[i]);
  }

  // YS values with F = .25, NSTEPS = 0, DELTA = 0.0
  std::cout << "THE TEST" << std::endl;
  {
    auto out = autopas::Smoothing::smoothLastPoint(obs, .25);
    // out should be an integer (or long)
    EXPECT_EQ(out, 6);
  }
}