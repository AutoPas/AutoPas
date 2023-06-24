/**
 * @file SmoothingTest.cpp
 * @author F. Gratl
 * @date 23/11/2020
 */

#include "SmoothingTest.h"

#include "autopas/tuning/utils/Evidence.h"
#include "autopas/tuning/utils/Smoothing.h"

TEST(SmoothingTest, lowessLastPoint) {
  std::vector<size_t> xvals{1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 8, 10, 12};
  std::vector<long> yvals{18, 2, 15, 6, 10, 4, 16, 11, 7, 3, 14, 17, 20, 12, 9, 13, 1, 8};

  std::vector<autopas::Evidence> obs;
  obs.reserve(xvals.size());
  for (size_t i = 0; i < xvals.size(); ++i) {
    obs.emplace_back(xvals[i], yvals[i]);
  }

  // YS values with F = .25, NSTEPS = 0, DELTA = 0.0
  {
    auto out = autopas::smoothing::smoothLastPoint(obs, (.25 * obs.size()));
    EXPECT_EQ(out, 6);
  }
}