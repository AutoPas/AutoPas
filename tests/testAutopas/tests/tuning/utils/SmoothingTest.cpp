/**
 * @file SmoothingTest.cpp
 * @author F. Gratl
 * @date 23/11/2020
 */

#include "SmoothingTest.h"

#include "autopas/tuning/searchSpace/Evidence.h"
#include "autopas/tuning/utils/Smoothing.h"

TEST(SmoothingTest, loessLastPoint) {
  std::vector<size_t> xvals{1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 8, 10, 12};
  std::vector<long> yvals{18000, 2000,  15000, 6000,  10000, 4000, 16000, 11000, 7000,
                          3000,  14000, 17000, 20000, 12000, 9000, 13000, 1000,  8000};

  std::vector<autopas::Evidence> obs;
  obs.reserve(xvals.size());
  for (size_t i = 0; i < xvals.size(); ++i) {
    obs.emplace_back(autopas::Evidence{xvals[i], 0, yvals[i]});
  }

  // YS values with F = .25, NSTEPS = 0, DELTA = 0.0
  {
    auto out = autopas::smoothing::smoothLastPoint(obs, (size_t)(.25 * (double)obs.size()));
    EXPECT_EQ(out, 5724);
  }
}