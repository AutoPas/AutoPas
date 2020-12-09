/**
 * @file OptimumSelectorTest.cpp
 * @author F. Gratl
 * @date 3/15/19
 */

#include "OptimumSelectorTest.h"

#include "autopas/selectors/OptimumSelector.h"

namespace OptimumSelectorTest {

TEST(OptimumSelectorTest, min) {
  std::vector<unsigned long> vals = {5, 6, 3, 1, 7};

  auto min = autopas::OptimumSelector::optimumValue(vals, autopas::SelectorStrategyOption::fastestAbs);

  EXPECT_EQ(1, min);
}

TEST(OptimumSelectorTest, mean) {
  std::vector<unsigned long> vals = {5, 6, 3, 1, 7};

  auto mean = autopas::OptimumSelector::optimumValue(vals, autopas::SelectorStrategyOption::fastestMean);

  EXPECT_EQ(4, mean);
}

TEST(OptimumSelectorTest, median) {
  std::vector<unsigned long> vals = {5, 6, 3, 1, 7};

  auto median = autopas::OptimumSelector::optimumValue(vals, autopas::SelectorStrategyOption::fastestMedian);

  EXPECT_EQ(5, median);
}

}  // end namespace OptimumSelectorTest
