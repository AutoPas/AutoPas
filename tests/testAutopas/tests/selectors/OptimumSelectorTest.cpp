/**
 * @file OptimumSelectorTest.cpp
 * @author F. Gratl
 * @date 3/15/19
 */

#include "OptimumSelectorTest.h"

#include "autopas/selectors/OptimumSelector.h"

TEST(OptimumSelectorTest, min) {
  std::vector<long> vals = {5, 6, 3, 1, 7};

  auto min = autopas::OptimumSelector::optimumValue(vals, autopas::SelectorStrategyOption::fastestAbs);

  EXPECT_EQ(1, min);
}

TEST(OptimumSelectorTest, mean) {
  std::vector<long> vals = {5, 6, 3, 1, 7};

  auto mean = autopas::OptimumSelector::optimumValue(vals, autopas::SelectorStrategyOption::fastestMean);

  EXPECT_EQ(4, mean);
}

TEST(OptimumSelectorTest, median) {
  std::vector<long> vals = {5, 6, 3, 1, 7};

  auto median = autopas::OptimumSelector::optimumValue(vals, autopas::SelectorStrategyOption::fastestMedian);

  EXPECT_EQ(5, median);
}

/**
 * Make sure no selector crashes on empty vectors
 */
TEST(OptimumSelectorTest, empty) {
  std::vector<long> vals{};

  for (const auto option : autopas::SelectorStrategyOption::getAllOptions()) {
    long value = 42;
    EXPECT_NO_THROW(value = autopas::OptimumSelector::optimumValue(vals, option)) << "Failed for option " << option;
    EXPECT_EQ(0, value) << "Failed for option " << option;
  }
}