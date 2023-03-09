/**
 * @file ArrayUtilsTest.cpp
 * @author tchipev
 * @date 19.01.18
 */

#include <gtest/gtest.h>

#include <random>

#include "autopas/utils/ArrayUtils.h"

using namespace autopas;

TEST(ArrayUtilsTest, teststatic_cast_array) {
  {
    std::array<long, 3> in({1l, 2l, 3l});
    auto out = utils::ArrayUtils::static_cast_array<unsigned long>(in);
    static_assert(std::is_same<decltype(out), std::array<unsigned long, 3>>::value, "Type mismatch");
  }

  {
    std::array<long, 3> in({1l, 2l, 3l});
    auto out = utils::ArrayUtils::static_cast_array<double>(in);
    static_assert(std::is_same<decltype(out), std::array<double, 3>>::value, "Type mismatch");
  }
}

TEST(ArrayUtilsTest, testto_string) {
  std::array<size_t, 0> emptyContainer{};
  EXPECT_EQ("[]", utils::ArrayUtils::to_string(emptyContainer));

  std::array<size_t, 3> arrayContainer{1, 2, 3};
  EXPECT_EQ("[1, 2, 3]", utils::ArrayUtils::to_string(arrayContainer));
  EXPECT_EQ("1x2x3", utils::ArrayUtils::to_string(arrayContainer, "x", {"", ""}));

  std::vector<size_t> vectorContainer{1, 2, 3};
  EXPECT_EQ("[1, 2, 3]", utils::ArrayUtils::to_string(vectorContainer));
  EXPECT_EQ("1x2x3", utils::ArrayUtils::to_string(vectorContainer, "x", {"", ""}));

  std::basic_string bStringContainer("123");
  EXPECT_EQ("[1, 2, 3]", utils::ArrayUtils::to_string(bStringContainer));
  EXPECT_EQ("1x2x3", utils::ArrayUtils::to_string(bStringContainer, "x", {"", ""}));
}

/**
 * Test utils::ArrayUtils::balanceVectors.
 * Creates a 2D Vector, rebalances it and checks if the lengths are what we expect.
 */
TEST(ArrayUtilsTest, testBalanceVectors) {
  constexpr size_t numVectors = 8;
  // Setup some vectors of random length between 0 and 10
  std::mt19937 g(42);
  std::uniform_int_distribution dist(0, 10);
  std::array<std::vector<size_t>, numVectors> vecvec{};
  for (auto &v : vecvec) {
    v.resize(dist(g));
  }

  // count the number of "inserted" elements
  const auto numElements = [&]() {
    size_t sum = 0;
    for (const auto &vec : vecvec) {
      sum += vec.size();
    }
    return sum;
  }();
  // calculate expectations
  const auto expectedPerVec = numElements / vecvec.size();
  const auto rest = numElements % vecvec.size();

  // actual test starts here
  utils::ArrayUtils::balanceVectors(vecvec);

  for (size_t i = 0; i < vecvec.size(); ++i) {
    // the first vectors might have one extra element if the total number was not divisible by the number of vectors
    EXPECT_EQ(vecvec[i].size(), i < rest ? expectedPerVec + 1 : expectedPerVec);
  }
}