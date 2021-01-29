/**
 * @file ArrayUtilsTest.cpp
 * @author tchipev
 * @date 19.01.18
 */

#include <gtest/gtest.h>

#include "autopas/utils/ArrayUtils.h"

using namespace autopas;

TEST(ArrayUtilsTest, teststatic_cast_array) {
  {
    std::array<long, 3> in({1l, 2l, 3l});
    auto out = utils::ArrayUtils::static_cast_array<uint64_t>(in);
    static_assert(std::is_same<decltype(out), std::array<uint64_t, 3>>::value, "Type mismatch");
  }

  {
    std::array<long, 3> in({1l, 2l, 3l});
    auto out = utils::ArrayUtils::static_cast_array<double>(in);
    static_assert(std::is_same<decltype(out), std::array<double, 3>>::value, "Type mismatch");
  }
}

TEST(ArrayUtilsTest, testto_string) {
  std::array<size_t, 0> emptyContainer;
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
