/**
 * @file ArrayMathTest.cpp
 * @author tchipev
 * @date 19.01.18
 */

#include <gtest/gtest.h>
#include "autopas/utils/ArrayMath.h"

using namespace autopas;

TEST(ArrayMathTest, testAddInt) {
  std::array<int, 3> a({1, 2, 3}), b({4, 5, 6});
  std::array<int, 3> correctResult({5, 7, 9});
  std::array<int, 3> result = ArrayMath::add(a, b);
  for (int d = 0; d < 3; ++d) {
    ASSERT_EQ(result[d], correctResult[d]);
  }
}

TEST(ArrayMathTest, testAddDouble) {
  std::array<double, 3> a({1.1, 2.2, 3.3}), b({4.4, 5.5, 6.6});
  std::array<double, 3> correctResult({5.5, 7.7, 9.9});
  std::array<double, 3> result = ArrayMath::add(a, b);
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], correctResult[d]);
  }

  result = ArrayMath::add(result, a);
  std::array<double, 3> correctResult2({6.6, 9.9, 13.2});
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], correctResult2[d]);
  }
}

TEST(ArrayMathTest, testMulDouble) {
  std::array<double, 3> a({1.1, 2.2, 3.3}), b({4.4, 5.5, 6.6});
  std::array<double, 3> correctResult({1.21, 4.84, 10.89});
  std::array<double, 3> result = ArrayMath::mul(a, a);
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], correctResult[d]);
  }

  correctResult = std::array<double, 3>({4.84, 12.1, 21.78});
  result = ArrayMath::mul(a, b);
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], correctResult[d]);
  }
}

TEST(ArrayMathTest, testDot) {
  std::array<double, 3> a({1.1, 2.2, 3.3}), b({4.4, 5.5, 6.6});
  double rd = ArrayMath::dot(a, a);
  ASSERT_DOUBLE_EQ(rd, 16.94);
  ASSERT_DOUBLE_EQ(38.720, ArrayMath::dot(a, b));
}

TEST(ArrayMathTest, testAddScalar) {
  std::array<double, 3> a({1.1, 2.2, 3.3});
  std::array<double, 3> correctResult({3.1, 4.2, 5.3});
  std::array<double, 3> result = ArrayMath::addScalar(a, 2.0);
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], correctResult[d]);
  }
}

TEST(ArrayMathTest, testSubScalar) {
  std::array<double, 3> a({1.1, 2.2, 3.3});
  std::array<double, 3> correctResult({-0.9, 0.2, 1.3});
  std::array<double, 3> result = ArrayMath::subScalar(a, 2.0);
  for (int d = 0; d < 3; ++d) {
    ASSERT_NEAR(result[d], correctResult[d], 1e-15);
  }
}

TEST(ArrayMathTest, testMulScalar) {
  std::array<double, 3> a({1.1, 2.2, 3.3});
  std::array<double, 3> correctResult({2.2, 4.4, 6.6});
  std::array<double, 3> result = ArrayMath::mulScalar(a, 2.0);
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], correctResult[d]);
  }
}

TEST(ArrayMathTest, teststatic_cast_array) {
  {
    std::array<long, 3> in({1l, 2l, 3l});
    auto out = ArrayMath::static_cast_array<unsigned long>(in);
    static_assert(std::is_same<decltype(out), std::array<unsigned long, 3>>::value, "Type mismatch");
  }

  {
    std::array<long, 3> in({1l, 2l, 3l});
    auto out = ArrayMath::static_cast_array<double>(in);
    static_assert(std::is_same<decltype(out), std::array<double, 3>>::value, "Type mismatch");
  }
}
