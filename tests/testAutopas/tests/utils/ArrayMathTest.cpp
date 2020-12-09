/**
 * @file ArrayMathTest.cpp
 * @author tchipev
 * @date 19.01.18
 */

#include <gtest/gtest.h>

#include "autopas/utils/ArrayMath.h"

namespace ArrayMathTest {

using namespace autopas;

TEST(ArrayMathTest, testAddInt) {
  std::array<int, 3> a({1, 2, 3}), b({4, 5, 6});
  std::array<int, 3> correctResult({5, 7, 9});
  std::array<int, 3> result = utils::ArrayMath::add(a, b);
  for (int d = 0; d < 3; ++d) {
    ASSERT_EQ(result[d], correctResult[d]);
  }
}

TEST(ArrayMathTest, testAddDouble) {
  std::array<double, 3> a({1.1, 2.2, 3.3}), b({4.4, 5.5, 6.6});
  std::array<double, 3> correctResult({5.5, 7.7, 9.9});
  std::array<double, 3> result = utils::ArrayMath::add(a, b);
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], correctResult[d]);
  }

  result = utils::ArrayMath::add(result, a);
  std::array<double, 3> correctResult2({6.6, 9.9, 13.2});
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], correctResult2[d]);
  }
}

TEST(ArrayMathTest, testMulDouble) {
  std::array<double, 3> a({1.1, 2.2, 3.3}), b({4.4, 5.5, 6.6});
  std::array<double, 3> correctResult({1.21, 4.84, 10.89});
  std::array<double, 3> result = utils::ArrayMath::mul(a, a);
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], correctResult[d]);
  }

  correctResult = std::array<double, 3>({4.84, 12.1, 21.78});
  result = utils::ArrayMath::mul(a, b);
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], correctResult[d]);
  }
}

TEST(ArrayMathTest, testDot) {
  std::array<double, 3> a({1.1, 2.2, 3.3}), b({4.4, 5.5, 6.6});
  double rd = utils::ArrayMath::dot(a, a);
  ASSERT_DOUBLE_EQ(rd, 16.94);
  ASSERT_DOUBLE_EQ(38.720, utils::ArrayMath::dot(a, b));
}

TEST(ArrayMathTest, testAddScalar) {
  std::array<double, 3> a({1.1, 2.2, 3.3});
  std::array<double, 3> correctResult({3.1, 4.2, 5.3});
  std::array<double, 3> result = utils::ArrayMath::addScalar(a, 2.0);
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], correctResult[d]);
  }
}

TEST(ArrayMathTest, testSubScalar) {
  std::array<double, 3> a({1.1, 2.2, 3.3});
  std::array<double, 3> correctResult({-0.9, 0.2, 1.3});
  std::array<double, 3> result = utils::ArrayMath::subScalar(a, 2.0);
  for (int d = 0; d < 3; ++d) {
    ASSERT_NEAR(result[d], correctResult[d], 1e-15);
  }
}

TEST(ArrayMathTest, testMulScalar) {
  std::array<double, 3> a({1.1, 2.2, 3.3});
  std::array<double, 3> correctResult({2.2, 4.4, 6.6});
  std::array<double, 3> result = utils::ArrayMath::mulScalar(a, 2.0);
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], correctResult[d]);
  }
}

TEST(ArrayMathTest, testL2Norm) {
  std::array<double, 3> a({1.1, -2.2, 3.3});
  const double correctResult = std::sqrt(16.94);
  const double result = utils::ArrayMath::L2Norm(a);
  ASSERT_DOUBLE_EQ(result, correctResult);
}

TEST(ArrayMathTest, testProd) {
  std::array<double, 3> a({1.1, -2.2, 3.3});
  const double correctResult = -7.986;
  const double result = utils::ArrayMath::prod(a);
  ASSERT_DOUBLE_EQ(result, correctResult);
}

TEST(ArrayMathTest, testNormalize) {
  std::array<double, 3> a({1.1, 2.2, 3.3});
  std::array<double, 3> correctResult({sqrt(14) / 14.0, sqrt(14) / 7.0, 3 / sqrt(14)});
  std::array<double, 3> result = utils::ArrayMath::normalize(a);
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], correctResult[d]);
  }
}

TEST(ArrayMathTest, testNormalizeNullVector) {
  std::array<double, 3> a({0.0, 0.0, 0.0});
  std::array<double, 3> result = utils::ArrayMath::normalize(a);
  for (int d = 0; d < 3; ++d) {
    ASSERT_TRUE(std::isnan(result[d]));
  }
}

}  // end namespace ArrayMathTest
