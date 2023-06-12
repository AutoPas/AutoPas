/**
 * @file ConstexprMathTest.cpp
 * @author D. Martin
 * @date 15.05.2023
 */

#include <gtest/gtest.h>

#include <cmath>

#include "autopas/utils/ConstexprMath.h"

using namespace autopas;

TEST(ConstexprMathTest, testFloatSqrtPosValue) {
  constexpr double x = 2.0;
  const double exp = std::sqrt(x);
  constexpr double res = utils::ConstexprMath::sqrt(x, 1e-16);
  ASSERT_DOUBLE_EQ(res, exp);
}

TEST(ConstexprMathTest, testFloatSqrtZero) {
  constexpr double x = 0;
  const double exp = std::sqrt(x);
  constexpr double res = utils::ConstexprMath::sqrt(x, 1e-16);
  ASSERT_DOUBLE_EQ(res, exp);
}

TEST(ConstexprMathTest, testFloatSqrtNegative) {
  constexpr double x = -1.;
  constexpr double res = utils::ConstexprMath::sqrt(x, 1e-16);
  ASSERT_TRUE(std::isnan(res));
}

TEST(ConstexprMathTest, testIntSqrtPos) {
  constexpr int x = 25;
  const int exp = std::sqrt(x);
  constexpr int res = utils::ConstexprMath::sqrt(x);
  ASSERT_EQ(res, exp);
}

TEST(ConstexprMathTest, testIntSqrtPos2) {
  constexpr int x = 17;
  const int exp = std::sqrt(x);
  constexpr int res = utils::ConstexprMath::sqrt(x);
  ASSERT_EQ(res, exp);
}

TEST(ConstexprMathTest, testIntSqrtZero) {
  constexpr int x = 0;
  const int exp = std::sqrt(x);
  constexpr int res = utils::ConstexprMath::sqrt(x);
  ASSERT_EQ(res, exp);
}

TEST(ConstexprMathTest, testIntSqrtNegative) {
  const int x = -1;
  EXPECT_ANY_THROW(utils::ConstexprMath::sqrt(x));
}