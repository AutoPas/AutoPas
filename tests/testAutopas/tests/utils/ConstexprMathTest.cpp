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
  double x = 2.0;
  double exp = std::sqrt(x);
  ASSERT_DOUBLE_EQ(utils::ConstexprMath::sqrt(x), exp);
}

TEST(ConstexprMathTest, testFloatSqrtZero) {
  double x = 0;
  double exp = std::sqrt(x);
  ASSERT_DOUBLE_EQ(utils::ConstexprMath::sqrt(x), exp);
}

TEST(ConstexprMathTest, testFloatSqrtNegative) {
  double x = -1.;
  ASSERT_TRUE(std::isnan(utils::ConstexprMath::sqrt(x)));
}

TEST(ConstexprMathTest, testIntSqrtPos) {
  int x = 25;
  int exp = std::sqrt(x);
  ASSERT_EQ(utils::ConstexprMath::sqrt(x), exp);
}

TEST(ConstexprMathTest, testIntSqrtPos2) {
  int x = 17;
  int exp = std::sqrt(x);
  ASSERT_EQ(utils::ConstexprMath::sqrt(x), exp);
}

TEST(ConstexprMathTest, testIntSqrtZero) {
  int x = 0;
  int exp = std::sqrt(x);
  ASSERT_EQ(utils::ConstexprMath::sqrt(x), exp);
}

TEST(ConstexprMathTest, testIntSqrtNegative) {
  int x = -1;
  EXPECT_ANY_THROW(utils::ConstexprMath::sqrt(x));
}