/**
 * @file MathTest.cpp.cc
 * @author F. Gratl
 * @date 23.04.21
 */

#include "MathTest.h"

#include "autopas/utils/Math.h"

TYPED_TEST_SUITE_P(MathTest);

TYPED_TEST_P(MathTest, safeAddTest) {
  TypeParam a = 200;
  TypeParam b = 4000;

  EXPECT_EQ(a + b, autopas::utils::Math::safeAdd(a, b)) << "Simple addition (+ +) failed.";
  EXPECT_EQ(std::numeric_limits<TypeParam>::max(),
            autopas::utils::Math::safeAdd(std::numeric_limits<TypeParam>::max(), std::numeric_limits<TypeParam>::max()))
      << "Overflowing addition produced unexpected result.";

  // tests with signed numbers
  // underflow can only happen with signed types
  if constexpr (std::is_signed_v<TypeParam>) {
    EXPECT_EQ(-a + b, autopas::utils::Math::safeAdd(-a, b)) << "Simple addition (- +) failed.";
    EXPECT_EQ(a + -b, autopas::utils::Math::safeAdd(a, -b)) << "Simple addition (+ -) failed.";
    if constexpr (std::is_integral_v<TypeParam>) {
      EXPECT_EQ(
          std::numeric_limits<TypeParam>::min(),
          autopas::utils::Math::safeAdd(std::numeric_limits<TypeParam>::min(), std::numeric_limits<TypeParam>::min()))
          << "Underflowing addition produced unexpected result.";
    } else if constexpr (std::is_floating_point_v<TypeParam>) {
      // for floating point numbers, numeric_limits::min() returns the positive value closest to zero
      // However, since float range is guaranteed to be symmetrical we can use -numeric_limits::max().
      EXPECT_EQ(
          -std::numeric_limits<TypeParam>::max(),
          autopas::utils::Math::safeAdd(-std::numeric_limits<TypeParam>::max(), -std::numeric_limits<TypeParam>::max()))
          << "Underflowing addition produced unexpected result.";
    }
  }
}

TYPED_TEST_P(MathTest, safeSubTest) {
  TypeParam a = 5000;
  TypeParam b = 200;

  EXPECT_EQ(a - b, autopas::utils::Math::safeSub(a, b)) << "Simple subtraction (+ +) failed.";

  // tests with signed numbers
  // underflow can only happen with signed types
  if constexpr (std::is_signed_v<TypeParam>) {
    EXPECT_EQ(-a - b, autopas::utils::Math::safeSub(-a, b)) << "Simple subtraction (- +) failed.";
    EXPECT_EQ(a - -b, autopas::utils::Math::safeSub(a, -b)) << "Simple subtraction (+ -) failed.";
    EXPECT_EQ(
        std::numeric_limits<TypeParam>::max(),
        autopas::utils::Math::safeSub(std::numeric_limits<TypeParam>::max(), -std::numeric_limits<TypeParam>::max()))
        << "Overflowing subtraction produced unexpected result.";
    if constexpr (std::is_integral_v<TypeParam>) {
      EXPECT_EQ(
          std::numeric_limits<TypeParam>::min(),
          autopas::utils::Math::safeSub(std::numeric_limits<TypeParam>::min(), std::numeric_limits<TypeParam>::max()))
          << "Underflowing subtraction produced unexpected result.";
    } else if constexpr (std::is_floating_point_v<TypeParam>) {
      // for floating point numbers, numeric_limits::min() returns the positive value closest to zero
      // However, since float range is guaranteed to be symmetrical we can use -numeric_limits::max().
      EXPECT_EQ(
          -std::numeric_limits<TypeParam>::max(),
          autopas::utils::Math::safeSub(-std::numeric_limits<TypeParam>::max(), std::numeric_limits<TypeParam>::max()))
          << "Underflowing subtraction produced unexpected result.";
    }
  } else {  // tests for unsigned numbers
    EXPECT_EQ(0, autopas::utils::Math::safeSub(b, a)) << "Unsigned Subtraction underflowing zero failed.";
  }
}

TYPED_TEST_P(MathTest, safeMulTest) {
  TypeParam a = 200;
  TypeParam b = 4000;

  EXPECT_EQ(a * b, autopas::utils::Math::safeMul(a, b)) << "Simple multiplication (+ +) failed.";
  EXPECT_EQ(std::numeric_limits<TypeParam>::max(),
            autopas::utils::Math::safeMul(std::numeric_limits<TypeParam>::max(), std::numeric_limits<TypeParam>::max()))
      << "Overflowing multiplication (+ +) produced unexpected result.";

  // tests with signed numbers
  // underflow can only happen with signed types
  if constexpr (std::is_signed_v<TypeParam>) {
    EXPECT_EQ(-a * b, autopas::utils::Math::safeMul(-a, b)) << "Simple multiplication (- +) failed.";
    EXPECT_EQ(a * -b, autopas::utils::Math::safeMul(a, -b)) << "Simple multiplication (+ -) failed.";
    EXPECT_EQ(
        std::numeric_limits<TypeParam>::max(),
        autopas::utils::Math::safeMul(-std::numeric_limits<TypeParam>::max(), -std::numeric_limits<TypeParam>::max()))
        << "Overflowing multiplication (- -) produced unexpected result.";
    if constexpr (std::is_integral_v<TypeParam>) {
      EXPECT_EQ(
          std::numeric_limits<TypeParam>::min(),
          autopas::utils::Math::safeMul(std::numeric_limits<TypeParam>::max(), std::numeric_limits<TypeParam>::min()))
          << "Underflowing multiplication (+ -) produced unexpected result.";
      EXPECT_EQ(
          std::numeric_limits<TypeParam>::min(),
          autopas::utils::Math::safeMul(std::numeric_limits<TypeParam>::min(), std::numeric_limits<TypeParam>::max()))
          << "Underflowing multiplication (- +) produced unexpected result.";
    } else if constexpr (std::is_floating_point_v<TypeParam>) {
      // for floating point numbers, numeric_limits::min() returns the positive value closest to zero
      // However, since float range is guaranteed to be symmetrical we can use -numeric_limits::max().
      EXPECT_EQ(
          -std::numeric_limits<TypeParam>::max(),
          autopas::utils::Math::safeMul(std::numeric_limits<TypeParam>::max(), -std::numeric_limits<TypeParam>::max()))
          << "Underflowing multiplication (+ -) produced unexpected result.";
      EXPECT_EQ(
          -std::numeric_limits<TypeParam>::max(),
          autopas::utils::Math::safeMul(-std::numeric_limits<TypeParam>::max(), std::numeric_limits<TypeParam>::max()))
          << "Underflowing multiplication (- +) produced unexpected result.";
    }
  }
}

REGISTER_TYPED_TEST_SUITE_P(MathTest, safeAddTest, safeSubTest, safeMulTest);

// All types that are at least of size int:
using TestedTypes = ::testing::Types<unsigned int, int, long, unsigned long, float, double>;
INSTANTIATE_TYPED_TEST_SUITE_P(GeneratedTyped, MathTest, TestedTypes);