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
  std::array<int, 3> result = utils::ArrayMath::add(a, b);
  for (int d = 0; d < 3; ++d) {
    ASSERT_EQ(result[d], a[d] + b[d]);
  }
}

TEST(ArrayMathTest, testAddDouble) {
  std::array<double, 3> a({1.1, 2.2, 3.3}), b({4.4, 5.5, 6.6});
  std::array<double, 3> result = utils::ArrayMath::add(a, b);
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], a[d] + b[d]);
  }

  result = utils::ArrayMath::add(result, a);
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], a[d] + b[d] + a[d]);
  }
}

TEST(ArrayMathTest, testAddOpInt) {
  using namespace autopas::utils::ArrayMath::literals;
  std::array<int, 3> a({1, 2, 3}), b({4, 5, 6});
  std::array<int, 3> result = a + b;
  for (int d = 0; d < 3; ++d) {
    ASSERT_EQ(result[d], a[d] + b[d]);
  }

  result += a;
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], a[d] + b[d] + a[d]);
  }
}

TEST(ArrayMathTest, testAddOpDouble) {
  using namespace autopas::utils::ArrayMath::literals;
  std::array<double, 3> a({1.1, 2.2, 3.3}), b({4.4, 5.5, 6.6});
  std::array<double, 3> result = a + b;
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], a[d] + b[d]);
  }

  result += a;
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], a[d] + b[d] + a[d]);
  }
}

TEST(ArrayMathTest, testMulDouble) {
  std::array<double, 3> a({1.1, 2.2, 3.3}), b({4.4, 5.5, 6.6});
  std::array<double, 3> result = utils::ArrayMath::mul(a, a);
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], a[d] * a[d]);
  }

  result = utils::ArrayMath::mul(a, b);
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], a[d] * b[d]);
  }
}

TEST(ArrayMathTest, testMulOpDouble) {
  using namespace autopas::utils::ArrayMath::literals;

  std::array<double, 3> a({1.1, 2.2, 3.3}), b({4.4, 5.5, 6.6});
  std::array<double, 3> result = a * a;
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], a[d] * a[d]);
  }

  result *= b;
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], a[d] * a[d] * b[d]);
  }
}

TEST(ArrayMathTest, testDot) {
  std::array<double, 3> a({1.1, 2.2, 3.3}), b({4.4, 5.5, 6.6});
  double rd = utils::ArrayMath::dot(a, a);
  ASSERT_DOUBLE_EQ(rd, a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
  ASSERT_DOUBLE_EQ(a[0] * b[0] + a[1] * b[1] + a[2] * b[2], utils::ArrayMath::dot(a, b));
}

TEST(ArrayMathTest, testAddScalar) {
  std::array<double, 3> a({1.1, 2.2, 3.3});
  double scalar = 2.0;
  std::array<double, 3> result = utils::ArrayMath::addScalar(a, scalar);
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], a[d] + scalar);
  }
}

TEST(ArrayMathTest, testAddScalarOp) {
  using namespace autopas::utils::ArrayMath::literals;

  std::array<double, 3> a({1.1, 2.2, 3.3});
  double scalar = 2.0;
  std::array<double, 3> result = a + scalar;
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], a[d] + scalar);
  }

  result += 2.0;
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], a[d] + scalar + scalar);
  }
}

TEST(ArrayMathTest, testSubScalar) {
  std::array<double, 3> a({1.1, 2.2, 3.3});
  double scalar = 2.0;
  std::array<double, 3> result = utils::ArrayMath::subScalar(a, scalar);
  for (int d = 0; d < 3; ++d) {
    ASSERT_NEAR(result[d], a[d] - scalar, 1e-15);
  }
}

TEST(ArrayMathTest, testSubScalarOp) {
  using namespace autopas::utils::ArrayMath::literals;

  std::array<double, 3> a({1.1, 2.2, 3.3});
  double scalar = 2.0;
  std::array<double, 3> result = a - scalar;
  for (int d = 0; d < 3; ++d) {
    ASSERT_NEAR(result[d], a[d] - scalar, 1e-15);
  }

  result -= scalar;
  for (int d = 0; d < 3; ++d) {
    ASSERT_NEAR(result[d], a[d] - scalar - scalar, 1e-15);
  }
}

TEST(ArrayMathTest, testMulScalar) {
  std::array<double, 3> a({1.1, 2.2, 3.3});
  double scalar = 2.0;
  std::array<double, 3> result = utils::ArrayMath::mulScalar(a, scalar);
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], a[d] * scalar);
  }
}

TEST(ArrayMathTest, testMulScalarOp) {
  using namespace autopas::utils::ArrayMath::literals;

  std::array<double, 3> a({1.1, 2.2, 3.3});
  double scalar = 2.0;
  std::array<double, 3> result = a * scalar;
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], a[d] * scalar);
  }

  result *= scalar;
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], a[d] * scalar * scalar);
  }
}

TEST(ArrayMathTest, testL2Norm) {
  std::array<double, 3> a({1.1, -2.2, 3.3});
  const double result = utils::ArrayMath::L2Norm(a);
  ASSERT_DOUBLE_EQ(result, std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]));
}

TEST(ArrayMathTest, testProd) {
  std::array<double, 3> a({1.1, -2.2, 3.3});
  const double result = utils::ArrayMath::prod(a);
  ASSERT_DOUBLE_EQ(result, a[0] * a[1] * a[2]);
}

TEST(ArrayMathTest, testNormalize) {
  std::array<double, 3> a({1.1, 2.2, 3.3});
  std::array<double, 3> result = utils::ArrayMath::normalize(a);
  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(result[d], a[d] * (1.0 / std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])));
  }
}

TEST(ArrayMathTest, testNormalizeNullVector) {
  std::array<double, 3> a({0.0, 0.0, 0.0});
  std::array<double, 3> result = utils::ArrayMath::normalize(a);
  for (int d = 0; d < 3; ++d) {
    ASSERT_TRUE(std::isnan(result[d]));
  }
}
