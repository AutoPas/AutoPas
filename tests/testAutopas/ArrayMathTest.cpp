/*
 * ArrayMathTest.cpp
 *
 *  Created on: 19 Jan 2018
 *      Author: tchipevn
 */

#include "autopas.h"
#include "gtest/gtest.h"

using namespace autopas;

TEST(ArrayMathTest, testAddInt) {
	std::array<int, 3> a({1, 2, 3}), b({4, 5, 6});
	std::array<int, 3> correctResult({5, 7, 9});
	std::array<int, 3> result = arrayMath::add(a,b);
	for (int d = 0; d < 3; ++d) {
		ASSERT_EQ(result[d], correctResult[d]);
	}
}

TEST(ArrayMathTest, testAddDouble) {
	std::array<double, 3> a({1.1, 2.2, 3.3}), b({4.4, 5.5, 6.6});
	std::array<double, 3> correctResult({5.5, 7.7, 9.9});
	std::array<double, 3> result = arrayMath::add(a,b);
	for (int d = 0; d < 3; ++d) {
		ASSERT_DOUBLE_EQ(result[d], correctResult[d]);
	}

	result = arrayMath::add(result, a);
	std::array<double, 3> correctResult2({6.6, 9.9, 13.2});
	for (int d = 0; d < 3; ++d) {
		ASSERT_DOUBLE_EQ(result[d], correctResult2[d]);
	}
}

TEST(ArrayMathTest, testMulDouble) {
	std::array<double, 3> a({1.1, 2.2, 3.3}), b({4.4, 5.5, 6.6});
	std::array<double, 3> correctResult({1.21, 4.84, 10.89});
	std::array<double, 3> result = arrayMath::mul(a,a);
	for (int d = 0; d < 3; ++d) {
		ASSERT_DOUBLE_EQ(result[d], correctResult[d]);
	}

	correctResult = std::array<double, 3>({4.84, 12.1, 21.78});
	result = arrayMath::mul(a,b);
	for (int d = 0; d < 3; ++d) {
		ASSERT_DOUBLE_EQ(result[d], correctResult[d]);
	}
}

TEST(ArrayMathTest, testDot) {
	std::array<double, 3> a({1.1, 2.2, 3.3}), b({4.4, 5.5, 6.6});
	double rd = arrayMath::dot(a, a);
	ASSERT_DOUBLE_EQ(rd, 16.94);
	ASSERT_DOUBLE_EQ(38.720, arrayMath::dot(a,b));
}

TEST(ArrayMathTest, testMulScalar) {
	std::array<double, 3> a({1.1, 2.2, 3.3});
	std::array<double, 3> correctResult({2.2, 4.4, 6.6});
	std::array<double, 3> result = arrayMath::mulScalar(a, 2.0);
	for (int d = 0; d < 3; ++d) {
		ASSERT_DOUBLE_EQ(result[d], correctResult[d]);
	}
}
