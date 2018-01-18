/*
 * arrayMath.h
 *
 *  Created on: 18 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_UTILS_ARRAYMATH_H_
#define SRC_UTILS_ARRAYMATH_H_

#include <array>

namespace autopas {

template<class T, std::size_t SIZE>
std::array<T, SIZE> arrayAdd(const std::array<T, SIZE>& a, const std::array<T, SIZE>& b) {
	std::array<T, SIZE> result;
	for (std::size_t d = 0; d < SIZE; ++d) {
		result[d] = a[d] + b[d];
	}
	return result;
}

template<class T, std::size_t SIZE>
std::array<T, SIZE> arraySub(const std::array<T, SIZE>& a, const std::array<T, SIZE>& b) {
	std::array<T, SIZE> result;
	for (std::size_t d = 0; d < SIZE; ++d) {
		result[d] = a[d] - b[d];
	}
	return result;
}

template<class T, std::size_t SIZE>
std::array<T, SIZE> arrayMul(const std::array<T, SIZE>& a, const std::array<T, SIZE>& b) {
	std::array<T, SIZE> result;
	for (std::size_t d = 0; d < SIZE; ++d) {
		result[d] = a[d] * b[d];
	}
	return result;
}

template<class T, std::size_t SIZE>
T arrayDot(const std::array<T, SIZE>& a, const std::array<T, SIZE>& b) {
	T result = static_cast<T>(0.0);
	for (std::size_t d = 0; d < SIZE; ++d) {
		result += a[d] * b[d];
	}
	return result;
}

} /* namespace autopas */




#endif /* SRC_UTILS_ARRAYMATH_H_ */
