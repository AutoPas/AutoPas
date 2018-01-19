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

namespace arrayMath {

template<class T, std::size_t SIZE>
std::array<T, SIZE> add(const std::array<T, SIZE>& a, const std::array<T, SIZE>& b) {
	std::array<T, SIZE> result;
	for (std::size_t d = 0; d < SIZE; ++d) {
		result[d] = a[d] + b[d];
	}
	return result;
}

template<class T, std::size_t SIZE>
std::array<T, SIZE> sub(const std::array<T, SIZE>& a, const std::array<T, SIZE>& b) {
	std::array<T, SIZE> result;
	for (std::size_t d = 0; d < SIZE; ++d) {
		result[d] = a[d] - b[d];
	}
	return result;
}

template<class T, std::size_t SIZE>
std::array<T, SIZE> mul(const std::array<T, SIZE>& a, const std::array<T, SIZE>& b) {
	std::array<T, SIZE> result;
	for (std::size_t d = 0; d < SIZE; ++d) {
		result[d] = a[d] * b[d];
	}
	return result;
}

template<class T, std::size_t SIZE>
std::array<T, SIZE> addScalar(const std::array<T, SIZE>& a, T s) {
	std::array<T, SIZE> result;
	for (std::size_t d = 0; d < SIZE; ++d) {
		result[d] = a[d] + s;
	}
	return result;
}

template<class T, std::size_t SIZE>
std::array<T, SIZE> mulScalar(const std::array<T, SIZE>& a, T s) {
	std::array<T, SIZE> result;
	for (std::size_t d = 0; d < SIZE; ++d) {
		result[d] = a[d] * s;
	}
	return result;
}

template<class T, std::size_t SIZE>
T dot(const std::array<T, SIZE>& a, const std::array<T, SIZE>& b) {
	T result = static_cast<T>(0.0);
	for (std::size_t d = 0; d < SIZE; ++d) {
		result += a[d] * b[d];
	}
	return result;
}

} /* namespace arrayMath */

} /* namespace autopas */




#endif /* SRC_UTILS_ARRAYMATH_H_ */
