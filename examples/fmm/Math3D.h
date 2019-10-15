/**
 * @file Math3D.h
 * @date 23.09.19
 * @author Joachim Marin
 */

#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <vector>

namespace Math3D {

constexpr auto maxFactorialParameter = 128;

/**
 * Must be called before using functions from Math3D.h.
 */
void initMath();

class MathOpt {
 public:
  static std::vector<double> factorialValue;
  static std::vector<double> doubleFactorialValue;
  static std::vector<double> legendreCache;
  static std::vector<std::vector<double>> getAValue;
  static int legendreLastM;
  static std::complex<double> sphericalCache;
  static std::vector<std::complex<double>> imaginaryPower;
};

/**
 * Returns a copy of a 3d vector in cartesian coordinates using spherical coordinates.
 * @param 3d vector in cartesian coordinates
 * @return copied 3d vector in spherical coordinates
 */
std::array<double, 3> toSpherical(const std::array<double, 3> &cartesian);

inline double doubleFactorial(int x) { return MathOpt::doubleFactorialValue.at(x); }

inline double factorial(int x) { return MathOpt::factorialValue.at(x); }

/**
 * Calculates the associated Legendre polynomial of order m and degree n at x.
 * @param m
 * @param n
 * @param x
 * @return
 */
double associatedLegendrePolynomial(int m, int n, double x);

/**
 * After associatedLegendrePolynomial has been called for a certain m and x. This function can be used to get the value
 * of the associated Legendre polynomial of same order m and same value x, but smaller n.
 * @param n Degree of the associated legendre polynomial.
 * @return
 */
inline double getLegendreCache(int n) { return MathOpt::legendreCache[n - MathOpt::legendreLastM]; }

/**
 * Calculates the spherical harmonics (Y) using associated legendre polynomials.
 * @param m (superscript)
 * @param n (subscript)
 * @param theta (first parameter)
 * @param phi (second parameter)
 * @return Result of the spherical harmonics as complex number.
 */
std::complex<double> sphericalHarmonics(int m, int n, double theta, double phi);

/**
 * If multiple spherical harmonics are evaluated for the same m and theta in a row. It is faster to first call
 * this function and then use sphericalHarmonicsCached.
 * @param m (superscript)
 * @param n maximum n (subscript) that will be used.
 * @param theta (first parameter)
 */
void sphericalHarmonicsBuildCache(int m, int n, double theta, double phi);
/**
 * If multiple spherical harmonics are evaluated for the same m and theta in a row. It is faster to first call
 * sphericalHarmonics and then use this function.
 * @param m (superscript)
 * @param n (subscript)
 * @param theta (first parameter)
 * @param phi (second parameter)
 * @return Result of the spherical harmonics as complex number.
 */
inline std::complex<double> sphericalHarmonicsCached(int m, int n, double theta, double phi) {
  return MathOpt::sphericalCache * std::sqrt(factorial(n - std::abs(m)) / factorial(n + std::abs(m))) *
         getLegendreCache(n);
}

/**
 * Returns A(m,n) as defined in 5.23.
 * @param m (superscript)
 * @param n (subscript)
 * @return A(m,n)
 */
inline double getA(int m, int n) { return MathOpt::getAValue[maxFactorialParameter / 2 + m][n]; }

/**
 * Returns a copy of the 3d vector 'point' using 'center' as origin.
 * @param The 3d vector in cartesian coordinates to convert.
 * @param The 3d vector in cartesian coordinates used as origin.
 * @return The converted 3d vector in cartesian coordinates.
 */
std::array<double, 3> subtract(const std::array<double, 3> &a, const std::array<double, 3> &b);
std::array<double, 3> add(const std::array<double, 3> &a, const std::array<double, 3> &b);
std::array<double, 3> mul(const std::array<double, 3> &a, double scalar);

/**
 * Returns std::pow(1i, exponent), where 1i is the imaginary unit.
 * @param exponent Exponent as integer.
 * @return
 */
inline std::complex<double> powI(int exponent) { return MathOpt::imaginaryPower[exponent % 4 + 4]; }

}