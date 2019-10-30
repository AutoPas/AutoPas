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

using Complex = std::complex<double>;

constexpr auto maxFactorialParameter = 128;

/**
 * Must be called before using functions from Math3D.h.
 */

class Math3D {
 public:
  Math3D() {
    Math3D::initialize();
  }

  void static initialize();

  /**
   * Returns a copy of a 3d vector in cartesian coordinates using spherical coordinates.
   * @param 3d vector in cartesian coordinates
   * @return copied 3d vector in spherical coordinates
   */
  std::array<double, 3> static toSpherical(const std::array<double, 3> &cartesian);

  inline double static doubleFactorial(int x) { return doubleFactorialValue.at(x); }

  inline double static factorial(int x) { return factorialValue.at(x); }

  /**
   * Calculates the associated Legendre polynomial of order m and degree n at x.
   * @param m
   * @param n
   * @param x
   * @return
   */
  double associatedLegendrePolynomial(int m, int n, double x);

  /**
   * After associatedLegendrePolynomial has been called for a certain m and x. This function can be used to get the
   * value of the associated Legendre polynomial of same order m and same value x, but smaller n.
   * @param n Degree of the associated legendre polynomial.
   * @return
   */
  inline double getLegendreCache(int n) { return legendreCache[n - legendreLastM]; }

  /**
   * Calculates the spherical harmonics (Y) using associated legendre polynomials.
   * @param m (superscript)
   * @param n (subscript)
   * @param theta (first parameter)
   * @param phi (second parameter)
   * @return Result of the spherical harmonics as complex number.
   */
  Complex sphericalHarmonics(int m, int n, double theta, double phi);

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
  inline Complex sphericalHarmonicsCached(int m, int n, double theta, double phi) {
    return sphericalCache * std::sqrt(factorial(n - std::abs(m)) / factorial(n + std::abs(m))) * getLegendreCache(n);
  }

  /**
   * Returns A(m,n) as defined in 5.23.
   * @param m (superscript)
   * @param n (subscript)
   * @return A(m,n)
   */
  static inline double getA(int m, int n) { return getAValue[maxFactorialParameter / 2 + m][n]; }

  /**
   * Returns std::pow(1i, exponent), where 1i is the imaginary unit.
   * @param exponent Exponent as integer.
   * @return
   */
  static inline Complex powI(int exponent) { return imaginaryPower[exponent % 4 + 4]; }

 private:
  static bool initialized;
  static std::vector<double> factorialValue;
  static std::vector<double> doubleFactorialValue;
  static std::vector<std::vector<double>> getAValue;
  static std::vector<Complex> imaginaryPower;

  int legendreLastM = 0;
  std::complex<double> sphericalCache;
  std::vector<double> legendreCache;
};
