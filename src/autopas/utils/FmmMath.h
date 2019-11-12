/**
 * @file FmmMath.h
 * @date 12.11.19
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

constexpr auto maxFactorialParameter = 128;

/**
 * Must be called before using functions from Math3D.h.
 */

template <typename floatType, typename intType>
class FmmMath {
  using Complex = std::complex<floatType>;

 public:
  FmmMath() { FmmMath<floatType, intType>::initialize(); }

  void static initialize() {
    if (FmmMath::initialized) {
      return;
    }

    // factorial
    FmmMath::factorialValue = std::vector<floatType>(maxFactorialParameter + 1);
    FmmMath::factorialValue[0] = 1;
    for (intType i = 1; i <= maxFactorialParameter; ++i) {
      FmmMath::factorialValue[i] = i * FmmMath::factorialValue[i - 1];
    }

    // double factorial
    FmmMath::doubleFactorialValue = std::vector<floatType>(maxFactorialParameter + 1);
    FmmMath::doubleFactorialValue[0] = 1;
    FmmMath::doubleFactorialValue[1] = 1;
    for (intType i = 2; i <= maxFactorialParameter; ++i) {
      FmmMath::doubleFactorialValue[i] = i * FmmMath::doubleFactorialValue[i - 2];
    }

    // getA
    FmmMath::getAValue = std::vector<std::vector<floatType>>(maxFactorialParameter + 1,
                                                             std::vector<floatType>(maxFactorialParameter / 2 + 1));
    for (intType m = -maxFactorialParameter / 2; m <= maxFactorialParameter / 2; ++m) {
      for (intType n = std::abs(m); n <= maxFactorialParameter / 2; ++n) {
        FmmMath::getAValue[maxFactorialParameter / 2 + m][n] = FmmMath::calculateA(m, n);
      }
    }

    // powI
    FmmMath::imaginaryPower = std::vector<Complex>(8);
    using namespace std::complex_literals;
    FmmMath::imaginaryPower[0] = 1;
    FmmMath::imaginaryPower[1] = 1i;
    FmmMath::imaginaryPower[2] = -1;
    FmmMath::imaginaryPower[3] = -1i;
    for (intType i = 0; i < 4; ++i) {
      FmmMath::imaginaryPower[i + 4] = FmmMath::imaginaryPower[i];
    }

    FmmMath::initialized = true;
  }

  /**
   * Returns a copy of a 3d vector in cartesian coordinates using spherical coordinates.
   * @param 3d vector in cartesian coordinates
   * @return copied 3d vector in spherical coordinates
   */
  std::array<floatType, 3> static toSpherical(const std::array<floatType, 3> &cartesian) {
    floatType x = cartesian[0];
    floatType y = cartesian[1];
    floatType z = cartesian[2];
    floatType theta;

    floatType rho = std::sqrt(x * x + y * y + z * z);
    floatType phi = std::atan2(y, x);
    if (rho > 0) {
      theta = std::acos(z / rho);
    } else {
      theta = 0;
    }

    if (__isnan(phi)) {
      std::cerr << "phi is nan" << std::endl;
    }

    if (__isnan(theta)) {
      std::cerr << "theta is nan" << std::endl;
    }

    return std::array<floatType, 3>({rho, theta, phi});
  }

  inline floatType static doubleFactorial(intType x) { return doubleFactorialValue.at(x); }

  inline floatType static factorial(intType x) { return factorialValue.at(x); }

  /**
   * Calculates the associated Legendre polynomial of order m and degree n at x.
   * @param m
   * @param n
   * @param x
   * @return
   */
  floatType associatedLegendrePolynomial(intType m, intType n, floatType x) {
    if (m > n) {
      std::cerr << "associatedLegendrePolynomial(" << m << "," << n << "," << x << ") is not defined for m > n"
                << std::endl;
    }
    if (x > 1 || x < -1) {
      std::cerr << "associatedLegendrePolynomial(" << m << "," << n << "," << x << ") is only defined for -1 <= x <= 1"
                << std::endl;
    }

    if (n < 0) {
      std::cerr << "associatedLegendrePolynomial(" << m << "," << n << "," << x << ") is not defined for n < 0"
                << std::endl;
    }

    // Doing the recurrence relation for m < 0 at the start ensures that m >= 0 afterwards.
    // Here this function is used, to check for the updated parameters again.
    // Then associatedLegendrePolynomial will only ever be called with correct parameters.
    if (m < 0) {
      return std::pow(-1, m) * FmmMath::factorial(n - m) / FmmMath::factorial(n + m) *
             associatedLegendrePolynomial(-m, n, x);
    }
    // build cache
    // The recurrence relation accesses associated legendre polynomials for degrees n-1 and n-2.
    // To ensure every degree is only calculated once, the degrees are calculated here from m to n.

    legendreCache = std::vector<floatType>(n - m + 2);

    legendreCache[0] = associatedLegendrePolynomialRec(m, m, x);
    legendreCache[1] = associatedLegendrePolynomialRec(m, m + 1, x);

    for (intType i = 2; i <= n - m; ++i) {
      intType cacheN = m + i;
      legendreCache[i] =
          (x * (2 * cacheN - 1) * legendreCache[i - 1] - (cacheN + m - 1) * legendreCache[i - 2]) / (cacheN - m);
    }
    legendreLastM = m;

    // return associatedLegendrePolynomialRec(m, n, x);
    return legendreCache[n - m];
  }

  /**
   * After associatedLegendrePolynomial has been called for a certain m and x. This function can be used to get the
   * value of the associated Legendre polynomial of same order m and same value x, but smaller n.
   * @param n Degree of the associated legendre polynomial.
   * @return
   */
  inline floatType getLegendreCache(intType n) { return legendreCache[n - legendreLastM]; }

  /**
   * Calculates the spherical harmonics (Y) using associated legendre polynomials.
   * @param m (superscript)
   * @param n (subscript)
   * @param theta (first parameter)
   * @param phi (second parameter)
   * @return Result of the spherical harmonics as complex number.
   */
  Complex sphericalHarmonics(intType m, intType n, floatType theta, floatType phi) {
    if (n < std::abs(m)) {
      std::cerr << "sphericalHarmonics(" << m << "," << n << "," << theta << "," << phi
                << ") is not defined for n < |m|" << std::endl;
    }

    using namespace std::complex_literals;
    return std::exp(1i * phi * m) *
           std::sqrt(FmmMath::factorial(n - std::abs(m)) / FmmMath::factorial(n + std::abs(m))) *
           associatedLegendrePolynomial(std::abs(m), n, std::cos(theta));
  }

  /**
   * If multiple spherical harmonics are evaluated for the same m and theta in a row. It is faster to first call
   * this function and then use sphericalHarmonicsCached.
   * @param m (superscript)
   * @param n maximum n (subscript) that will be used.
   * @param theta (first parameter)
   */
  void sphericalHarmonicsBuildCache(intType m, intType n, floatType theta, floatType phi) {
    associatedLegendrePolynomial(std::abs(m), n, std::cos(theta));
    using namespace std::complex_literals;
    sphericalCache = std::exp(1i * static_cast<floatType>(m) * phi);
  }
  /**
   * If multiple spherical harmonics are evaluated for the same m and theta in a row. It is faster to first call
   * sphericalHarmonics and then use this function.
   * @param m (superscript)
   * @param n (subscript)
   * @param theta (first parameter)
   * @param phi (second parameter)
   * @return Result of the spherical harmonics as complex number.
   */
  inline Complex sphericalHarmonicsCached(intType m, intType n, floatType theta, floatType phi) {
    return sphericalCache * std::sqrt(factorial(n - std::abs(m)) / factorial(n + std::abs(m))) * getLegendreCache(n);
  }

  /**
   * Returns A(m,n) as defined in 5.23.
   * @param m (superscript)
   * @param n (subscript)
   * @return A(m,n)
   */
  static inline floatType getA(intType m, intType n) { return getAValue[maxFactorialParameter / 2 + m][n]; }

  /**
   * Returns std::pow(1i, exponent), where 1i is the imaginary unit.
   * @param exponent Exponent as integer.
   * @return
   */
  static inline Complex powI(intType exponent) { return imaginaryPower[exponent % 4 + 4]; }

 private:
  // Recurrence relation from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2634295/#APP1
  floatType associatedLegendrePolynomialRec(intType m, intType n, floatType x) {
    floatType ret = 1.0;
    if (n == m) {
      // special case n=m=0
      if (n == 0) {
        return ret;
      }

      if (n % 2 == 1) {
        ret = -1;
      }
      ret *= FmmMath<floatType, intType>::doubleFactorial(2 * m - 1);
      ret *= std::pow(1.0 - x * x, 0.5 * m);
      return ret;
    }
    if (n == m + 1) {
      return x * (2 * m + 1) * associatedLegendrePolynomialRec(m, m, x);
    }

    ret = x * (2 * n - 1) * associatedLegendrePolynomialRec(m, n - 1, x) -
          (n + m - 1) * associatedLegendrePolynomialRec(m, n - 2, x);
    ret /= (n - m);
    return ret;
  }

  static floatType calculateA(intType m, intType n) {
    if (n - m < 0 || n + m < 0) {
      std::cerr << "getA(" << m << "," << n << ") is not defined for n - m < 0 or n + m < 0" << std::endl;
    }

    floatType result = std::pow(-1, n) / std::sqrt(FmmMath<floatType, intType>::factorial(n - m) *
                                                   FmmMath<floatType, intType>::factorial(n + m));

    assert(!__isnan(result));

    return result;
  }

  static bool initialized;
  static std::vector<floatType> factorialValue;
  static std::vector<floatType> doubleFactorialValue;
  static std::vector<std::vector<floatType>> getAValue;
  static std::vector<Complex> imaginaryPower;

  intType legendreLastM = 0;
  Complex sphericalCache;
  std::vector<floatType> legendreCache;
};

template <typename floatType, typename intType>
std::vector<floatType> FmmMath<floatType, intType>::factorialValue;

template <typename floatType, typename intType>
std::vector<floatType> FmmMath<floatType, intType>::doubleFactorialValue;

template <typename floatType, typename intType>
std::vector<std::vector<floatType>> FmmMath<floatType, intType>::getAValue;

template <typename floatType, typename intType>
std::vector<std::complex<floatType>> FmmMath<floatType, intType>::imaginaryPower;

template <typename floatType, typename intType>
bool FmmMath<floatType, intType>::initialized = false;
