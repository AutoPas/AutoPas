/**
 * @file Math3D.h
 * @date 23.09.19
 * @author Joachim Marin
 */

#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <complex>


/**
 * Returns a copy of a 3d vector in cartesian coordinates using spherical coordinates.
 * @param 3d vector in cartesian coordinates
 * @return copied 3d vector in spherical coordinates
 */
std::array<double, 3> toSpherical(const std::array<double, 3> &cartesian) {
  double x = cartesian[0];
  double y = cartesian[1];
  double z = cartesian[2];
  double r = std::sqrt(x * x + y * y + z * z);
  double alpha = 0;
  double beta = 0;
  if (x != 0.0) {
    alpha = std::atan(y / x);
  } else {
    if (y > 0) {
      alpha = M_PI_2;
    } else {
      alpha = -M_PI_2;
    }
  }
  if (r > 0) {
    beta = std::acos(z / r);
  }

  if (__isnan(alpha)) {
    std::cerr << "alpha is nan" << std::endl;
  }

  if (__isnan(beta)) {
    std::cerr << "beta is nan" << std::endl;
  }

  return std::array<double, 3>({r, alpha, beta});
}

/*std::array<double, 3> toCartesian(std::array<double, 3> spherical) {
  double r = spherical[0];
  double alpha = spherical[1];
  double beta = spherical[2];
  double x = r * std::sin(beta) * std::cos(alpha);
  double y = r * std::sin(beta) * std::sin(alpha);
  double z = r * std::cos(beta);
  return std::array<double, 3>({x, y, z});
}*/

int doubleFactorialRec(int x) {
  return (x <= 1) ? 1 : doubleFactorialRec(x - 2);
}

// Parameter is only checked at the first call.
int doubleFactorial(int x) {
  if (x < 0) {
    std::cerr << "doubleFactorial(" << x << ") is not defined" << std::endl;
  }
  return doubleFactorialRec(x);
}

int factorialRec(int x) {
  return (x <= 1) ? 1 : factorialRec(x - 1);
}

// Parameter is only checked at the first call.
int factorial(int x) {
  if (x < 0) {
    std::cerr << "factorial(" << x << ") is not defined" << std::endl;
  }
  return factorialRec(x);
}

// Recurrence relation from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2634295/#APP1
double associatedLegendrePolynomialRec(int m, int n, double x) {

  double ret = 1;
  if (n == m) {
    // special case n=m=0
    if (n == 0) {
      return ret;
    }

    if (n % 2 == 1) {
      ret = -1;
    }
    ret *= doubleFactorial(2 * m - 1);
    ret *= std::pow(1.0 - x * x, 0.5 * m);
    return ret;
  }
  if (n == m - 1) {
    return x * (2 * m + 1) * associatedLegendrePolynomialRec(m, m, x);
  }

  ret = x * (2 * n - 1) * associatedLegendrePolynomialRec(m, n - 1, x) -
        (n + m - 1) * associatedLegendrePolynomialRec(m, n - 2, x);
  ret /= (n - m);
  return ret;
}


double associatedLegendrePolynomial(int m, int n, double x) {
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
    return std::pow(-1, m) * factorial(n - m) / factorial(n + m) * associatedLegendrePolynomial(-m, n, x);
  }
  return associatedLegendrePolynomialRec(m, n, x);
}

/**
 * Calculates the spherical harmonics (Y) using associated legendre polynomials.
 * @param m (superscript)
 * @param n (subscript)
 * @param alpha (first parameter)
 * @param beta (second parameter)
 * @return Result of the spherical harmonics as complex number.
 */
std::complex<double> sphericalHarmonics(int m, int n, double alpha, double beta) {

  if (n < std::abs(m)) {
    std::cerr << "sphericalHarmonics(" << m << "," << n << "," << alpha << "," << beta << ") is not defined for n < |m|"
              << std::endl;
  }

  std::complex<double> ret = 1;
  double root = std::sqrt(factorial(n - std::abs(m)) / factorial(n + std::abs(m)));

  ret *= root;

  using namespace std::complex_literals;
  ret *= std::exp(1i * double(m) * beta);

  ret *= associatedLegendrePolynomial(std::abs(m), n, std::cos(alpha));

  return ret;
}

/**
 * Returns a copy of the 3d vector 'point' using 'center' as origin.
 * @param The 3d vector in cartesian coordinates to convert.
 * @param The 3d vector in cartesian coordinates used as origin.
 * @return The converted 3d vector in cartesian coordinates.
 */
std::array<double, 3> center(const std::array<double, 3> &point, const std::array<double, 3> &center) {
  return std::array<double, 3>({point[0] - center[0], point[1] - center[1], point[2] - center[2]});
}