/**
 * @file Math3D.cpp
 * @date 23.09.19
 * @author Joachim Marin
 */

#include "Math3D.h"

namespace Math3D {

std::vector<double> MathOpt::factorialValue;
std::vector<double> MathOpt::doubleFactorialValue;
std::vector<double> MathOpt::legendreCache;
std::complex<double> MathOpt::sphericalCache;
std::vector<std::vector<double>> MathOpt::getAValue;
int MathOpt::legendreLastM;
std::vector<std::complex<double>> MathOpt::imaginaryPower;

std::array<double, 3> toSpherical(const std::array<double, 3> &cartesian) {
  double x = cartesian[0];
  double y = cartesian[1];
  double z = cartesian[2];
  double theta;

  double rho = std::sqrt(x * x + y * y + z * z);
  double phi = std::atan2(y, x);
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

  return std::array<double, 3>({rho, theta /*alpha*/, phi /*beta*/});
}

// Recurrence relation from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2634295/#APP1
double associatedLegendrePolynomialRec(int m, int n, double x) {
  double ret = 1.0;
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
  if (n == m + 1) {
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
  // build cache
  // The recurrence relation accesses associated legendre polynomials for degrees n-1 and n-2.
  // To ensure every degree is only calculated once, the degrees are calculated here from m to n.

  MathOpt::legendreCache = std::vector<double>(n - m + 2);

  MathOpt::legendreCache[0] = associatedLegendrePolynomialRec(m, m, x);
  MathOpt::legendreCache[1] = associatedLegendrePolynomialRec(m, m + 1, x);

  for (int i = 2; i <= n - m; ++i) {
    int cacheN = m + i;
    MathOpt::legendreCache[i] =
        (x * (2 * cacheN - 1) * MathOpt::legendreCache[i - 1] - (cacheN + m - 1) * MathOpt::legendreCache[i - 2]) /
        (cacheN - m);
  }
  MathOpt::legendreLastM = m;

  // return associatedLegendrePolynomialRec(m, n, x);
  return MathOpt::legendreCache[n - m];
}

std::complex<double> sphericalHarmonics(int m, int n, double theta, double phi) {
  if (n < std::abs(m)) {
    std::cerr << "sphericalHarmonics(" << m << "," << n << "," << theta << "," << phi << ") is not defined for n < |m|"
              << std::endl;
  }

  using namespace std::complex_literals;
  return std::exp(1i * double(m) * phi) * std::sqrt(factorial(n - std::abs(m)) / factorial(n + std::abs(m))) *
         associatedLegendrePolynomial(std::abs(m), n, std::cos(theta));
}

void sphericalHarmonicsBuildCache(int m, int n, double theta, double phi) {
  associatedLegendrePolynomial(std::abs(m), n, std::cos(theta));
  using namespace std::complex_literals;
  MathOpt::sphericalCache = std::exp(1i * static_cast<double>(m) * phi);
}

double calculateA(int m, int n) {
  if (n - m < 0 || n + m < 0) {
    std::cerr << "getA(" << m << "," << n << ") is not defined for n - m < 0 or n + m < 0" << std::endl;
  }

  double result = std::pow(-1, n) / std::sqrt(factorial(n - m) * factorial(n + m));

  assert(!__isnan(result));

  return result;
}

std::array<double, 3> subtract(const std::array<double, 3> &a, const std::array<double, 3> &b) {
  return std::array<double, 3>({a[0] - b[0], a[1] - b[1], a[2] - b[2]});
}
std::array<double, 3> add(const std::array<double, 3> &a, const std::array<double, 3> &b) {
  return std::array<double, 3>({a[0] + b[0], a[1] + b[1], a[2] + b[2]});
};
std::array<double, 3> mul(const std::array<double, 3> &a, double scalar) {
  return std::array<double, 3>({a[0] * scalar, a[1] * scalar, a[2] * scalar});
};

void initMath() {
  // factorial
  MathOpt::factorialValue = std::vector<double>(maxFactorialParameter + 1);
  MathOpt::factorialValue[0] = 1;
  for (int i = 1; i <= maxFactorialParameter; ++i) {
    MathOpt::factorialValue[i] = i * MathOpt::factorialValue[i - 1];
  }

  // double factorial
  MathOpt::doubleFactorialValue = std::vector<double>(maxFactorialParameter + 1);
  MathOpt::doubleFactorialValue[0] = 1;
  MathOpt::doubleFactorialValue[1] = 1;
  for (int i = 2; i <= maxFactorialParameter; ++i) {
    MathOpt::doubleFactorialValue[i] = i * MathOpt::doubleFactorialValue[i - 2];
  }

  // getA
  MathOpt::getAValue =
      std::vector<std::vector<double>>(maxFactorialParameter + 1, std::vector<double>(maxFactorialParameter / 2 + 1));
  for (int m = -maxFactorialParameter / 2; m <= maxFactorialParameter / 2; ++m) {
    for (int n = std::abs(m); n <= maxFactorialParameter / 2; ++n) {
      MathOpt::getAValue[maxFactorialParameter / 2 + m][n] = calculateA(m, n);
    }
  }

  // powI
  MathOpt::imaginaryPower = std::vector<std::complex<double>>(8);
  using namespace std::complex_literals;
  MathOpt::imaginaryPower[0] = 1;
  MathOpt::imaginaryPower[1] = 1i;
  MathOpt::imaginaryPower[2] = -1;
  MathOpt::imaginaryPower[3] = -1i;
  for (int i = 0; i < 4; ++i) {
    MathOpt::imaginaryPower[i + 4] = MathOpt::imaginaryPower[i];
  }
}
}