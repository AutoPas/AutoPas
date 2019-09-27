/**
 * @file Math3D.cpp
 * @date 23.09.19
 * @author Joachim Marin
 */

#include "Math3D.h"

class MathOpt {
 public:
  static std::vector<double> factorialValue;
  static std::vector<double> doubleFactorialValue;
};
std::vector<double> MathOpt::factorialValue;
std::vector<double> MathOpt::doubleFactorialValue;

std::array<double, 3> toSpherical(const std::array<double, 3> &cartesian) {
  double x = cartesian[0];
  double y = cartesian[1];
  double z = cartesian[2];
  double rho;
  double phi;
  double theta;

  rho = std::sqrt(x * x + y * y + z * z);
  phi = std::atan2(y, x);
  if (rho > 0) {
    theta = std::acos(z / rho);
  } else {
    theta = 0;
  }

  /*double oldPhi = 0;
  if (x != 0.0) {
    oldPhi = std::atan(y / x);
  } else {
    if (y > 0) {
      oldPhi = M_PI_2;
    } else {
      oldPhi = -M_PI_2;
    }
  }

  if (phi != oldPhi) {
    std::cerr << "toSpherical(" << x << "," << y << "," << z << ") phi = " << phi << ", oldPhi = " << oldPhi
              << std::endl;
  }*/

  if (__isnan(phi)) {
    std::cerr << "phi is nan" << std::endl;
  }

  if (__isnan(theta)) {
    std::cerr << "theta is nan" << std::endl;
  }

  /*double r = std::sqrt(x * x + y * y + z * z);
  double phi = 0;
  double theta = 0;
  if (x != 0.0) {
    phi = std::atan(y / x);
  } else {
    if (y > 0) {
      phi = M_PI_2;
    } else {
      phi = -M_PI_2;
    }
  }
  if (r > 0) {
    theta = std::acos(z / r);
  }

  if (__isnan(phi)) {
    std::cerr << "phi is nan" << std::endl;
  }

  if (__isnan(theta)) {
    std::cerr << "theta is nan" << std::endl;
  }*/

  return std::array<double, 3>({rho, theta /*alpha*/, phi /*beta*/});
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

// Parameter is only checked at the first call.
double doubleFactorial(int x) {
  if (x < 0) {
    std::cerr << "doubleFactorial(" << x << ") is not defined" << std::endl;
  }
  return MathOpt::doubleFactorialValue.at(x);
}

// Parameter is only checked at the first call.
double factorial(int x) {
  if (x < 0) {
    std::cerr << "factorial(" << x << ") is not defined" << std::endl;
  }
  return MathOpt::factorialValue.at(x);
  // return factorialRec(x);
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
  std::vector<double> cache(n - m + 2);
  cache[0] = associatedLegendrePolynomialRec(m, m, x);
  cache[1] = associatedLegendrePolynomialRec(m, m + 1, x);

  for (int i = 2; i <= n - m; ++i) {
    int cacheN = m + i;
    cache[i] = (x * (2 * cacheN - 1) * cache[i - 1] - (cacheN + m - 1) * cache[i - 2]) / (cacheN - m);
  }

  // return associatedLegendrePolynomialRec(m, n, x);
  return cache[n - m];
}

std::complex<double> sphericalHarmonics(int m, int n, double theta, double phi) {
  if (n < std::abs(m)) {
    std::cerr << "sphericalHarmonics(" << m << "," << n << "," << theta << "," << phi << ") is not defined for n < |m|"
              << std::endl;
  }

  std::complex<double> ret = 1;
  double root = std::sqrt(factorial(n - std::abs(m)) / factorial(n + std::abs(m)));

  ret *= root;

  using namespace std::complex_literals;
  ret *= std::exp(1i * double(m) * phi);

  ret *= associatedLegendrePolynomial(std::abs(m), n, std::cos(theta));

  return ret;
}

double getA(int m, int n) {
  if (n - m < 0 || n + m < 0) {
    std::cerr << "getA(" << m << "," << n << ") is not defined for n - m < 0 or n + m < 0" << std::endl;
  }

  double result = std::pow(-1, n) / std::sqrt(factorial(n - m) * factorial(n + m));

  /*std::cout << std::pow(-1, n) << "/"
            << "sqrt(" << factorial(n - m) << "*" << factorial(n + m) << ")" << std::endl;*/

  assert(!__isnan(result));

  return result;
}

std::array<double, 3> subtract(const std::array<double, 3> &a, const std::array<double, 3> &b) {
  return std::array<double, 3>({a[0] - b[0], a[1] - b[1], a[2] - b[2]});
}

void initMath() {
  // factorial
  MathOpt::factorialValue = std::vector<double>(100);
  MathOpt::factorialValue[0] = 1;
  for (int i = 1; i < 100; ++i) {
    MathOpt::factorialValue[i] = i * MathOpt::factorialValue[i - 1];
  }

  // double factorial
  MathOpt::doubleFactorialValue = std::vector<double>(100);
  MathOpt::doubleFactorialValue[0] = 1;
  MathOpt::doubleFactorialValue[1] = 1;
  for (int i = 2; i < 100; ++i) {
    MathOpt::doubleFactorialValue[i] = i * MathOpt::doubleFactorialValue[i - 2];
  }
}