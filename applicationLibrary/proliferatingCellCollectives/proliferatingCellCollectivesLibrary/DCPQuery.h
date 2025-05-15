
/**
 * @file SpherocylinderCell.cpp
 * @author Manuel Lerchner (Source: https://github.com/flatironinstitute/SimToolbox/blob/main/Collision/DCPQuery.hpp)
 * @date 14/05/2025
 */

#include <cmath>
#include <limits>

#include "autopas/utils/ArrayMath.h"

/**
 * @brief functor for minimal segment-segment distance query in 3D space
 *
 * this object must be thread private in multi-threading environment
 * @tparam N dimension
 * @tparam Real floating-point type
 * @tparam Vector 3D vector type
 */
class DCPQuery {
 public:
  struct Result {
    double distance, sqrDistance;
    double parameter[2];
    std::array<double, 3> closest[2];
  };

  /**
   * @brief functor for computing minimal segment-segment distance in 3D space
   *
   * @param P0 end point 0 of segment P
   * @param P1 end point 1 of segment P
   * @param Q0 end point 0 of segment Q
   * @param Q1 end point 1 of segment Q
   * @param Ploc result point on P
   * @param Qloc result point on Q
   * @param s \f$s\in[0,1]\f$ describing Ploc on P
   * @param t \f$t\in[0,1]\f$ describing Qloc on Q
   * @return Real computed minimal distance
   */
  double operator()(std::array<double, 3> const &P0, std::array<double, 3> const &P1, std::array<double, 3> const &Q0,
                    std::array<double, 3> const &Q1, std::array<double, 3> &Ploc, std::array<double, 3> &Qloc,
                    double &s, double &t);

  double operator()(std::array<double, 3> const &P0, std::array<double, 3> const &P1, std::array<double, 3> const &Q0,
                    std::array<double, 3> const &Q1, std::array<double, 3> &Ploc, std::array<double, 3> &Qloc);

 private:
  // Compute the root of h(z) = h0 + slope*z and clamp it to the interval
  // [0,1].  It is required that for h1 = h(1), either (h0 < 0 and h1 > 0)
  // or (h0 > 0 and h1 < 0).
  double GetClampedRoot(double slope, double h0, double h1);

  // Compute the intersection of the line dR/ds = 0 with the domain [0,1]^2.
  // The direction of the line dR/ds is conjugate to (1,0), so the algorithm
  // for minimization is effectively the conjugate gradient algorithm for a
  // quadratic function.
  void ComputeIntersection(double const sValue[2], int const classify[2], int edge[2], double end[2][2]);

  // Compute the location of the minimum of R on the segment of intersection
  // for the line dR/ds = 0 and the domain [0,1]^2.
  void ComputeMinimumParameters(int const edge[2], double const end[2][2], double parameter[2]);

  // The coefficients of R(s,t), not including the constant term.
  double mA, mB, mC, mD, mE;

  // dR/ds(i,j) at the four corners of the domain
  double mF00, mF10, mF01, mF11;

  // dR/dt(i,j) at the four corners of the domain
  double mG00, mG10, mG01, mG11;
};

double DCPQuery::operator()(std::array<double, 3> const &P0, std::array<double, 3> const &P1,
                            std::array<double, 3> const &Q0, std::array<double, 3> const &Q1,
                            std::array<double, 3> &Ploc, std::array<double, 3> &Qloc) {
  double s, t = 0;
  return (*this)(P0, P1, Q0, Q1, Ploc, Qloc, s, t);
}

double DCPQuery::operator()(std::array<double, 3> const &P0, std::array<double, 3> const &P1,
                            std::array<double, 3> const &Q0, std::array<double, 3> const &Q1,
                            std::array<double, 3> &Ploc, std::array<double, 3> &Qloc, double &s, double &t) {
  using namespace autopas::utils::ArrayMath;

  Result result;

  // The code allows degenerate line segments; that is, P0 and P1 can be
  // the same point or Q0 and Q1 can be the same point.  The quadratic
  // function for squared distance between the segment is
  //   R(s,t) = a*s^2 - 2*b*s*t + c*t^2 + 2*d*s - 2*e*t + f
  // for (s,t) in [0,1]^2 where
  //   a = Dot(P1-P0,P1-P0), b = Dot(P1-P0,Q1-Q0), c = Dot(Q1-Q0,Q1-Q0),
  //   d = Dot(P1-P0,P0-Q0), e = Dot(Q1-Q0,P0-Q0), f = Dot(P0-Q0,P0-Q0)
  auto P1mP0 = P1 - P0;
  auto Q1mQ0 = Q1 - Q0;
  auto P0mQ0 = P0 - Q0;
  mA = dot(P1mP0, P1mP0);
  mB = dot(P1mP0, Q1mQ0);
  mC = dot(Q1mQ0, Q1mQ0);
  mD = dot(P1mP0, P0mQ0);
  mE = dot(Q1mQ0, P0mQ0);

  mF00 = mD;
  mF10 = mF00 + mA;
  mF01 = mF00 - mB;
  mF11 = mF10 - mB;

  mG00 = -mE;
  mG10 = mG00 - mB;
  mG01 = mG00 + mC;
  mG11 = mG10 + mC;

  if (mA > 0 && mC > 0) {
    // Compute the solutions to dR/ds(s0,0) = 0 and dR/ds(s1,1) = 0.  The
    // location of sI on the s-axis is stored in classifyI (I = 0 or 1).  If
    // sI <= 0, classifyI is -1.  If sI >= 1, classifyI is 1.  If 0 < sI < 1,
    // classifyI is 0.  This information helps determine where to search for
    // the minimum point (s,t).  The fij values are dR/ds(i,j) for i and j in
    // {0,1}.

    double sValue[2];
    sValue[0] = GetClampedRoot(mA, mF00, mF10);
    sValue[1] = GetClampedRoot(mA, mF01, mF11);

    int classify[2];
    for (int i = 0; i < 2; ++i) {
      if (sValue[i] <= 0) {
        classify[i] = -1;
      } else if (sValue[i] >= 1) {
        classify[i] = +1;
      } else {
        classify[i] = 0;
      }
    }

    if (classify[0] == -1 && classify[1] == -1) {
      // The minimum must occur on s = 0 for 0 <= t <= 1.
      result.parameter[0] = 0;
      result.parameter[1] = GetClampedRoot(mC, mG00, mG01);
    } else if (classify[0] == +1 && classify[1] == +1) {
      // The minimum must occur on s = 1 for 0 <= t <= 1.
      result.parameter[0] = 1;
      result.parameter[1] = GetClampedRoot(mC, mG10, mG11);
    } else {
      // The line dR/ds = 0 intersects the domain [0,1]^2 in a
      // nondegenerate segment.  Compute the endpoints of that segment,
      // end[0] and end[1].  The edge[i] flag tells you on which domain
      // edge end[i] lives: 0 (s=0), 1 (s=1), 2 (t=0), 3 (t=1).
      int edge[2];
      double end[2][2];
      ComputeIntersection(sValue, classify, edge, end);

      // The directional derivative of R along the segment of
      // intersection is
      //   H(z) = (end[1][1]-end[1][0])*dR/dt((1-z)*end[0] + z*end[1])
      // for z in [0,1].  The formula uses the fact that dR/ds = 0 on
      // the segment.  Compute the minimum of H on [0,1].
      ComputeMinimumParameters(edge, end, result.parameter);
    }
  } else {
    if (mA > 0) {
      // The Q-segment is degenerate (Q0 and Q1 are the same point) and
      // the quadratic is R(s,0) = a*s^2 + 2*d*s + f and has (half)
      // first derivative F(t) = a*s + d.  The closest P-point is
      // interior to the P-segment when F(0) < 0 and F(1) > 0.
      result.parameter[0] = GetClampedRoot(mA, mF00, mF10);
      result.parameter[1] = 0;
    } else if (mC > 0) {
      // The P-segment is degenerate (P0 and P1 are the same point) and
      // the quadratic is R(0,t) = c*t^2 - 2*e*t + f and has (half)
      // first derivative G(t) = c*t - e.  The closest Q-point is
      // interior to the Q-segment when G(0) < 0 and G(1) > 0.
      result.parameter[0] = 0;
      result.parameter[1] = GetClampedRoot(mC, mG00, mG01);
    } else {
      // P-segment and Q-segment are degenerate.
      result.parameter[0] = 0;
      result.parameter[1] = 0;
    }
  }

  result.closest[0] = P0 * (1 - result.parameter[0]) + P1 * result.parameter[0];
  result.closest[1] = Q0 * (1 - result.parameter[1]) + Q1 * result.parameter[1];
  auto diff = result.closest[0] - result.closest[1];
  result.sqrDistance = dot(diff, diff);
  result.distance = sqrt(result.sqrDistance);
  Ploc = result.closest[0];
  Qloc = result.closest[1];
  s = result.parameter[0];
  t = result.parameter[1];
  return result.distance;
}

inline double DCPQuery::GetClampedRoot(double slope, double h0, double h1) {
  // slope = h1-h0
  // h0 and h1 should have different sign, return the zero point
  constexpr double eps = std::numeric_limits<double>::epsilon();
  if (fabs(slope - (h1 - h0)) >= 10 * eps * std::max(abs(h1), abs(h0))) {
    AutoPasLog(DEBUG, fmt::format("Inconsistent parameters in GetClampedRoot: "
                                  "slope = %f, h1 = %f, h0 = %f",
                                  slope, h1, h0));
    return 0.5; 
  }

  double r;

  if (std::abs(h0) < eps && std::abs(h1) < eps) {
    // tiny slope, h0 \approx h1, distance almost a constant, choose mid point
    r = 0.5;
  } else if (h0 < 0) {
    if (h1 > 0) {
      // r = -h0 / slope; // need better accuracy
      // clamp r between [0,1]
      r = std::min(std::max(-h0 / slope, 0.0), 1.0);
    } else {
      r = 1;
    }
  } else {
    r = 0;
  }
  return r;
}

inline void DCPQuery::ComputeIntersection(double const sValue[2], int const classify[2], int edge[2],
                                          double end[2][2]) {
  // The divisions are theoretically numbers in [0,1].  Numerical rounding
  // errors might cause the result to be outside the interval.  When this
  // happens, it must be that both numerator and denominator are nearly
  // zero.  The denominator is nearly zero when the segments are nearly
  // perpendicular.  The numerator is nearly zero when the P-segment is
  // nearly degenerate (mF00 = a is small).  The choice of 0.5 should not
  // cause significant accuracy problems.
  //
  // NOTE:  You can use bisection to recompute the root or even use
  // bisection to compute the root and skip the division.  This is generally
  // slower, which might be a problem for high-performance applications.

  if (classify[0] < 0) {
    edge[0] = 0;
    end[0][0] = 0;
    end[0][1] = mF00 / mB;
    if (end[0][1] < 0 || end[0][1] > 1) {
      end[0][1] = 0.5;
    }

    if (classify[1] == 0) {
      edge[1] = 3;
      end[1][0] = sValue[1];
      end[1][1] = 1;
    } else  // classify[1] > 0
    {
      edge[1] = 1;
      end[1][0] = 1;
      end[1][1] = mF10 / mB;
      if (end[1][1] < 0 || end[1][1] > 1) {
        end[1][1] = 0.5;
      }
    }
  } else if (classify[0] == 0) {
    edge[0] = 2;
    end[0][0] = sValue[0];
    end[0][1] = 0;

    if (classify[1] < 0) {
      edge[1] = 0;
      end[1][0] = 0;
      end[1][1] = mF00 / mB;
      if (end[1][1] < 0 || end[1][1] > 1) {
        end[1][1] = 0.5;
      }
    } else if (classify[1] == 0) {
      edge[1] = 3;
      end[1][0] = sValue[1];
      end[1][1] = 1;
    } else {
      edge[1] = 1;
      end[1][0] = 1;
      end[1][1] = mF10 / mB;
      if (end[1][1] < 0 || end[1][1] > 1) {
        end[1][1] = 0.5;
      }
    }
  } else  // classify[0] > 0
  {
    edge[0] = 1;
    end[0][0] = 1;
    end[0][1] = mF10 / mB;
    if (end[0][1] < 0 || end[0][1] > 1) {
      end[0][1] = 0.5;
    }

    if (classify[1] == 0) {
      edge[1] = 3;
      end[1][0] = sValue[1];
      end[1][1] = 1;
    } else {
      edge[1] = 0;
      end[1][0] = 0;
      end[1][1] = mF00 / mB;
      if (end[1][1] < 0 || end[1][1] > 1) {
        end[1][1] = 0.5;
      }
    }
  }
}

inline void DCPQuery::ComputeMinimumParameters(int const edge[2], double const end[2][2], double parameter[2]) {
  constexpr double eps = std::numeric_limits<double>::epsilon();
  double delta = end[1][1] - end[0][1];
  double h0 = delta * ((-mB * end[0][0] - mE) + mC * end[0][1]);  // source of rounding error
  double h1 = delta * ((-mB * end[1][0] - mE) + mC * end[1][1]);

  if (std::abs(h0) < std::abs(mC) * eps && std::abs(h1) < std::abs(mC) * eps) {
    double z = 0.5;
    double omz = 1 - z;
    parameter[0] = omz * end[0][0] + z * end[1][0];
    parameter[1] = omz * end[0][1] + z * end[1][1];
  } else if (h0 >= 0) {
    if (edge[0] == 0) {
      parameter[0] = 0;
      parameter[1] = GetClampedRoot(mC, mG00, mG01);
    } else if (edge[0] == 1) {
      parameter[0] = 1;
      parameter[1] = GetClampedRoot(mC, mG10, mG11);
    } else {
      parameter[0] = end[0][0];
      parameter[1] = end[0][1];
    }
  } else {
    if (h1 <= 0) {
      if (edge[1] == 0) {
        parameter[0] = 0;
        parameter[1] = GetClampedRoot(mC, mG00, mG01);
      } else if (edge[1] == 1) {
        parameter[0] = 1;
        parameter[1] = GetClampedRoot(mC, mG10, mG11);
      } else {
        parameter[0] = end[1][0];
        parameter[1] = end[1][1];
      }
    } else  // h0 < 0 and h1 > 0
    {
      double z = GetClampedRoot(h1 - h0, h0, h1);
      double omz = 1 - z;
      parameter[0] = omz * end[0][0] + z * end[1][0];
      parameter[1] = omz * end[0][1] + z * end[1][1];
    }
  }

  return;
}