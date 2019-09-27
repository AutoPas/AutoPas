/**
 * @file Math3D.h
 * @date 23.09.19
 * @author Joachim Marin
 */

#pragma once

#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <vector>

void initMath();

/**
 * Returns a copy of a 3d vector in cartesian coordinates using spherical coordinates.
 * @param 3d vector in cartesian coordinates
 * @return copied 3d vector in spherical coordinates
 */
std::array<double, 3> toSpherical(const std::array<double, 3> &cartesian);

double doubleFactorial(int x);

double factorial(int x);

double associatedLegendrePolynomial(int m, int n, double x);

/**
 * Calculates the spherical harmonics (Y) using associated legendre polynomials.
 * @param m (superscript)
 * @param n (subscript)
 * @param theta (first parameter)
 * @param phi (second parameter)
 * @return Result of the spherical harmonics as complex number.
 */
std::complex<double> sphericalHarmonics(int m, int n, double theta, double phi);

double getA(int m, int n);

/**
 * Returns a copy of the 3d vector 'point' using 'center' as origin.
 * @param The 3d vector in cartesian coordinates to convert.
 * @param The 3d vector in cartesian coordinates used as origin.
 * @return The converted 3d vector in cartesian coordinates.
 */
std::array<double, 3> subtract(const std::array<double, 3> &a, const std::array<double, 3> &b);