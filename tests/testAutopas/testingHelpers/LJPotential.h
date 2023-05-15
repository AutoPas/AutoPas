/**
 * @file LJPotential.h
 * @author D. Martin
 * @date 13.05.23
 *
 * A simple reference implementation of the lennard jones potential to calculate expected forces
 */

#pragma once
#include "autopas/utils/ConstexprMath.h"

/**
 * Calculates the potential energy between particle i and j using the Lennard Jones 12-6 potential.
 * @param i coordinate of the first particle
 * @param j coordinate of the second particle
 * @param cutoff the cutoff distance in wich we consider interactions
 * @param sigma sigma value for particles
 * @param epsilon epsilon value for particles
 * @return potential energy
 */
constexpr double calculateLJPotential(std::array<double, 3> i, std::array<double, 3> j, double cutoff, double sigma,
                                      double epsilon) {
  using namespace autopas::utils::ArrayMath::literals;

  // the vector from particle j to i
  auto imj = i - j;

  // the distance between both particles
  auto r = autopas::utils::ConstexprMath::sqrt(imj[0] * imj[0] + imj[1] * imj[1] + imj[2] * imj[2]);

  // if the distance between the two particles is larger then cutoff, we don't
  // consider this interaction
  if (r > cutoff) {
    return 0;
  }

  auto sdr = sigma / r;
  // (sigma / r)^6
  auto lj6 = sdr * sdr * sdr * sdr * sdr * sdr;
  // (sigma / r)^12
  auto lj12 = lj6 * lj6;
  // the lennard jones potential
  auto v = 4.0 * epsilon * (lj12 - lj6);

  return v;
}

/**
 * Calculates the force exerted by particle j on particle i using the 12-6 potential of Lennard Jones.
 * @param i coordinate of the first particle
 * @param j coordinate of the second particle
 * @param cutoff the cutoff distance in wich we consider interactions
 * @param sigma sigma value for particles
 * @param epsilon epsilon value for particles
 * @return The force exerted by particle j on particle i
 */
constexpr std::array<double, 3> calculateLJForce(std::array<double, 3> i, std::array<double, 3> j, double cutoff,
                                                 double sigma, double epsilon) {
  using namespace autopas::utils::ArrayMath::literals;

  // the vector from particle j to i
  auto imj = i - j;

  // the distance between both particles
  auto r = autopas::utils::ConstexprMath::sqrt(imj[0] * imj[0] + imj[1] * imj[1] + imj[2] * imj[2]);

  // if the distance between the two prticles is larger thn cutoff, we don't
  // consider this interaction
  if (r > cutoff) {
    return {0, 0, 0};
  }

  // r^6
  auto r6 = r * r * r * r * r * r;
  // sigma^6
  auto sigma6 = sigma * sigma * sigma * sigma * sigma * sigma;
  // sigma^6 / r^7
  auto dlj6 = sigma6 / (r6 * r);
  // sigma^12 / r^13
  auto dlj12 = (sigma6 * sigma6) / (r6 * r6 * r);
  // the derivative with respect to r of the lennard jones potential
  auto dUr = 48. * epsilon * (dlj12 - 0.5 * dlj6);

  // the forces in x, y and z direction
  auto fx = (imj[0] / r) * dUr;
  auto fy = (imj[1] / r) * dUr;
  auto fz = (imj[2] / r) * dUr;

  auto f = std::array<double, 3>{fx, fy, fz};

  return f;
}

/**
 * Calculates the virial between particle i and j using the Lennard Jones 12-6 potential.
 * @param i coordinate of the first particle
 * @param j coordinate of the second particle
 * @param cutoff the cutoff distance in wich we consider interactions
 * @param sigma sigma value for particles
 * @param epsilon epsilon value for particles
 * @return virial
 */
constexpr std::array<double, 3> calculateLJVirial(std::array<double, 3> i, std::array<double, 3> j, double cutoff,
                                                  double sigma, double epsilon) {
  using namespace autopas::utils::ArrayMath::literals;

  // first we need the forces
  auto f = calculateLJForce(i, j, cutoff, sigma, epsilon);
  // the vector from particle j to i
  auto imj = i - j;
  auto virial = std::array<double, 3>{imj[0] * f[0], imj[1] * f[1], imj[2] * f[2]};
  return virial;
}

/**
 * Calculates the sum of all components of the virial between particle i and j using the Lennard Jones 12-6 potential.
 * @param i coordinate of the first particle
 * @param j coordinate of the second particle
 * @param cutoff the cutoff distance in wich we consider interactions
 * @param sigma sigma value for particles
 * @param epsilon epsilon value for particles
 * @return sum of all three components of the virial vector
 */
constexpr double calculateLJVirialTotal(std::array<double, 3> i, std::array<double, 3> j, double cutoff, double sigma,
                                        double epsilon) {
  auto virial = calculateLJVirial(i, j, cutoff, sigma, epsilon);
  return virial[0] + virial[1] + virial[2];
}