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
 * @param posI coordinate of the first particle
 * @param posJ coordinate of the second particle
 * @param cutoff the cutoff distance in wich we consider interactions
 * @param sigma sigma value for particles
 * @param epsilon epsilon value for particles
 * @return potential energy
 */
constexpr double calculateLJPotential(const std::array<double, 3> &posI, const std::array<double, 3> &posJ,
                                      double cutoff, double sigma, double epsilon) {
  using namespace autopas::utils::ArrayMath::literals;

  // the vector from particle j to i
  const auto posIMinusPosJ = posI - posJ;

  // the distance between both particles
  const auto dist =
      autopas::utils::ConstexprMath::sqrt(autopas::utils::ArrayMath::dot(posIMinusPosJ, posIMinusPosJ), 1e-16);

  // if the distance between the two particles is larger then cutoff, we don't
  // consider this interaction
  if (dist > cutoff) {
    return 0;
  }

  // (sigma / dist)
  const auto sdr = sigma / dist;
  // (sigma / dist)^6
  const auto lj6 = sdr * sdr * sdr * sdr * sdr * sdr;
  // (sigma / dist)^12
  const auto lj12 = lj6 * lj6;
  // the lennard jones potential
  const auto potentialEnergy = 4.0 * epsilon * (lj12 - lj6);

  return potentialEnergy;
}

/**
 * Calculates the force exerted by particle j on particle i using the 12-6 potential of Lennard Jones.
 * @param posI coordinate of the first particle
 * @param posJ coordinate of the second particle
 * @param cutoff the cutoff distance in wich we consider interactions
 * @param sigma sigma value for particles
 * @param epsilon epsilon value for particles
 * @return The force exerted by particle j on particle i
 */
constexpr std::array<double, 3> calculateLJForce(const std::array<double, 3> &posI, const std::array<double, 3> &posJ,
                                                 double cutoff, double sigma, double epsilon) {
  using namespace autopas::utils::ArrayMath::literals;

  // the vector from particle j to i
  const auto posIMinusPosJ = posI - posJ;

  // the distance between both particles
  const auto dist =
      autopas::utils::ConstexprMath::sqrt(autopas::utils::ArrayMath::dot(posIMinusPosJ, posIMinusPosJ), 1e-16);

  // if the distance between the two prticles is larger thn cutoff, we don't
  // consider this interaction
  if (dist > cutoff) {
    return {0, 0, 0};
  }

  // r^6
  const auto r6 = dist * dist * dist * dist * dist * dist;
  // sigma^6
  const auto sigma6 = sigma * sigma * sigma * sigma * sigma * sigma;
  // sigma^6 / r^7
  const auto dlj6 = sigma6 / (r6 * dist);
  // sigma^12 / r^13
  const auto dlj12 = (sigma6 * sigma6) / (r6 * r6 * dist);
  // the derivative with respect to r of the lennard jones potential
  const auto dUr = 48. * epsilon * (dlj12 - 0.5 * dlj6);

  // the forces in x, y and z direction
  const auto fx = (posIMinusPosJ[0] / dist) * dUr;
  const auto fy = (posIMinusPosJ[1] / dist) * dUr;
  const auto fz = (posIMinusPosJ[2] / dist) * dUr;

  const auto force = std::array<double, 3>{fx, fy, fz};

  return force;
}

/**
 * Calculates the virial between particle i and j using the Lennard Jones 12-6 potential.
 * @param posI coordinate of the first particle
 * @param posJ coordinate of the second particle
 * @param cutoff the cutoff distance in wich we consider interactions
 * @param sigma sigma value for particles
 * @param epsilon epsilon value for particles
 * @return virial
 */
constexpr std::array<double, 3> calculateLJVirial(const std::array<double, 3> &posI, const std::array<double, 3> &posJ,
                                                  double cutoff, double sigma, double epsilon) {
  using namespace autopas::utils::ArrayMath::literals;

  // first we need the forces
  const auto force = calculateLJForce(posI, posJ, cutoff, sigma, epsilon);
  // the vector from particle j to i
  const auto posIMinusPosJ = posI - posJ;
  const auto virial =
      std::array<double, 3>{posIMinusPosJ[0] * force[0], posIMinusPosJ[1] * force[1], posIMinusPosJ[2] * force[2]};
  return virial;
}

/**
 * Calculates the sum of all components of the virial between particle i and j using the Lennard Jones 12-6 potential.
 * @param posI coordinate of the first particle
 * @param posJ coordinate of the second particle
 * @param cutoff the cutoff distance in wich we consider interactions
 * @param sigma sigma value for particles
 * @param epsilon epsilon value for particles
 * @return sum of all three components of the virial vector
 */
constexpr double calculateLJVirialTotal(const std::array<double, 3> &posI, const std::array<double, 3> &posJ,
                                        double cutoff, double sigma, double epsilon) {
  const auto virial = calculateLJVirial(posI, posJ, cutoff, sigma, epsilon);
  return virial[0] + virial[1] + virial[2];
}