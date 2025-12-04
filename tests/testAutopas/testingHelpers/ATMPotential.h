/**
 * @file ATMPotential.h
 * @author muehlhaeusser
 * @date 29.08.23
 *
 * A simple reference implementation of the axilrod teller muto potential to calculate expected forces
 */

#pragma once
#include <tuple>

#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ConstexprMath.h"

/**
 * Calculates the potential energy between particles i, j and k using the Axilrod Teller Muto potential.
 * Here we use the original formula with trigonometric functions.
 * @param posI coordinate of the first particle
 * @param posJ coordinate of the second particle
 * @param posK coordinate of the third particle
 * @param cutoff the cutoff distance in which we consider interactions
 * @param nu nu value for particles
 * @return potential energy
 */
constexpr double calculateATMPotential(const std::array<double, 3> &posI, const std::array<double, 3> &posJ,
                                       const std::array<double, 3> &posK, double cutoff, double nu) {
  using namespace autopas::utils::ArrayMath::literals;

  // the distance vectors
  const auto posIToPosJ = posI - posJ;
  const auto posIToPosK = posI - posK;
  const auto posJToPosK = posJ - posK;

  // distances between particles
  const auto distIJ =
      autopas::utils::ConstexprMath::sqrt(autopas::utils::ArrayMath::dot(posIToPosJ, posIToPosJ), 1e-16);
  const auto distIK =
      autopas::utils::ConstexprMath::sqrt(autopas::utils::ArrayMath::dot(posIToPosK, posIToPosK), 1e-16);
  const auto distJK =
      autopas::utils::ConstexprMath::sqrt(autopas::utils::ArrayMath::dot(posJToPosK, posJToPosK), 1e-16);

  // if the distance between the two particles is larger than cutoff, we don't
  // consider this interaction
  if (distIJ > cutoff or distIK > cutoff or distJK > cutoff) {
    return 0;
  }

  // calculate the cosines between the particles
  // cos_i = (r_ki * r_ji) / (|r_ki| |r_ji|)

  const auto cosI = autopas::utils::ArrayMath::dot(posIToPosK, posIToPosJ) / (distIK * distIJ);
  const auto cosJ = -autopas::utils::ArrayMath::dot(posJToPosK, posIToPosJ) / (distJK * distIJ);
  const auto cosK = autopas::utils::ArrayMath::dot(posIToPosK, posJToPosK) / (distIK * distJK);

  const auto rrr = distIJ * distJK * distIK;
  const auto rrr3 = rrr * rrr * rrr;

  // the axilrod teller potential
  const auto potentialEnergy = nu * (3 * cosI * cosJ * cosK + 1.) / rrr3;

  return potentialEnergy;
}

/**
 * Calculates the forces exerted on three particles using the Axilrod-Teller-Muto potential.
 * @param posI coordinate of the first particle
 * @param posJ coordinate of the second particle
 * @param posK coordinate of the third particle
 * @param cutoff the cutoff distance in which we consider interactions
 * @param nu Axilrod-Teller-Muto Factor
 * @return The forces exerted on particle i, particle j, particle k
 */
constexpr std::array<std::array<double, 3>, 3> calculateATMForce(const std::array<double, 3> &posI,
                                                                 const std::array<double, 3> &posJ,
                                                                 const std::array<double, 3> &posK, double cutoff,
                                                                 double nu) {
  using namespace autopas::utils::ArrayMath::literals;

  // the distance vectors
  const auto posJToPosI = posJ - posI;
  const auto posIToPosK = posI - posK;
  const auto posKToPosJ = posK - posJ;

  // distances between particles
  const auto distJISquared = autopas::utils::ArrayMath::dot(posJToPosI, posJToPosI);
  const auto distJI = autopas::utils::ConstexprMath::sqrt(distJISquared, 1e-16);
  const auto distIKSquared = autopas::utils::ArrayMath::dot(posIToPosK, posIToPosK);
  const auto distIK = autopas::utils::ConstexprMath::sqrt(distIKSquared, 1e-16);
  const auto distKJSquared = autopas::utils::ArrayMath::dot(posKToPosJ, posKToPosJ);
  const auto distKJ = autopas::utils::ConstexprMath::sqrt(distKJSquared, 1e-16);

  // if the distance between the two particles is larger than cutoff, we don't
  // consider this interaction
  if (distJI > cutoff or distIK > cutoff or distKJ > cutoff) {
    return {std::array<double, 3>{0., 0., 0.}, std::array<double, 3>{0., 0., 0.}, std::array<double, 3>{0., 0., 0.}};
  }

  // Dot products that are the numerators of the cosine formula
  const double cosI = autopas::utils::ArrayMath::dot(posJToPosI, posIToPosK);
  const double cosJ = autopas::utils::ArrayMath::dot(posKToPosJ, posJToPosI);
  const double cosK = autopas::utils::ArrayMath::dot(posIToPosK, posKToPosJ);

  const double cos6 = cosI * cosJ * cosK;
  const double distSquaredAll = distJI * distJI * distIK * distIK * distKJ * distKJ;
  const double dist5All = distSquaredAll * distSquaredAll * std::sqrt(distSquaredAll);

  auto forceI = std::array<double, 3>{0., 0., 0.};
  auto forceJ = std::array<double, 3>{0., 0., 0.};
  auto forceK = std::array<double, 3>{0., 0., 0.};

  // loop over all dimensions
  for (size_t i = 0; i < 3; i++) {
    forceI[i] = posKToPosJ[i] * cosI * (cosJ - cosK) +
                posJToPosI[i] * (cosJ * cosK - distKJSquared * distIKSquared + 5.0 * cos6 / distJISquared) +
                posIToPosK[i] * (-cosJ * cosK + distJISquared * distKJSquared - 5.0 * cos6 / distIKSquared);

    forceJ[i] = posIToPosK[i] * cosJ * (cosK - cosI) +
                posKToPosJ[i] * (cosI * cosK - distJISquared * distIKSquared + 5.0 * cos6 / distKJSquared) +
                posJToPosI[i] * (-cosI * cosK + distKJSquared * distIKSquared - 5.0 * cos6 / distJISquared);

    forceK[i] = posJToPosI[i] * cosK * (cosI - cosJ) +
                posIToPosK[i] * (cosI * cosJ - distJISquared * distKJSquared + 5.0 * cos6 / distIKSquared) +
                posKToPosJ[i] * (-cosI * cosJ + distJISquared * distIKSquared - 5.0 * cos6 / distKJSquared);
  }
  forceI *= 3.0 * nu / dist5All;
  forceJ *= 3.0 * nu / dist5All;
  forceK *= 3.0 * nu / dist5All;

  return {forceI, forceJ, forceK};
}

/**
 * Calculates the virial between three particles i,j,k from the Axilrod-Teller-Muto potential.
 * @param posI coordinate of the first particle
 * @param posJ coordinate of the second particle
 * @param posK coordinate of the third particle
 * @param cutoff the cutoff distance in which we consider interactions
 * @param nu Axilrod-Teller-Muto Factor
 * @return virial
 */
constexpr std::array<std::array<double, 3>, 3> calculateATMVirials(const std::array<double, 3> &posI,
                                                                   const std::array<double, 3> &posJ,
                                                                   const std::array<double, 3> &posK, double cutoff,
                                                                   double nu) {
  using namespace autopas::utils::ArrayMath::literals;

  // first we need the forces
  const auto forces = calculateATMForce(posI, posJ, posK, cutoff, nu);
  const auto forceI = forces[0];
  const auto forceJ = forces[1];
  const auto forceK = forces[2];

  const auto virialI = forceI * (posI * 2. - posJ - posK) / 3.;
  const auto virialJ = forceJ * (posJ * 2. - posI - posK) / 3.;
  const auto virialK = forceK * (posK * 2. - posI - posJ) / 3.;

  return {virialI, virialJ, virialK};
}

/**
 * Calculates the sum of all components of the virial between particle i, j, k using the Axilrod-Teller-Muto potential.
 * @param posI coordinate of the first particle
 * @param posJ coordinate of the second particle
 * @param posK coordinate of the third particle
 * @param cutoff the cutoff distance in which we consider interactions
 * @param nu Axilrod-Teller-Muto Factor
 * @return sum of all three components of the virial vector
 */
constexpr double calculateATMVirialTotal(const std::array<double, 3> &posI, const std::array<double, 3> &posJ,
                                         const std::array<double, 3> &posK, double cutoff, double nu) {
  using namespace autopas::utils::ArrayMath::literals;
  const auto [virialI, virialJ, virialK] = calculateATMVirials(posI, posJ, posK, cutoff, nu);
  const auto virialSum = virialI + virialJ + virialK;
  return virialSum[0] + virialSum[1] + virialSum[2];
}

/**
 * Returns the sum of all components of the virial for each particle i, j, k individually using the Axilrod-Teller-Muto
 * potential.
 * @param posI coordinate of the first particle
 * @param posJ coordinate of the second particle
 * @param posK coordinate of the third particle
 * @param cutoff the cutoff distance in which we consider interactions
 * @param nu Axilrod-Teller-Muto Factor
 * @return sum of all three components of the virial vector
 */
constexpr std::array<double, 3> calculateATMVirialTotalPerParticle(const std::array<double, 3> &posI,
                                                                   const std::array<double, 3> &posJ,
                                                                   const std::array<double, 3> &posK, double cutoff,
                                                                   double nu) {
  using namespace autopas::utils::ArrayMath::literals;
  const auto [virialI, virialJ, virialK] = calculateATMVirials(posI, posJ, posK, cutoff, nu);
  const auto virialSumI = virialI[0] + virialI[1] + virialI[2];
  const auto virialSumJ = virialJ[0] + virialJ[1] + virialJ[2];
  const auto virialSumK = virialK[0] + virialK[1] + virialK[2];
  return {virialSumI, virialSumJ, virialSumK};
}