/**
 * @file NumParticlesEstimator.cpp
 * @author F. Gratl
 * @date 06.03.23
 */

#include "NumParticlesEstimator.h"

size_t autopas::utils::NumParticlesEstimator::estimateNumHalosUniform(size_t numParticles,
                                                                      const std::array<double, 3> &boxMin,
                                                                      const std::array<double, 3> &boxMax,
                                                                      double haloWidth) {
  using autopas::utils::ArrayMath::addScalar;
  using autopas::utils::ArrayMath::sub;

  // calculate global densities as an estimator for how many particles might end up in the halo
  const auto boxLength = sub(boxMax, boxMin);
  const auto boxVolume = boxLength[0] * boxLength[1] * boxLength[2];
  const auto haloBoxLength = addScalar(boxLength, 2. * haloWidth);
  const auto haloVolume = haloBoxLength[0] * haloBoxLength[1] * haloBoxLength[2];
  // assumption: particle density is the same in the box and the halo. Mathematically speaking:
  // haloVolume / boxVolume == numHaloParticles / numParticles
  return static_cast<size_t>(static_cast<double>(numParticles) * (haloVolume / boxVolume));
}
