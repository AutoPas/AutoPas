/**
 * @file NumParticlesEstimator.h
 * @author F. Gratl
 * @date 06.03.23
 */

#pragma once

#include <array>

#include "autopas/utils/ArrayMath.h"

namespace autopas::utils::NumParticlesEstimator {

/**
 * Given a number of particles and the dimensions of a box, estimate the number of halo particles.
 * Assumptions:
 *   - Uniform particle distribution.
 *   - Cuboid particle container.
 *   - Halo width is the same in all dimensions.
 * @param numParticles
 * @param boxMin
 * @param boxMax
 * @param haloWidth
 * @return Estimated number of halo particles.
 */
size_t estimateNumHalosUniform(size_t numParticles, const std::array<double, 3> &boxMin,
                               const std::array<double, 3> &boxMax, double haloWidth);

}  // namespace autopas::utils::NumParticlesEstimator
