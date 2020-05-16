/**
 * @file loadEstimators.h
 *
 * @date 24 Apr 2020
 * @author fischerv
 */

#pragma once

#include <array>
#include <set>
#include <vector>

#include "autopas/options/ContainerOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

/**
 * Collection of functions for estimating the load required to update a specific region within a containers
 */
namespace autopas::loadEstimators {

/**
 * Sums up the squared number of particles for all cells within region
 *
 * @param cells
 * @param cellsPerDimension
 * @param lowerCorner lower boundary indices for region
 * @param upperCorner upper boundary indices for region
 * @return estimated load for given region
 */
template <class ParticleCell>
unsigned long squaredParticlesPerCell(const std::vector<ParticleCell> &cells,
                                      const std::array<unsigned long, 3> &cellsPerDimension,
                                      const std::array<unsigned long, 3> &lowerCorner,
                                      const std::array<unsigned long, 3> &upperCorner) {
  unsigned long sum = 0;
  for (unsigned long x = lowerCorner[0]; x <= upperCorner[0]; x++) {
    for (unsigned long y = lowerCorner[1]; y <= upperCorner[1]; y++) {
      for (unsigned long z = lowerCorner[2]; z <= upperCorner[2]; z++) {
        auto load =
            cells[autopas::utils::ThreeDimensionalMapping::threeToOneD(x, y, z, cellsPerDimension)].numParticles();
        sum += load * load;
      }
    }
  }
  return sum;
}

/**
 * returns set of load estimators compatible with container.
 *
 * @param container
 * @return compatible load estimators
 */
static const std::set<autopas::LoadEstimatorOption> allCompatibleLoadEstimators(autopas::ContainerOption container) {
  switch (container) {
    case ContainerOption::linkedCells: {
      return std::set<autopas::LoadEstimatorOption>{LoadEstimatorOption::none,
                                                    LoadEstimatorOption::squaredParticlesPerCell};
    }
    case ContainerOption::verletListsCells: {
      return std::set<autopas::LoadEstimatorOption>{LoadEstimatorOption::none,
                                                    LoadEstimatorOption::squaredParticlesPerCell};
    }
    default: {
      return std::set<autopas::LoadEstimatorOption>{};
    }
  }
}

/**
 * returns whether or not the given traversal uses load estimation.
 *
 * @param traversal
 * @return
 */
static bool usesLoadEstimator(autopas::TraversalOption traversal) {
  switch (traversal) {
    case TraversalOption::BalancedSliced: {
      return true;
    }
    case TraversalOption::BalancedSlicedVerlet: {
      return true;
    }
    default: {
      return false;
    }
  }
}

/**
 * If traversal uses load estimation, returns all load estimators in allowedOptions,
 * that are compatible with container, but always allows none if intersection is empty.
 *
 * @param container
 * @param traversal
 * @param allowedOptions
 * @return applicable traversals or {none}
 */
static const std::set<autopas::LoadEstimatorOption> getApplicableLoadEstimators(
    autopas::ContainerOption container, autopas::TraversalOption traversal,
    std::set<autopas::LoadEstimatorOption> allowedOptions) {
  if (usesLoadEstimator(traversal)) {
    auto compatible = allCompatibleLoadEstimators(container);
    std::set<autopas::LoadEstimatorOption> intersection;
    std::set_intersection(allowedOptions.begin(), allowedOptions.end(), compatible.begin(), compatible.end(),
                          std::inserter(intersection, intersection.begin()));
    if (not intersection.empty()) {
      return intersection;
    }
  }
  return std::set<autopas::LoadEstimatorOption>{LoadEstimatorOption::none};
}

}  // namespace autopas::loadEstimators
