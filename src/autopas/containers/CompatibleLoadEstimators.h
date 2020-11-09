/**
 * @file CompatibleLoadEstimators.h
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

/**
 * Collection of functions for estimating the load required to update a specific region within a containers.
 */
namespace autopas::loadEstimators {

/**
 * Returns set of load estimators compatible with the container.
 *
 * @param container
 * @return compatible load estimators
 */
static std::set<autopas::LoadEstimatorOption> allCompatibleLoadEstimators(autopas::ContainerOption container) {
  switch (container) {
    case ContainerOption::linkedCells: {
      return std::set<autopas::LoadEstimatorOption>{LoadEstimatorOption::none,
                                                    LoadEstimatorOption::squaredParticlesPerCell};
    }
    case ContainerOption::verletListsCells: {
      return std::set<autopas::LoadEstimatorOption>{LoadEstimatorOption::none,
                                                    LoadEstimatorOption::squaredParticlesPerCell,
                                                    LoadEstimatorOption::neighborListLength};
    }
    case ContainerOption::verletClusterLists: {
      return std::set<autopas::LoadEstimatorOption>{LoadEstimatorOption::none, LoadEstimatorOption::neighborListLength};
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
    case TraversalOption::lc_sliced_balanced:
      [[fallthrough]];
    case TraversalOption::vcl_sliced_balanced:
      [[fallthrough]];
    case TraversalOption::vlc_sliced_balanced: {
      return true;
    }
    default: {
      return false;
    }
  }
}

/**
 * If traversal uses load estimation, returns all load estimators in allowedOptions,
 * that are compatible with the container, but always allows none if the intersection is empty.
 *
 * @param container
 * @param traversal
 * @param allowedOptions
 * @return applicable traversals or {none}
 */
static std::set<autopas::LoadEstimatorOption> getApplicableLoadEstimators(
    autopas::ContainerOption container, autopas::TraversalOption traversal,
    const std::set<autopas::LoadEstimatorOption> allowedOptions) {
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
