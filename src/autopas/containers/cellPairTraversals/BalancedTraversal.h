/**
 * @file BalancedTraversal.h
 *
 * @date 09 May 2020
 * @author fischerv
 */

#pragma once

#include <array>
#include <functional>
#include <utility>

namespace autopas {

/**
 * Base class for traversals utilising load balancing
 */
class BalancedTraversal {
 public:
  /**
   * Type signature for load estimators
   */
  using EstimatorFunction =
      std::function<unsigned long(const std::array<unsigned long, 3> &, const std::array<unsigned long, 3> &,
                                  const std::array<unsigned long, 3> &)>;

  /**
   * Setter for load estimation algorithm.
   *
   * @param loadEstimator
   */
  void setLoadEstimator(EstimatorFunction loadEstimator) { _loadEstimator = loadEstimator; }

 protected:
  /**
   * Algorithm to use for estimating load
   *
   * parameters: cellsPerDimension, lowerCorner, upperCorner
   */
  EstimatorFunction _loadEstimator;
};
}  // namespace autopas
