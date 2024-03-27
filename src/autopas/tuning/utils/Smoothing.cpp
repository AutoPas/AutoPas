/**
 * @file Smoothing.cpp
 * @author F. Gratl
 * @date 23/11/2020
 */

#include "Smoothing.h"

#include "autopas/tuning/searchSpace/Evidence.h"
#include "autopas/utils/Math.h"

/**
 * Calculates the weights for the k-neighbors of the last point in the vector using the tri-cube function.
 * If the sum of weights is <= 0 all neighbors are at the same x coordinate as the point and smoothing will not change
 * the point. This case is indicated by the bool in the return tuple.
 *
 * @param points Sorted collection of observations.
 * @param pointsPerEstimation Number of neighbors to consider for smoothing.
 * @param maxDistFromIntervalEdge Largest distance between the point that shall be fitted and its neighbors.
 * @return Tuple of a vector containing the weights for the neighbors
 * and a bool indicating whether fitting is unnecessary.
 */
std::tuple<std::vector<double>, bool> calculateWeightsSimple(const std::vector<autopas::Evidence> &points,
                                                             size_t pointsPerEstimation,
                                                             size_t maxDistFromIntervalEdge) {
  // since we will only smooth the last point there is no outer loop and indexToFit shall be fixed
  const size_t indexToFit = points.size() - 1;
  const size_t firstIndex = indexToFit - pointsPerEstimation + 1;
  // initialize all weights with 0
  std::vector<double> weights(pointsPerEstimation);

  const auto &xi = points[indexToFit].iteration;

  // Define thresholds for shortcuts: If residuals are beyond these values, they
  // are assumed to be 0, respectively 1.
  auto maxDistFromIntervalEdgeHigh = maxDistFromIntervalEdge * .999;
  auto maxDistFromIntervalEdgeLow = maxDistFromIntervalEdge * .001;

  double sumOfWeights = 0.0;

  // compute all weighs that can be assumed to be non-zero
  for (size_t j = 0; j < weights.size(); ++j) {
    auto xj = points[j + firstIndex].iteration;
    // residual = |xi - xj|
    const size_t residual = xj > xi ? xj - xi : xi - xj;

    // if x are too far apart ...
    if (residual <= maxDistFromIntervalEdgeHigh) {
      if (residual > maxDistFromIntervalEdgeLow) {
        // tri-cube function
        weights[j] = autopas::utils::Math::pow<3>(
            1.0 - autopas::utils::Math::pow<3>(static_cast<double>(residual) / maxDistFromIntervalEdge));
      } else {
        weights[j] = 1.0;
      }

      sumOfWeights += weights[j];
    } else if (xj > xi) {
      // break because from here on weights will only be more 0
      break;
    }
  }

  bool fitIsOk = true;
  // if all residuals were 0 fitting is not necessary
  if (sumOfWeights <= .0) {
    fitIsOk = false;
  } else {
    // normalize weights
    for (auto &weight : weights) {
      weight /= sumOfWeights;
    }
  }
  return std::make_tuple(weights, fitIsOk);
}

/**
 * Calculates the smoothed y-value of the last point in the vector.
 * The fitted value is the sum of projections of the y-values of the neighbors in the chosen interval.
 * Each projection is the sum of the respective weight and the proportion of the residuals of the weighted sum of
 * squared residuals.
 *
 * @param points
 * @param pointsPerEstimation
 * @param weights
 * @return
 */
double calculateYFitSimple(const std::vector<autopas::Evidence> &points, size_t pointsPerEstimation,
                           const std::vector<double> &weights) {
  // since we will only smooth the last point there is no outer loop and indexToFit shall be fixed
  const size_t indexToFit = points.size() - 1;
  const size_t firstIndex = indexToFit - pointsPerEstimation + 1;
  std::vector<double> projections = weights;

  // weighted center of x
  double sumWeightedX = 0.;
  for (size_t j = 0; j < weights.size(); ++j) {
    sumWeightedX += weights[j] * points[j + firstIndex].iteration;
  }

  double weightedDistFromCenterXSquare = 0.;
  std::vector<size_t> deviations(weights.size());
  for (size_t j = 0; j < weights.size(); ++j) {
    deviations[j] = points[j + firstIndex].iteration - sumWeightedX;
    weightedDistFromCenterXSquare += weights[j] * deviations[j] * deviations[j];
  }

  // threshold whether points are not too clumped up
  size_t pointsRange = points.back().iteration - points.front().iteration;
  if (weightedDistFromCenterXSquare > 1e-6 * pointsRange * pointsRange) {
    // here indexToFit is always at the end of the interval
    auto distIndexToCenter = deviations.back();
    double distDivSqDev = distIndexToCenter / weightedDistFromCenterXSquare;

    for (size_t j = 0; j < weights.size(); ++j) {
      projections[j] = weights[j] * (1. + distDivSqDev * deviations[j]);
    }
  }

  double yFitted = 0.;
  for (size_t j = 0; j < projections.size(); ++j) {
    yFitted += projections[j] * points[j + firstIndex].value;
  }

  return yFitted;
}

long autopas::smoothing::smoothLastPoint(const std::vector<Evidence> &points, size_t pointsPerEstimation) {
  // if one or no points are used for smoothing do nothing
  if (pointsPerEstimation <= 2) {
    return points.back().value;
  }
  // if there are not enough points to smooth do nothing
  if (points.size() < 2) {
    if (points.size() == 1) {
      return points[0].value;
    }
    return 0;
  }

  // do not try to use more points than there are.
  pointsPerEstimation = std::min(pointsPerEstimation, points.size());

  // since we will only smooth the last point there is no outer loop and indexToFit shall be fixed
  size_t indexToFit = points.size() - 1;
  // find neighborhood
  const auto firstIndex = indexToFit - pointsPerEstimation + 1;

  // maxDistFromIntervalEdge = max(xi - xFirst, xLast - xi)
  auto maxDistFromIntervalEdge = std::max(points[indexToFit].iteration - points[firstIndex].iteration, 0ul);

  // Calculate weights
  auto [weights, fitOk] = calculateWeightsSimple(points, pointsPerEstimation, maxDistFromIntervalEdge);

  // either apply fit or take original datapoint
  if (fitOk) {
    return std::round(calculateYFitSimple(points, pointsPerEstimation, weights));
  } else {
    return points[indexToFit].value;
  }
}
