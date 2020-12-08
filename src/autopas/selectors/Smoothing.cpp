/**
 * @file Smoothing.cpp
 * @author F. Gratl
 * @date 23/11/2020
 */

#include "Smoothing.h"

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
std::tuple<std::vector<double>, bool> calculateWeightsSimple(const std::vector<std::pair<size_t, size_t>> &points,
                                                             size_t pointsPerEstimation,
                                                             size_t maxDistFromIntervalEdge) {
  const size_t i = points.size() - 1;
  const size_t firstIndex = i - pointsPerEstimation + 1;
  // initialize all weights with 0
  std::vector<double> weights(pointsPerEstimation);

  const auto &xi = points[i].first;

  // define thresholds for shortcuts: if residuals are beyond these values the
  // are assumed to be 0 respectively 1
  auto maxDistFromIntervalEdgeHigh = maxDistFromIntervalEdge * .999;
  auto maxDistFromIntervalEdgeLow = maxDistFromIntervalEdge * .001;

  double sumOfWeights = 0.0;

  // compute all weighs that can be assumed to be non-zero
  for (size_t j = 0; j < weights.size(); ++j) {
    auto xj = points[j + firstIndex].first;
    // residual = |xi - xj|
    size_t residual = xj > xi ? xj - xi : xi - xj;

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
 * Each porjection is the sum of the respective weight and the residuals proportion of the weighted sum of squared
 * residuals.
 *
 * @param points
 * @param pointsPerEstimation
 * @param weights
 * @return
 */
double calculateYFitSimple(const std::vector<std::pair<size_t, size_t>> &points, size_t pointsPerEstimation,
                           const std::vector<double> &weights) {
  const size_t i = points.size() - 1;
  const size_t firstIndex = i - pointsPerEstimation + 1;
  std::vector<double> projections = weights;

  // weighted center of x
  double sumWeightedX = 0.;
  for (size_t j = 0; j < weights.size(); ++j) {
    sumWeightedX += weights[j] * points[j + firstIndex].first;
  }

  double weightedDistFromCenterXSquare = 0.;
  std::vector<size_t> deviations(weights.size());
  for (size_t j = 0; j < weights.size(); ++j) {
    deviations[j] = points[j + firstIndex].first - sumWeightedX;
    weightedDistFromCenterXSquare += weights[j] * deviations[j] * deviations[j];
  }

  // threshold whether points are not too clumped up
  size_t pointsRange = points.back().first - points.front().first;
  if (weightedDistFromCenterXSquare > 1e-6 * pointsRange * pointsRange) {
    // here i is always at the end of the interval
    auto distIToCenter = deviations.back();
    double distDivSqDev = distIToCenter / weightedDistFromCenterXSquare;

    for (size_t j = 0; j < weights.size(); ++j) {
      projections[j] = weights[j] * (1. + distDivSqDev * deviations[j]);
    }
  }

  double yFittedI = 0.;
  for (size_t j = 0; j < projections.size(); ++j) {
    yFittedI += projections[j] * points[j + firstIndex].second;
  }

  return yFittedI;
}

size_t autopas::smoothing::smoothLastPoint(const std::vector<std::pair<size_t, size_t>> &points,
                                           size_t pointsPerEstimation) {
  // if one or no points are used for smoothing do nothing
  if (pointsPerEstimation <= 2) {
    return points.back().second;
  }
  // if there are not enough points to smooth do nothing
  if (points.size() < 2) {
    if (points.size() == 1) {
      return points[0].second;
    }
    return 0;
  }

  // do not try to use more points than there are.
  pointsPerEstimation = std::min(pointsPerEstimation, points.size());

  // only fit last point
  size_t i = points.size() - 1;
  // find neighborhood
  const auto firstIndex = i - pointsPerEstimation + 1;

  // maxDistFromIntervalEdge = max(xi - xFirst, xLast - xi)
  auto maxDistFromIntervalEdge = std::max(points[i].first - points[firstIndex].first, 0ul);

  // Calculate weights
  auto [weights, fitOk] = calculateWeightsSimple(points, pointsPerEstimation, maxDistFromIntervalEdge);

  // either apply fit or take original datapoint
  if (fitOk) {
    return std::round(calculateYFitSimple(points, pointsPerEstimation, weights));
  } else {
    return points[i].second;
  }
}
