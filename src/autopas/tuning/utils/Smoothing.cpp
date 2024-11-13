/**
 * @file Smoothing.cpp
 * @author F. Gratl
 * @date 23/11/2020
 */

#include "Smoothing.h"

#include "autopas/tuning/searchSpace/Evidence.h"
#include "autopas/utils/Math.h"

std::tuple<std::vector<double>, bool> autopas::smoothing::calculateWeightsSimple(
    const std::vector<autopas::Evidence> &points, size_t pointsPerEstimation, size_t maxDistFromIntervalEdge) {
  // since we will only smooth the last point there is no outer loop of LOESS and indexToFit shall be fixed
  const auto indexToFit = points.size() - 1;
  const auto firstIndex = indexToFit - pointsPerEstimation + 1;
  // initialize all weights with 0
  std::vector<double> weights(pointsPerEstimation);

  // get iteration number of the point being smoothed
  const auto &xi = static_cast<long>(points[indexToFit].iteration);

  double sumOfWeights = 0.0;

  // compute all weights that can be assumed to be non-zero
  for (size_t j = 0; j < pointsPerEstimation; ++j) {
    // set weight to zero
    weights[j] = 0.;
    const auto xj = static_cast<long>(points[j + firstIndex].iteration);
    // residual = |xi - xj|
    const size_t residual = std::abs(xi - xj);

    // check xi and xj are within maxDistFromIntervalEdge
    if (residual <= maxDistFromIntervalEdge) {
      // tri-cube function
      weights[j] = autopas::utils::Math::pow<3>(
          1.0 -
          autopas::utils::Math::pow<3>(static_cast<double>(residual) / static_cast<double>(maxDistFromIntervalEdge)));
    }
    sumOfWeights += weights[j];
  }

  bool fitIsOk = true;
  // if all residuals were 0 fitting is not necessary because all data comes from the same iteration
  if (utils::Math::isNearAbs(sumOfWeights, .0, 1e-12)) {
    fitIsOk = false;
  } else {
    // normalize weights
    for (auto &weight : weights) {
      weight /= sumOfWeights;
    }
  }
  return std::make_tuple(weights, fitIsOk);
}

double autopas::smoothing::calculateYFitSimple(const std::vector<autopas::Evidence> &points, size_t pointsPerEstimation,
                                               const std::vector<double> &weights) {
  // since we will only smooth the last point there is no outer loop and indexToFit shall be fixed
  const size_t indexToFit = points.size() - 1;
  const size_t firstIndex = indexToFit - pointsPerEstimation + 1;
  std::vector<double> projections = weights;

  // weighted center of x
  double sumWeightedX = 0.;
  for (size_t j = 0; j < weights.size(); ++j) {
    sumWeightedX += weights[j] * static_cast<double>(points[j + firstIndex].iteration);
  }

  double weightedDistFromCenterXSquare = 0.;
  std::vector<double> deviations(weights.size());
  for (size_t j = 0; j < weights.size(); ++j) {
    deviations[j] = static_cast<double>(points[j + firstIndex].iteration) - sumWeightedX;
    weightedDistFromCenterXSquare += weights[j] * deviations[j] * deviations[j];
  }

  // threshold whether points are not too clumped up
  const auto pointsRange = static_cast<double>(points.back().iteration - points.front().iteration);
  if (weightedDistFromCenterXSquare > 1e-6 * pointsRange * pointsRange) {
    // here indexToFit is always at the end of the interval
    const auto distIndexToCenter = deviations.back();
    const auto distDivSqDev = distIndexToCenter / weightedDistFromCenterXSquare;

    for (size_t j = 0; j < weights.size(); ++j) {
      projections[j] = weights[j] * (1. + distDivSqDev * deviations[j]);
    }
  }

  double yFitted = 0.;
  for (size_t j = 0; j < projections.size(); ++j) {
    yFitted += projections[j] * static_cast<double>(points[j + firstIndex].value);
  }

  return yFitted;
}

long autopas::smoothing::smoothLastPoint(const std::vector<Evidence> &points, size_t pointsPerEstimation) {
  // smoothening requires at least two points
  if (pointsPerEstimation < 2 or points.size() == 1) {
    return points.back().value;
  }
  // special case: there are no points so return 0 instead of points.back().value
  if (points.empty()) {
    return 0;
  }

  // do not try to use more points than there are.
  pointsPerEstimation = std::min(pointsPerEstimation, points.size());

  // get the index of the point to be smoothed (the last point)
  const auto indexToFit = points.size() - 1;
  // get the index of the first point in the neighborhood around indexToFit used for smoothening
  const auto firstIndex = indexToFit - pointsPerEstimation + 1;

  // maxDistFromIntervalEdge = xi - xFirst. By distance, we refer to the distance in iteration number.
  const auto maxDistFromIntervalEdge = points[indexToFit].iteration - points[firstIndex].iteration;

  // Calculate weights
  const auto [weights, fitOk] = calculateWeightsSimple(points, pointsPerEstimation, maxDistFromIntervalEdge);

  // either apply fit or take original datapoint
  if (fitOk) {
    return std::lround(calculateYFitSimple(points, pointsPerEstimation, weights));
  } else {
    return points[indexToFit].value;
  }
}
