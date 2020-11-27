/**
 * @file Smoothing.cpp
 * @author F. Gratl
 * @date 23/11/2020
 */

#include "Smoothing.h"

#include "autopas/utils/Math.h"

/**
 * TODO
 * @param points
 * @param pointsPerEstimation
 * @param maxDistFromIntervalEdge
 * @return
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
 * TODO
 * @param points
 * @param pointsPerEstimation
 * @param maxDistFromIntervalEdge
 * @param weights
 * @return
 */
double calculateYFitSimple(const std::vector<std::pair<size_t, size_t>> &points, size_t pointsPerEstimation,
                           const size_t maxDistFromIntervalEdge, const std::vector<double> &weights) {
  const size_t i = points.size() - 1;
  const size_t firstIndex = i - pointsPerEstimation + 1;
  std::vector<double> projections = weights;

  if (maxDistFromIntervalEdge > 0) {
    // weighted center of x
    // TODO: make this a size_t ?
    double sumWeightedX = 0.;
    for (size_t j = 0; j < weights.size(); ++j) {
      sumWeightedX += weights[j] * points[j + firstIndex].first;
    }

    double weightedDistFromCenterXSquare = 0.;
    for (size_t j = 0; j < weights.size(); ++j) {
      auto deviation = points[j + firstIndex].first - sumWeightedX;
      weightedDistFromCenterXSquare += weights[j] * deviation * deviation;
    }

    // threshold whether points are not too clumped up
    size_t pointsRange = points.back().first - points.front().first;
    if (weightedDistFromCenterXSquare > 1e-6 * pointsRange * pointsRange) {
      double distIToCenter = points[i].first - sumWeightedX;
      double distDivSqDev = distIToCenter / weightedDistFromCenterXSquare;

      for (size_t j = 0; j < weights.size(); ++j) {
        double deviation = points[j + firstIndex].first - sumWeightedX;
        projections[j] = weights[j] * (1. + distDivSqDev * deviation);
      }
    }
  }

  double yFittedI = 0.;
  for (size_t j = 0; j < projections.size(); ++j) {
    yFittedI += projections[j] * points[j + firstIndex].second;
  }

  return yFittedI;
}

size_t autopas::smoothing::smoothLastPoint(const std::vector<std::pair<size_t, size_t>> &points, double span) {
  if (span <= 0 or span > 1.0) {
    throw std::runtime_error("span should be 0 < span <= 1");
  }

  if (points.size() < 2) {
    if (points.size() == 1) {
      return points[0].second;
    }
    return 0;
  }

  // pick at least two points and not more than total number of points
  size_t pointsPerEstimation = std::max(static_cast<size_t>(span * points.size()), 2ul);

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
    return std::round(calculateYFitSimple(points, pointsPerEstimation, maxDistFromIntervalEdge, weights));
  } else {
    return points[i].second;
  }
}
