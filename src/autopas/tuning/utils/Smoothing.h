/**
 * @file Smoothing.h
 * @author F. Gratl
 * @date 23/11/2020
 */

#pragma once

#include <cstddef>
#include <tuple>
#include <vector>

#include "autopas/tuning/searchSpace/Evidence.h"

namespace autopas::smoothing {

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
                                                             size_t maxDistFromIntervalEdge);

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
                           const std::vector<double> &weights);

/**
 * Calculates the smoothed y value for the last point in the given points according to the LOESS algorithm.
 *
 * @note This function operates with unsigned values and does not protect against under-/overflow!
 *
 * @param points
 * @param pointsPerEstimation Number of points to take into account for smoothing.
 * @return
 */
long smoothLastPoint(const std::vector<Evidence> &points, size_t pointsPerEstimation);

}  // namespace autopas::smoothing
