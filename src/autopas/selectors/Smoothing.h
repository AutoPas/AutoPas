/**
 * @file Smoothing.h
 * @author F. Gratl
 * @date 23/11/2020
 */

#pragma once

#include <cstdlib>
#include <tuple>
#include <vector>

namespace autopas::smoothing {

/**
 * Calculates the smoothed y value for the last point in the given points according to the LOESS algorithm.
 *
 * @note This function operates with unsigned values and does not protect against under-/overflow!
 *
 * @param points
 * @param pointsPerEstimation Number of points to take into account for smoothing.
 * @return
 */
long smoothLastPoint(const std::vector<std::pair<size_t, long>> &points, size_t pointsPerEstimation);

}  // namespace autopas::smoothing
