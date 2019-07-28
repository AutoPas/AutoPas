/**
 * @file VerletClustersTraversalInterface.h
 * @author humig
 * @date 20.06.19
 */

#pragma once

/**
 * Some math operations for verlet clusters lists.
 */
namespace autopas::VerletClusterMaths {
/**
 * Converts 2d tower position to the index in the vector.
 * @param x x-coordinate of the tower
 * @param y y-coordinate of the tower
 * @param towersPerDim the towers per dimension of the container
 * @return index in vector
 */
static inline size_t index1D(const size_t x, const size_t y, std::array<size_t, 3> towersPerDim) {
  return x + y * towersPerDim[0];
}

}  // namespace autopas::VerletClusterMaths
