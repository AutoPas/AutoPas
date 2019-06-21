/**
 * @file VerletClustersTraversalInterface.h
 * @author humig
 * @date 20.06.19
 */

#pragma once

namespace autopas {

/**
 * Some math operations for verlet clusters lists.
 */
class VerletClusterMaths {
 public:
  /**
   * the index type to access the particle cells
   */
  using index_t = std::size_t;
  /**
   * Converts 2d grid position to the index in the vector.
   * @param x x-position in grid
   * @param y y-position in grid
   * @param cellsPerDim the cells per dimension of the container
   * @return index in vector
   */
  static index_t index1D(const index_t x, const index_t y, std::array<index_t, 3> cellsPerDim) {
    return x + y * cellsPerDim[0];
  }
};

}  // namespace autopas
