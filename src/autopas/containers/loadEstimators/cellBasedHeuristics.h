/**
 * @file cellBasedHeuristics.h
 *
 * @date 24 Apr 2020
 * @author fischerv
 */

#pragma once

#include <array>
#include <vector>

#include "autopas/utils/ThreeDimensionalMapping.h"

/**
 * Collection of functions for estimating the load required to update a specific region within a containers
 */
namespace autopas::loadEstimators {

/**
 * Possible heuristics
 */
enum class CellBasedHeuristic {
  none = 0,
  squaredCellSize = 1,
};

/**
 * Sums up the squared number of particles for all cells within region
 *
 * @param cells
 * @param cellsPerDimension
 * @param lowerCorner lower boundary indices for region
 * @param upperCorner upper boundary indices for region
 * @return estimated load for given region
 */
template <class ParticleCell>
unsigned long squaredCellSize(const std::vector<ParticleCell> &cells,
                              const std::array<unsigned long, 3> &cellsPerDimension,
                              const std::array<unsigned long, 3> &lowerCorner,
                              const std::array<unsigned long, 3> &upperCorner) {
  unsigned long sum = 0;
  for (unsigned long x = lowerCorner[0]; x <= upperCorner[0]; x++) {
    for (unsigned long y = lowerCorner[1]; y <= upperCorner[1]; y++) {
      for (unsigned long z = lowerCorner[2]; z <= upperCorner[2]; z++) {
        auto load =
            cells[autopas::utils::ThreeDimensionalMapping::threeToOneD(x, y, z, cellsPerDimension)].numParticles();
        sum += load * load;
      }
    }
  }
  return sum;
}

/**
 * Estimate the load for the given region within a cell based container.
 *
 * @param heuristic
 * @param cells
 * @param cellsPerDimension
 * @param lowerCorner lower boundary indices for region
 * @param upperCorner upper boundary indices for region
 * @return estimated load for given region
 */
template <class ParticleCell>
unsigned long estimateCellBasedLoad(const CellBasedHeuristic &heuristic, const std::vector<ParticleCell> &cells,
                                    const std::array<unsigned long, 3> &cellsPerDimension,
                                    const std::array<unsigned long, 3> &lowerCorner,
                                    const std::array<unsigned long, 3> &upperCorner) {
  switch (heuristic) {
    case CellBasedHeuristic::none:
      return 1;
    case CellBasedHeuristic::squaredCellSize:
      return squaredCellSize(cells, cellsPerDimension, lowerCorner, upperCorner);
  }
  return 1;
}
}  // namespace autopas::loadEstimators
