/**
 * @file LoadEstimators.h
 *
 * @date 24 Apr 2020
 * @author fischerv
 */

#pragma once

#include <array>
#include <vector>

#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsHelpers.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

/**
 * Collection of functions for estimating the load required to update a specific region within a containers.
 */
namespace autopas::loadEstimators {

/**
 * Sums up the squared number of particles for all cells within region.
 *
 * @param cells
 * @param cellsPerDimension
 * @param lowerCorner lower boundary indices for region
 * @param upperCorner upper boundary indices for region
 * @return estimated load for given region
 */
template <class ParticleCell>
unsigned long squaredParticlesPerCell(const std::vector<ParticleCell> &cells,
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
 * Sums up the lengths of the verlet neighbor lists of all particles within region.
 *
 * @param neighborLists
 * @param cellsPerDimension
 * @param lowerCorner lower boundary indices for region
 * @param upperCorner upper boundary indices for region
 * @return estimated load for given region
 */
template <class Particle>
unsigned long neighborListLength(
    const typename autopas::VerletListsCellsHelpers<Particle>::NeighborListsType &neighborLists,
    const std::array<unsigned long, 3> &cellsPerDimension, const std::array<unsigned long, 3> &lowerCorner,
    const std::array<unsigned long, 3> &upperCorner) {
  unsigned long sum = 0;
  for (unsigned long x = lowerCorner[0]; x <= upperCorner[0]; x++) {
    for (unsigned long y = lowerCorner[1]; y <= upperCorner[1]; y++) {
      for (unsigned long z = lowerCorner[2]; z <= upperCorner[2]; z++) {
        auto cellIndex = autopas::utils::ThreeDimensionalMapping::threeToOneD(x, y, z, cellsPerDimension);
        unsigned long cellLoad = 0;
        for (auto &list : neighborLists[cellIndex]) {
          cellLoad += list.second.size();
        }
        sum += cellLoad;
      }
    }
  }
  return sum;
}

template <class Particle>
unsigned long neighborListLength(
    const typename autopas::VerletListsCellsHelpers<Particle>::PairwiseNeighborListsType &neighborLists,
    const std::array<unsigned long, 3> &cellsPerDimension, const std::array<unsigned long, 3> &lowerCorner,
    const std::array<unsigned long, 3> &upperCorner) {
  unsigned long sum = 0;
  for (unsigned long x = lowerCorner[0]; x <= upperCorner[0]; x++) {
    for (unsigned long y = lowerCorner[1]; y <= upperCorner[1]; y++) {
      for (unsigned long z = lowerCorner[2]; z <= upperCorner[2]; z++) {
        auto cellIndex = autopas::utils::ThreeDimensionalMapping::threeToOneD(x, y, z, cellsPerDimension);
        unsigned long cellLoad = 0;
        for (auto &list : neighborLists[cellIndex]) {
          for (size_t index = 0; index < list.size(); index++) {
            cellLoad += list[index].second.size();
          }
        }
        sum += cellLoad;
      }
    }
  }
  return sum;
}

}  // namespace autopas::loadEstimators
