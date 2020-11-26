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
 * Helper function for calculating the neighbor list length for the Verlet lists cells neighbor list.
 *
 * @param neighborLists
 * @param cellIndex the index of the current cell being processed
 * @return estimated load for current cell
 */
template <class Particle>
unsigned long neighborListLengthImpl(
    const typename autopas::VerletListsCellsHelpers<Particle>::NeighborListsType &neighborLists,
    unsigned long cellIndex) {
  unsigned long cellLoad = 0;
  for (auto &list : neighborLists[cellIndex]) {
    cellLoad += list.second.size();
  }

  return cellLoad;
}

/**
 * Helper function for calculating the neighbor list length for pairwise Verlet lists.
 *
 * @param neighborLists
 * @param cellIndex the index of the current cell being processed
 * @return estimated load for current cell
 */
template <class Particle>
unsigned long neighborListLengthImpl(
    const typename autopas::VerletListsCellsHelpers<Particle>::PairwiseNeighborListsType &neighborLists,
    unsigned long cellIndex) {
  unsigned long cellLoad = 0;
  for (auto &list : neighborLists[cellIndex]) {
    for (size_t index = 0; index < list.size(); index++) {
      cellLoad += list[index].second.size();
    }
  }

  return cellLoad;
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
template <class Particle, class NeighborList>
unsigned long neighborListLength(NeighborList &neighborLists, const std::array<unsigned long, 3> &cellsPerDimension,
                                 const std::array<unsigned long, 3> &lowerCorner,
                                 const std::array<unsigned long, 3> &upperCorner) {
  unsigned long sum = 0;
  for (unsigned long x = lowerCorner[0]; x <= upperCorner[0]; x++) {
    for (unsigned long y = lowerCorner[1]; y <= upperCorner[1]; y++) {
      for (unsigned long z = lowerCorner[2]; z <= upperCorner[2]; z++) {
        auto cellIndex = autopas::utils::ThreeDimensionalMapping::threeToOneD(x, y, z, cellsPerDimension);
        sum += neighborListLengthImpl<Particle>(neighborLists.getAoSNeighborList(), cellIndex);
      }
    }
  }
  return sum;
}

}  // namespace autopas::loadEstimators
