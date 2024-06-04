/**
 * @file VerletListsCellsHelpers.h
 * @author nguyen
 * @date 30.08.18
 */

#pragma once

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

namespace autopas::VerletListsCellsHelpers {
/**
 * Cell wise verlet lists for neighbors from all adjacent cells: For every cell, a vector of pairs.
 * Each pair maps a particle to a vector of its neighbors.
 *
 * @note From a content view, this is similar to an vector<unstructured_map<Particle*, std::vector<Particle *>>>.
 * However since we need to access all keys sequentially during the force computation this is faster,
 * even when the lookup of keys is slower.
 *
 * @tparam Particle
 */
template <class Particle>
using AllCellsNeighborListsType = std::vector<std::vector<std::pair<Particle *, std::vector<Particle *>>>>;

/**
 * Pairwise verlet lists: For every cell a vector, for every neighboring cell a vector of particle-neighborlist pairs.
 * Each pair maps a particle to a vector of its neighbor particles.
 * Cells<NeighboringCells<Particle,NeighborParticles>>
 *
 * @tparam Particle
 */
template <class Particle>
using PairwiseNeighborListsType = std::vector<std::vector<std::vector<std::pair<Particle *, std::vector<Particle *>>>>>;

/**
 * Indicates which build functor should be used for the generation of the neighbor list.
 * To be passed to the generator functors in the neighbor lists.
 */
enum class VLCBuildType {
  aosBuild,
  soaBuild,
};

/**
 * Helper Struct to bundle all information about base step offsets.
 */
struct BaseStepOffsets {
  /**
   * Offset from the base cell for cell 1.
   */
  int offset1;
  /**
   * Offset from the base cell for cell 2.
   */
  int offset2;
  /**
   * Estimation factor on the fraction of particles that will end up needing a neighbor list in the base cell.
   */
  double listSizeEstimateFactor;

  /**
   * Equality operator.
   * @param rhs
   * @return Equality of all members.
   */
  bool operator==(const BaseStepOffsets &rhs) const;

  /**
   * Inequality operator.
   * @param rhs
   * @return Inequality of at least one member.
   */
  bool operator!=(const BaseStepOffsets &rhs) const;
};

/**
 *
 * @param cellsPerDim
 * @return
 */
std::vector<BaseStepOffsets> buildC08BaseStep(const std::array<int, 3> &cellsPerDim);

/**
 * Simple heuristic to calculate the average number of particles per verlet list
 * assuming particles are evenly distributed in the domain box.
 *
 * @param numParticles
 * @param boxSize Size of the simulation box.
 * @param interactionLength Cutoff + skin.
 * @param correctionFactor Correction factor multiplied with the result.
 * @return numParticles * (list volume as fraction of box volume) * correctionFactor.
 */
size_t estimateListLength(size_t numParticles, const std::array<double, 3> &boxSize, double interactionLength,
                          double correctionFactor);

/**
 * Function to estimate the number of neighbor lists for one base step.
 *
 * @tparam Cells
 * @param baseCellIndex
 * @param useNewton3
 * @param cells
 * @param offsetsC08
 * @return
 */
template <class Cells>
size_t estimateNumLists(size_t baseCellIndex, bool useNewton3, Cells &cells,
                        const std::vector<BaseStepOffsets> &offsetsC08) {
  size_t estimate = 0;
  std::map<int, double> offsetFactors{};
  std::map<int, double> offsetFactorsNoN3{};
  for (const auto [offsetA, offsetB, factor] : offsetsC08) {
    offsetFactors[offsetA] = std::max(offsetFactors[offsetA], factor);
    offsetFactorsNoN3[offsetB] = std::max(offsetFactors[offsetB], factor);
  }
  for (const auto &[offset, factor] : offsetFactors) {
    estimate += cells[baseCellIndex + offset].size() * factor;
  }
  if (not useNewton3) {
    for (const auto &[offset, factor] : offsetFactorsNoN3) {
      estimate += cells[baseCellIndex + offset].size() * factor;
    }
  }
  return estimate;
};
}  // namespace autopas::VerletListsCellsHelpers
