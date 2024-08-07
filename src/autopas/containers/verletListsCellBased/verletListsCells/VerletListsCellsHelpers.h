/**
 * @file VerletListsCellsHelpers.h
 * @author nguyen
 * @date 30.08.18
 */

#pragma once

#include <array>
#include <cstddef>
#include <unordered_map>
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
 * Builds the list of offsets from the base cell for the c08 base step.
 * A offset pair are two cell indices relative to the base cell index that have to interact.
 *
 * The third tuple entry is an estimation factor on the fraction of particles that will end up needing a neighbor list
 * in the base cell.
 * This depends on the relative position of the two cells: Whether they are the same, or share a face, edge, or corner.
 *
 * @note This function is implemented for CSF==1, meaning preallocation weight factors are intended for CSF==1.
 *       For CSF > 1 the implementation works but the factors might be suboptimal.
 *       For CSF < 1 more factors (and different) are needed if more neighbor cells interact with the base cell.
 *
 * @param cellsPerDim Number of cells per dimension including halo.
 * @return Vector of tuples<offset1, offset2, listEstimateFactor>
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
 * @param baseCellIndex Cell index for which the estimate is made.
 * @param useNewton3 Whether or not the traversal that uses the lists employs Newton3.
 * @param cells Reference to the vector of cells.
 * @param offsetsC08 Vector of BaseStepOffsets.
 * @return An estimate of the number of lists that will be needed in the base cell.
 */
template <class Cells>
size_t estimateNumLists(size_t baseCellIndex, bool useNewton3, const Cells &cells,
                        const std::vector<BaseStepOffsets> &offsetsC08) {
  // First for every cell, find its biggest factor.
  // Meaning, find out what is the closest interaction type this cell is involved in.
  std::unordered_map<int, double> offsetFactors{};
  std::unordered_map<int, double> offsetFactorsNoN3{};
  for (const auto [offsetA, offsetB, factor] : offsetsC08) {
    offsetFactors[offsetA] = std::max(offsetFactors[offsetA], factor);
    offsetFactorsNoN3[offsetB] = std::max(offsetFactors[offsetB], factor);
  }
  // The estimate is constructed by summing the involved cells' sizes weighted by their factors.
  size_t estimate = 0;
  for (const auto &[offset, factor] : offsetFactors) {
    estimate += cells[baseCellIndex + offset].size() * factor;
  }
  // For the non Newton3 case, lists have to be created for particles that otherwise would already be covered.
  if (not useNewton3) {
    for (const auto &[offset, factor] : offsetFactorsNoN3) {
      estimate += cells[baseCellIndex + offset].size() * factor;
    }
  }
  return estimate;
};
}  // namespace autopas::VerletListsCellsHelpers
