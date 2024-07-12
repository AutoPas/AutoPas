/**
 * @file LCCellHandler.h
 * @author J. Schuhmacher
 * @date 11.07.2024
 */

#pragma once

#include <array>
#include <optional>
#include <tuple>
#include <vector>

#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/Math.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas::internal {

/**
 * Type Alias for the C08 base step:
 * A triplet consisting of:
 *  - offset of first cell
 *  - offset of second cell
 *  - sorting direction a.k.a. the normalized vector between cell1 and cell2 connecting their centers. This is used
 *      in the {@link autopas::internal::CellFunctor} for AoS processing of Cell Pairs, ultimatley in
 *      {@link autopas::SortedCellView} for building a projection order of particles to early stop the processing.
 */
using C08CellOffsetPair = std::tuple<unsigned long, unsigned long, std::array<double, 3>>;

/**
 * Type Alias for the C08 base step.
 * This is a vector consisting of triplets of {@link autopas::internal::C08CellOffsetPair}
 */
using C08CellOffsetPairVector = std::vector<C08CellOffsetPair>;

/**
 * Type Alias for the C08 base step.
 * This is a vector consisting of cell indices required for the C08 base step in one-dimensional notation.
 */
using C08CellOffsetVector = std::vector<int>;

/**
 * Represents the interaction directions between cell pairs in the C08 base step
 */
enum class C08CellDirection : int {
  /** The origin or the base cell in the four-cell-square **/
  FRONT_LEFT,
  /** The cell to the back from the origin **/
  BACK_LEFT,
  /** The cell to the right from the origin **/
  FRONT_RIGHT,
  /** The cell diagonally to the origin **/
  BACK_RIGHT,
};

/**
 * Array containing all four enum values of {@link autopas::internal::C08CellDirection}.
 */
constexpr inline std::array<C08CellDirection, 4> ALL_DIRECTIONS{
    {C08CellDirection::FRONT_LEFT, C08CellDirection::BACK_LEFT, C08CellDirection::FRONT_RIGHT,
     C08CellDirection::BACK_RIGHT}};

/**
 * Calculates the offset for cell2 in the pair, given a direction from the base cell, a z index and the overlap.
 * @param overlap1 the overlap of the interacting cells (plus one)
 * @param direction one of the four directions
 * @param z the z index of the spatial dimension
 * @return offset of cell2
 */
constexpr int calculateOffset2(int overlap1, const C08CellDirection &direction, int z);

/**
 * Calculates the multipliers (zero or one) for a given direction.
 * This is used as a helper to calculate the `distVec` in {@link autopas::internal::computePairwiseCellOffsetsC08}.
 * @param direction one of the four directions
 * @return pair of multipliers: first one for x dimension, second one for y dimension
 */
constexpr std::pair<int, int> toDirectionMultiplier(const C08CellDirection &direction);

/**
 * Returns true if the cell-interaction in the given direction, with the given base cell coordinates and overlap
 * needs to be included into the pairwise cell offsets
 * @param direction the direction from base cell to inetracting cell
 * @param overlap the overlap (calculated from interactionLength divided by cellLength)
 * @param x the x offset of the base cell
 * @param y the y offset of the base cell
 * @param z the z offset of the base cell
 * @return true if the pair needs to be included
 */
constexpr bool includeCase(const C08CellDirection &direction, const std::array<int, 3> &overlap, int x, int y, int z);

/**
 * Computes the sorting direction between two cells from center of cell1 to center of cell2.
 * @param cellsPerDimension array containing the number of cells per dimension
 * @param distVec1 the vector to cell1 from the origin
 * @param offset2 the offset of cell2
 * @return a vector containing the sorting direction, i.e. the vector from center of cell1 to center of cell2
 */
std::array<double, 3> computeSortingDirection(const std::array<int, 3> &cellsPerDimension,
                                              const std::array<double, 3> &distVec1, int offset2);

/**
 * Computes the cell offsets for the C08 base step based on the given number of cells per dimension and overlap.
 * @param cellsPerDimension the number of cells per dimension
 * @param overlap the extent of overlap between neighboring cells
 * @return a vector of unsigned long values representing the computed cell offsets in 1D coordinates
 */
C08CellOffsetVector computeCellOffsetsC08(const std::array<unsigned long, 3> &cellsPerDimension,
                                          const std::array<int, 3> &overlap);

/**
 * Computes the cell pair offsets for the C08 base step and the normalized vector between pair of cell-centers,
 * which is later used for early stopping the evaluation of the pairwise cell interactions due to being out-of-reach.
 * @param cellsPerDimension the number of cells per dimension
 * @param cellLength the length of a cell in CellBlock3D.
 * @param interactionLength the interaction length consisting of cutoff + skin
 * @return vector of triplets consisting of: offset of first cell, offset of second cell, sorting direction
 *    (find details in {@link autopas::C08StepOffsetVector}
 */
C08CellOffsetPairVector computePairwiseCellOffsetsC08(const std::array<unsigned long, 3> &cellsPerDimension,
                                                      const std::array<double, 3> &cellLength,
                                                      double interactionLength);

}  // namespace autopas::internal
