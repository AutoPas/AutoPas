/**
 * @file LCC08CellHandlerUtility.h
 * @author J. Schuhmacher
 * @date 11.07.2024
 */

#pragma once

#include <array>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/Math.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas::LCC08CellHandlerUtility {

/**
 * Type Alias for the C08 base step containing cell offsets. An offset is the distance from a base cell to another cell
 * in one-dimensional coordinates.
 * A triplet consisting of:
 *  - offset of first cell
 *  - offset of second cell
 *  - sorting direction a.k.a. the normalized vector between cell1 and cell2 connecting their centers. This is used
 *      in the CellFunctor for AoS processing of Cell Pairs, ultimately in
 *      autopas::SortedCellView for building a projection order of particles to early stop the processing.
 */
using OffsetPairSorting = std::tuple<unsigned long, unsigned long, std::array<double, 3>>;

/**
 * Type Alias for the C08 base step containing cell offsets for cell triplets. An offset is the distance from a base
 * cell to another cell in one-dimensional coordinates. It is a tuple consisting of:
 *  - offset of the first cell
 *  - offset of the second cell
 *  - offset of the third cell
 *  - sorting direction a.k.a. the normalized vector between cell1 and cell2 connecting their centers. This is used
 *      in the CellFunctor for AoS processing of Cell Triplets, ultimately in
 *      autopas::SortedCellView for building a projection order of particles to early stop the processing.
 */
using OffsetTripletSorting = std::tuple<unsigned long, unsigned long, unsigned long, std::array<double, 3>>;

/**
 * An offset is the distance from a base cell to another cell
 * in one dimensional coordinates. Hence, this is a pair consisting of.
 *  - offset of first cell
 *  - offset of second cell
 */
using OffsetPair = std::pair<unsigned long, unsigned long>;

/**
 * An offset is the distance from a base cell to another cell
 * in one-dimensional coordinates. Hence, this is a triplet consisting of.
 *  - offset of the first cell
 *  - offset of the second cell
 *  - offset of the third cell
 */
using OffsetTriplet = std::tuple<unsigned long, unsigned long, unsigned long>;

/**
 * A vector of OffsetPairs
 */
using OffsetPairVector = std::vector<OffsetPair>;

/**
 * A vector of OffsetTriplets
 */
using OffsetTripletVector = std::vector<OffsetTriplet>;

/**
 * Compile Time Modes for the function autopas::LCC08CellHandlerUtility::computePairwiseCellOffsetsC08
 *
 * @note In case of a new mode, this also requires the explciit instantation of the new template
 * in LCC08CellHandlerUtility.cpp and a modificaqtion to {@link OffsetPairType}
 */
enum class C08OffsetMode {
  /** Returns the C08 base step cell pairs without sorting */
  noSorting = 0,
  /** Returns the C08 base step cell pairs with sorting directions (for SortedView projection) */
  sorting = 1,
  /** Returns the C08 base step cell pairs adapted to C04, i.e. two-dimensions resolved on X-axis */
  c04NoSorting = 2,
};

/**
 * Template Magic Parameter Alias, which links the types {@link OffsetPairSorting}, {@link OffsetPair} and {@link
 * OffsetPairVector}
 */
template <C08OffsetMode Mode>
using OffsetPairType = std::vector<
    std::conditional_t<Mode == C08OffsetMode::sorting, OffsetPairSorting,
                       std::conditional_t<Mode == C08OffsetMode::c04NoSorting, OffsetPairVector, OffsetPair>>>;

/**
 * Template Magic Parameter Alias, which links the types {@link OffsetTripletSorting}, {@link OffsetTriplet} and {@link
 * OffsetTripletVector}
 */
template <C08OffsetMode Mode>
using OffsetTripletType = std::vector<
    std::conditional_t<Mode == C08OffsetMode::sorting, OffsetTripletSorting,
                       std::conditional_t<Mode == C08OffsetMode::c04NoSorting, OffsetTripletVector, OffsetTriplet>>>;

namespace internal {

/**
 * Represents the interaction directions between cell pairs in the C08 base step
 */
enum class C08CellDirection : int {
  /** The origin or the base cell in the four-cell-square **/
  frontLeft,
  /** The cell to the back from the origin **/
  backLeft,
  /** The cell to the right from the origin **/
  frontRight,
  /** The cell diagonally to the origin **/
  backRight,
};

/**
 * Array containing all four enum values of {@link C08CellDirection}.
 */
constexpr inline std::array<C08CellDirection, 4> ALL_DIRECTIONS{
    {C08CellDirection::frontLeft, C08CellDirection::backLeft, C08CellDirection::frontRight,
     C08CellDirection::backRight}};

/**
 * Error message thrown in case C08CellDirection was extended but the extension was not included in the switch-statement
 * @todo c++20: make this an std::string
 */
constexpr inline char ENUM_EXTENSION_EXCEPTION[]{
    "Enum C08CellDirection was extended, but its assciated switch-case statements was not!"};

/**
 * Helper function for autopas::LCC08CellHandlerUtility::computePairwiseCellOffsetsC08.
 * This function basically translates a direction, like backLeft to the corresponding vector pointing
 * towards this cell relativly starting from the base cell. We treat frontLeft as base cell
 * So, e.g. frontLeft --> (0, 0) since we are good
 * Alternativley, e.g. backLeft --> (0, 1) pointing in positive y direction
 *
 * The values might need to be scaled given the overlap (so e.g. pointing over multiple cells into a direction)
 * @param direction one of the four directions
 * @return pair of multipliers: first one for x dimension, second one for y dimension
 */
constexpr std::pair<int, int> toMaskXY(const C08CellDirection &direction);

/**
 * Returns true if the cell-interaction in the given direction, with the given base cell coordinates (x, y, z) and
 * overlap needs to be included into the pairwise cell offsets. This function basically masks certain cell combinations,
 * which do not need to be included given their relative position/ or given that another pair already includes
 * the interaction.
 * @param direction the direction from base cell to inetracting cell
 * @param overlap the overlap (calculated from interactionLength divided by cellLength)
 * @param x the x offset of the base cell
 * @param y the y offset of the base cell
 * @param z the z offset of the base cell
 * @return true if the pair needs to be included
 */
constexpr bool includeCellPair(const C08CellDirection &direction, const std::array<int, 3> &overlap, int x, int y,
                               int z);

/**
 * Computes the sorting direction between two cells from center of cell1 to center of cell2 using the 3D indices
 * of the cells while incoperating the cellLength (required in case of less regular cuboid cells)
 * @param offset1Vector the cartesianOffset of cell1
 * @param offset2Vector the cartesianOffset of cell2
 * @param cellLength the cell length in all three dimensions
 * @return normalized vector containing the sorting direction, i.e. the vector from cell1 to cell2
 */
std::array<double, 3> computeSortingDirection(const std::array<double, 3> &offset1Vector,
                                              const std::array<double, 3> &offset2Vector,
                                              const std::array<double, 3> &cellLength);

}  // namespace internal

/**
 * Computes the cell pair offsets for the C08 base step and the normalized vector between pair of cell-centers,
 * which is later used for early stopping the evaluation of the pairwise cell interactions due to being out-of-reach.
 * @tparam Mode Determines the concret return type (see {@link C08OffsetMode}
 * @param cellsPerDimension the number of cells per dimension
 * @param cellLength the length of a cell in CellBlock3D.
 * @param interactionLength the interaction length consisting of cutoff + skin
 * @return depending on template parameters a
 *  - vector containing cell offset pairs
 *  - vector containg cell offsets + sorting/ vector between cell centers triplets
 *  - vector of vector containing cell offsets (pre-sorted after X dimension)
 */
template <C08OffsetMode Mode>
OffsetPairType<Mode> computePairwiseCellOffsetsC08(const std::array<unsigned long, 3> &cellsPerDimension,
                                                   const std::array<double, 3> &cellLength, double interactionLength);

/**
 * Computes the cell triplet offsets for the C08 base step.
 * If the Mode is `sorting`, a normalized vector connecting the base cell and the second cell is additionally returned.
 * @tparam Mode Determines the concrete return type (see {@link C08OffsetMode}
 * @param cellsPerDimension the number of cells per dimension
 * @param cellLength the length of a cell in CellBlock3D.
 * @param interactionLength the interaction length consisting of cutoff + skin
 * @return depending on template parameters a
 *  - vector containing cell offset triplets
 *  - vector containing cell offsets + sorting/vector between cell centers of the base cell and second cell
 *  - vector of vector containing cell offsets (pre-sorted after X dimension)
 */
template <C08OffsetMode Mode>
OffsetTripletType<Mode> computeTriwiseCellOffsetsC08(const std::array<unsigned long, 3> &cellsPerDimension,
                                                     const std::array<double, 3> &cellLength, double interactionLength);

/**
 * @copydoc autopas::LCC08CellHandlerUtility::computeTriwiseCellOffsetsC08()
 * @note This method uses another method to compute cell triplet offsets, which is faster for smaller cell sizes.
 * Currently, this function is not used/tested.
 */
template <C08OffsetMode Mode>
OffsetTripletType<Mode> computeTriwiseCellOffsetsC08Optimized(const std::array<unsigned long, 3> &cellsPerDimension,
                                                              const std::array<double, 3> &cellLength,
                                                              double interactionLength);

}  // namespace autopas::LCC08CellHandlerUtility
