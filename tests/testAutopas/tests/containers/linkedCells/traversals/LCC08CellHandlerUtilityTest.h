/**
 * @file LCC08CellHandlerUtilityTest.h
 * @author J. Schuhmacher
 * @date 22.07.24
 */

#pragma once

#include <gmock/gmock.h>

#include <array>
#include <ranges>
#include <vector>

#include "AutoPasTestBase.h"
#include "autopas/containers/linkedCells/traversals/LCC08CellHandlerUtility.h"

/**
 * Test cases checking the C08 cell offset generation and interaction-cell pair generation
 */
class LCC08CellHandlerUtilityTest : public AutoPasTestBase {
 public:
  LCC08CellHandlerUtilityTest() = default;

  ~LCC08CellHandlerUtilityTest() override = default;

 protected:
  /**
   * The cells per dimension utilized over the tests
   */
  constexpr static inline std::array<unsigned long, 3> CELLS_PER_DIMENSION{
      12,
      12,
      12,
  };

  /**
   * The cell length utilized over the tests (the interaction length is the variable to control overlap)
   */
  constexpr static inline std::array<double, 3> CELL_LENGTH{1, 1, 1};

  /**
   * We us this for floating point comparisons, where the precision does not really matter (just the first numbers
   * should be equal)
   */
  constexpr static inline double TEST_EPSILON{1e-3};

  /**
   * Converts a vector of cell offset pairs into a sorted vector of cell offset differences.
   * So, each pair is flattend to a single difference. This enables testing of the cell pairs
   * in an agnostic fahsion independent of sorting order, order of offsets inside the pair, and the concret
   * pairs (when only the correct interaction is required)
   * @tparam Mode one of c08CellPairs modes (not c04CellPairs, as two-dimensionl inputs aren't supported)
   * @param offsetPairs the vector of offsets pairs/ triplets
   * @return vector of offset-differences
   */
  template <autopas::LCC08CellHandlerUtility::C08OffsetMode Mode>
  static std::vector<unsigned long> transformAndSortOffsetPairs(
      const autopas::LCC08CellHandlerUtility::OffsetPairType<Mode> &offsetPairs) {
    static_assert(Mode != autopas::LCC08CellHandlerUtility::C08OffsetMode::c04NoSorting,
                  "A two dimensional vector is not supported!");
    std::vector<unsigned long> pairOffsetsDifferences;
    pairOffsetsDifferences.reserve(offsetPairs.size());
    std::transform(offsetPairs.begin(), offsetPairs.end(), std::back_inserter(pairOffsetsDifferences),
                   [](const auto &tuple) {
                     return std::abs(static_cast<int>(std::get<0>(tuple)) - static_cast<int>(std::get<1>(tuple)));
                   });
    std::sort(pairOffsetsDifferences.begin(), pairOffsetsDifferences.end());
    return pairOffsetsDifferences;
  }

  /**
   * Sorts a vector of cell offset triplets as well as the triplets themselves.
   * @tparam Mode one of c08CellPairs modes (not c04CellPairs, as two-dimensional inputs aren't supported)
   * @param offsetTriplets the vector of offsets triplets
   * @return sorted vector of offset-triplets
   */
  template <autopas::LCC08CellHandlerUtility::C08OffsetMode Mode>
  static std::vector<std::tuple<unsigned long, unsigned long, unsigned long>> sortOffsetTriplets(
      const autopas::LCC08CellHandlerUtility::OffsetTripletType<Mode> &offsetTriplets) {
    static_assert(Mode != autopas::LCC08CellHandlerUtility::C08OffsetMode::c04NoSorting,
                  "A two dimensional vector is not supported!");

    std::vector<std::tuple<unsigned long, unsigned long, unsigned long>> sortedTripletOffsets;
    sortedTripletOffsets.reserve(offsetTriplets.size());

    std::transform(
        offsetTriplets.begin(), offsetTriplets.end(), std::back_inserter(sortedTripletOffsets), [](const auto &tuple) {
          std::array<unsigned long, 3> tmpArray = {std::get<0>(tuple), std::get<1>(tuple), std::get<2>(tuple)};
          std::ranges::sort(tmpArray);
          return std::make_tuple(tmpArray[0], tmpArray[1], tmpArray[2]);
        });
    std::ranges::sort(sortedTripletOffsets);
    return sortedTripletOffsets;
  }

  static std::vector<std::tuple<long, long, long>> generateC18Triplets(long overlap, double interactionLength);
};
