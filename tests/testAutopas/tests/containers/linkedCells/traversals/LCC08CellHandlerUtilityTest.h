/**
 * @file LCC08CellHandlerUtilityTest.h
 * @author J. Schuhmacher
 * @date 22.07.24
 */

#pragma once

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <array>
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
    std::vector<unsigned long> pairOffsetsDiffercnes;
    pairOffsetsDiffercnes.reserve(offsetPairs.size());
    std::transform(offsetPairs.begin(), offsetPairs.end(), std::back_inserter(pairOffsetsDiffercnes),
                   [](const auto &tuple) {
                     return std::abs(static_cast<int>(std::get<0>(tuple)) - static_cast<int>(std::get<1>(tuple)));
                   });
    std::sort(pairOffsetsDiffercnes.begin(), pairOffsetsDiffercnes.end());
    return pairOffsetsDiffercnes;
  }
};
