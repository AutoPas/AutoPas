/**
 * @file LCC08CellHandlerUtilityTest.cpp
 * @author J. Schuhmacher
 * @date 22.07.24
 */

#include "LCC08CellHandlerUtilityTest.h"

using testing::ContainerEq;
using testing::DoubleNear;
using testing::Eq;
using testing::Pointwise;

using autopas::LCC08CellHandlerUtility::C08OffsetMode;
using autopas::LCC08CellHandlerUtility::computePairwiseCellOffsetsC08;

/*
 * The given cell length and interaction length lead to an overlap of one.
 * Ergo, we have 2x2x2 cells and shall have 14 interaction pairs in total between the cells.
 * We test that the correct offset-distance are correct (details see below)
 */
TEST_F(LCC08CellHandlerUtilityTest, ComputePairwiseCellOffsetsC08Test_1x1x1) {
  constexpr double interactionLength{1.0};
  // Calculated by hand & checked with the legacy implementation (Commit ID: d560a7075)
  constexpr std::array<unsigned long, 14> expectedPairOffsetDifferences{
      0, 1, 11, 12, 13, 131, 132, 133, 143, 144, 145, 155, 156, 157,
  };
  constexpr std::array<std::pair<unsigned long, unsigned long>, 14> expectedPairOffsets{
      std::make_pair(0, 0),    std::make_pair(144, 0),  std::make_pair(156, 0), std::make_pair(145, 0),
      std::make_pair(157, 0),  std::make_pair(0, 12),   std::make_pair(1, 12),  std::make_pair(144, 12),
      std::make_pair(145, 12), std::make_pair(0, 1),    std::make_pair(144, 1), std::make_pair(156, 1),
      std::make_pair(0, 13),   std::make_pair(144, 13),
  };

  const auto actualOffsetPairs =
      computePairwiseCellOffsetsC08<C08OffsetMode::noSorting>(CELLS_PER_DIMENSION, CELL_LENGTH, interactionLength);
  // Ensure the correct amount of interaction pairs
  ASSERT_EQ(actualOffsetPairs.size(), 14);

  // Transform the offset pairs to offset differences and sort them in-order
  // This way, the test is agnostic towards sorting order, offset-pair-order, the concret pairs
  // E.g. (0, 1) or (1, 0) would both valid. Here it is just tested as 1
  // E.g. (0, 1) or (12, 13) would both be valid (if the pattern is applied everywhere the same, it'ls like applying
  // the same interaction, but always shifted in y += 1). Here it is just tested as 1
  std::vector<unsigned long> actualPairOffsetsDiffercnes =
      transformAndSortOffsetPairs<C08OffsetMode::noSorting>(actualOffsetPairs);
  ASSERT_THAT(actualPairOffsetsDiffercnes, Pointwise(Eq(), expectedPairOffsetDifferences));

  // This test case is more senstive, it will fail in case the ordering of the output-pairs is wrong
  // If this fails, your implementation is not necessarily wrong - but different
  ASSERT_THAT(actualOffsetPairs, Pointwise(Eq(), expectedPairOffsets));
}

TEST_F(LCC08CellHandlerUtilityTest, ComputePairwiseCellOffsetsC08Test_1x1x1_Sorting) {
  constexpr double interactionLength{1.0};
  // Computed with the legacy implementation (assuming it has been correct), Commit ID: d560a7075
  constexpr std::array<std::array<double, 3>, 14> expectedSortingVectors{{
      {0.57735, 0.57735, 0.57735},
      {0, 0, -1},
      {0, -0.707107, -0.707107},
      {-0.707107, 0, -0.707107},
      {-0.57735, -0.57735, -0.57735},
      {0, 1, 0},
      {-0.707107, 0.707107, 0},
      {0, 0.707107, -0.707107},
      {-0.57735, 0.57735, -0.57735},
      {1, 0, 0},
      {0.707107, 0, -0.707107},
      {0.57735, -0.57735, -0.57735},
      {0.707107, 0.707107, 0},
      {0.57735, 0.57735, -0.57735},
  }};

  const auto actualOffsetTriplets =
      computePairwiseCellOffsetsC08<C08OffsetMode::sorting>(CELLS_PER_DIMENSION, CELL_LENGTH, interactionLength);
  // Ensure the correct amount of interaction pairs
  ASSERT_EQ(actualOffsetTriplets.size(), 14);

  std::vector<std::array<double, 3>> actualSortingVectors{};
  actualSortingVectors.reserve(14);
  std::transform(actualOffsetTriplets.begin(), actualOffsetTriplets.end(), std::back_inserter(actualSortingVectors),
                 [](const auto &tuple) { return std::get<2>(tuple); });
  for (size_t i = 0; i < expectedSortingVectors.size(); ++i) {
    ASSERT_THAT(actualSortingVectors[i], Pointwise(DoubleNear(TEST_EPSILON), expectedSortingVectors[i]));
  }
}

/*
 * The given cell length and interaction length lead to an overlap of two (3x3x3).
 * Ergo, we have 3x3x3 cells and shall have 63 interaction pairs in total between the cells.
 */
TEST_F(LCC08CellHandlerUtilityTest, ComputePairwiseCellOffsetsC08Test_2x2x2) {
  constexpr double interactionLength{2.0};
  // Checked by Hand & Computed with the legacy implementation (Commit ID: d560a7075)
  constexpr std::array<unsigned long, 63> expectedPairOffsetDifferences{
      0,   1,   2,   10,  11,  12,  13,  14,  22,  23,  24,  25,  26,  118, 119, 120, 121, 122, 130, 131, 132,
      133, 134, 142, 143, 144, 145, 146, 154, 155, 156, 157, 158, 166, 167, 168, 169, 170, 262, 263, 264, 265,
      266, 274, 275, 276, 277, 278, 286, 287, 288, 289, 290, 298, 299, 300, 301, 302, 310, 311, 312, 313, 314,
  };

  const auto actualOffsetPairs =
      computePairwiseCellOffsetsC08<C08OffsetMode::sorting>(CELLS_PER_DIMENSION, CELL_LENGTH, interactionLength);
  // Ensure the correct amount of interaction pairs
  ASSERT_EQ(actualOffsetPairs.size(), 63);

  // Flatten to offset differences (explaination, see ComputePairwiseCellOffsetsC08Test_1x1x1 test case)
  std::vector<unsigned long> actualPairOffsetsDiffercnes =
      transformAndSortOffsetPairs<C08OffsetMode::sorting>(actualOffsetPairs);
  ASSERT_THAT(actualPairOffsetsDiffercnes, Pointwise(Eq(), expectedPairOffsetDifferences));
}

TEST_F(LCC08CellHandlerUtilityTest, ComputePairwiseCellOffsetsC08Test_2x2x2_Sorting) {
  constexpr double interactionLength{2.0};
  // Computed with the legacy implementation (assuming it has been correct), Commit ID: d560a7075
  constexpr std::array<std::array<double, 3>, 63> expectedSortingVectors{{
      {0.57735, 0.57735, 0.57735},
      {0, 0, -1},
      {0, -0.894427, -0.447214},
      {-0.894427, 0, -0.447214},
      {-0.666667, -0.666667, -0.333333},
      {0, 0, -1},
      {0, -0.707107, -0.707107},
      {-0.707107, 0, -0.707107},
      {-0.57735, -0.57735, -0.57735},
      {0, 1, 0},
      {-0.894427, 0.447214, 0},
      {0, 0.707107, -0.707107},
      {0, -0.707107, -0.707107},
      {-0.816497, 0.408248, -0.408248},
      {-0.816497, -0.408248, -0.408248},
      {0, 0.447214, -0.894427},
      {0, -0.447214, -0.894427},
      {-0.666667, 0.333333, -0.666667},
      {-0.666667, -0.333333, -0.666667},
      {0, 1, 0},
      {-0.707107, 0.707107, 0},
      {0, 0.894427, -0.447214},
      {-0.666667, 0.666667, -0.333333},
      {0, 0.707107, -0.707107},
      {-0.57735, 0.57735, -0.57735},
      {1, 0, 0},
      {0.707107, 0, -0.707107},
      {0.408248, -0.816497, -0.408248},
      {-0.707107, 0, -0.707107},
      {-0.408248, -0.816497, -0.408248},
      {0.447214, 0, -0.894427},
      {0.333333, -0.666667, -0.666667},
      {-0.447214, 0, -0.894427},
      {-0.333333, -0.666667, -0.666667},
      {0.707107, 0.707107, 0},
      {-0.707107, 0.707107, 0},
      {0.57735, 0.57735, -0.57735},
      {0.57735, -0.57735, -0.57735},
      {-0.57735, 0.57735, -0.57735},
      {-0.57735, -0.57735, -0.57735},
      {0.408248, 0.408248, -0.816497},
      {0.408248, -0.408248, -0.816497},
      {-0.408248, 0.408248, -0.816497},
      {-0.408248, -0.408248, -0.816497},
      {0.447214, 0.894427, 0},
      {-0.447214, 0.894427, 0},
      {0.408248, 0.816497, -0.408248},
      {-0.408248, 0.816497, -0.408248},
      {0.333333, 0.666667, -0.666667},
      {-0.333333, 0.666667, -0.666667},
      {1, 0, 0},
      {0.894427, 0, -0.447214},
      {0.666667, -0.666667, -0.333333},
      {0.707107, 0, -0.707107},
      {0.57735, -0.57735, -0.57735},
      {0.894427, 0.447214, 0},
      {0.816497, 0.408248, -0.408248},
      {0.816497, -0.408248, -0.408248},
      {0.666667, 0.333333, -0.666667},
      {0.666667, -0.333333, -0.666667},
      {0.707107, 0.707107, 0},
      {0.666667, 0.666667, -0.333333},
      {0.57735, 0.57735, -0.57735},
  }};

  const auto actualOffsetTriplets =
      computePairwiseCellOffsetsC08<C08OffsetMode::sorting>(CELLS_PER_DIMENSION, CELL_LENGTH, interactionLength);
  // Ensure the correct amount of interaction pairs
  ASSERT_EQ(actualOffsetTriplets.size(), 63);

  std::vector<std::array<double, 3>> actualSortingVectors{};
  actualSortingVectors.reserve(63);
  std::transform(actualOffsetTriplets.begin(), actualOffsetTriplets.end(), std::back_inserter(actualSortingVectors),
                 [](const auto &tuple) { return std::get<2>(tuple); });
  for (size_t i = 0; i < expectedSortingVectors.size(); ++i) {
    ASSERT_THAT(actualSortingVectors[i], Pointwise(DoubleNear(TEST_EPSILON), expectedSortingVectors[i]));
  }
}

/*
 * The given cell length and interaction length lead to an overlap of three (4x4x4).
 * Ergo, we have 4x4x4 cells and shall have 168 interaction pairs in total between the cells.
 */
TEST_F(LCC08CellHandlerUtilityTest, ComputePairwiseCellOffsetsC08Test_3x3x3) {
  constexpr double interactionLength{3.0};
  // Calculated with the legacy implementation (Commit ID: d560a7075)
  constexpr std::array<unsigned long, 168> expectedPairOffsetDifferences{
      0,   1,   2,   3,   9,   10,  11,  12,  13,  14,  15,  21,  22,  23,  24,  25,  26,  27,  33,  34,  35,
      36,  37,  38,  39,  105, 106, 107, 108, 109, 110, 111, 117, 118, 119, 120, 121, 122, 123, 129, 130, 131,
      132, 133, 134, 135, 141, 142, 143, 144, 145, 146, 147, 153, 154, 155, 156, 157, 158, 159, 165, 166, 167,
      168, 169, 170, 171, 177, 178, 179, 180, 181, 182, 183, 249, 250, 251, 252, 253, 254, 255, 261, 262, 263,
      264, 265, 266, 267, 273, 274, 275, 276, 277, 278, 279, 285, 286, 287, 288, 289, 290, 291, 297, 298, 299,
      300, 301, 302, 303, 309, 310, 311, 312, 313, 314, 315, 321, 322, 323, 324, 325, 326, 327, 394, 395, 396,
      397, 398, 405, 406, 407, 408, 409, 410, 411, 417, 418, 419, 420, 421, 422, 423, 429, 430, 431, 432, 433,
      434, 435, 441, 442, 443, 444, 445, 446, 447, 453, 454, 455, 456, 457, 458, 459, 466, 467, 468, 469, 470,
  };

  const auto actualOffsetPairs =
      computePairwiseCellOffsetsC08<C08OffsetMode::sorting>(CELLS_PER_DIMENSION, CELL_LENGTH, interactionLength);
  // Ensure the correct amount of interaction pairs
  EXPECT_EQ(actualOffsetPairs.size(), 168);

  // Flatten to offset differences (explaination, see ComputePairwiseCellOffsetsC08Test_1x1x1 test case)
  std::vector<unsigned long> actualPairOffsetsDiffercnes =
      transformAndSortOffsetPairs<C08OffsetMode::sorting>(actualOffsetPairs);
  ASSERT_THAT(actualPairOffsetsDiffercnes, Pointwise(Eq(), expectedPairOffsetDifferences));
}