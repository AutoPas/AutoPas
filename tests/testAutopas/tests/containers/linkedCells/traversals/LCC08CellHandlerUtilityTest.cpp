/**
 * @file LCC08CellHandlerUtilityTest.cpp
 * @author J. Schuhmacher
 * @date 22.07.24
 */

#include "LCC08CellHandlerUtilityTest.h"

#include "autopas/utils/StringUtils.h"

using testing::ContainerEq;
using testing::Eq;
using testing::Pointwise;

using autopas::LCC08CellHandlerUtility::C08OffsetMode;
using autopas::LCC08CellHandlerUtility::computePairwiseCellOffsetsC08;
using autopas::LCC08CellHandlerUtility::internal::computeRelativeCellOffsets;

/*
 * If the base cell has zero overlap with other cells dues to a overlap = interactionLength / cellLength = 0
 * ==> Expcected Cell Offsets for potential interactions, only self { 0 }
 */
TEST_F(LCC08CellHandlerUtilityTest, ComputeCellOffsetsC08Test_0x0x0) {
  constexpr std::array<int, 3> overlap{0, 0, 0};
  const std::vector<int> actualOffsets = computeRelativeCellOffsets(CELLS_PER_DIMENSION, overlap);
  const std::vector<int> expectedOffsets{0};
  ASSERT_THAT(actualOffsets, ContainerEq(expectedOffsets));
}

/*
 * If the base cell has an overlap of one with other cells
 * ==> Expected Cell Offets for potential interactions 2x2x2, with a 12 cells per dimension (for 1D index), it is
 *     {0, 144, 12, 156, 1, 145, 13, 157}
 */
TEST_F(LCC08CellHandlerUtilityTest, ComputeCellOffsetsC08Test_1x1x1) {
  constexpr std::array<int, 3> overlap{1, 1, 1};
  const std::vector<int> actualOffsets = computeRelativeCellOffsets(CELLS_PER_DIMENSION, overlap);
  const std::vector<int> expectedOffsets{0, 144, 12, 156, 1, 145, 13, 157};
  ASSERT_THAT(actualOffsets, ContainerEq(expectedOffsets));
}

/*
 * If the base cell has an overlap of two wither neighboring cells
 * ==> Expected Cell Offsets 3x3x3, i.e. an array of size 27, again with 12 cells per dimension (for 1D index)
 */
TEST_F(LCC08CellHandlerUtilityTest, ComputeCellOffsetsC08Test_2x2x2) {
  constexpr std::array<int, 3> overlap{2, 2, 2};
  const std::vector<int> actualOffsets = computeRelativeCellOffsets(CELLS_PER_DIMENSION, overlap);
  const std::vector<int> expectedOffsets{
      0,   144, 288, 12,  156, 300, 24,  168, 312, 1,   145, 289, 13,  157,
      301, 25,  169, 313, 2,   146, 290, 14,  158, 302, 26,  170, 314,
  };
  ASSERT_THAT(actualOffsets, ContainerEq(expectedOffsets));
}

/*
 * The given cell length and interaction length lead to an overlap of one.
 * Ergo, we have 2x2x2 cells and shall have 14 interaction pairs in total between the cells.
 * We test that the correct offset-distance are correct (details see below)
 */
TEST_F(LCC08CellHandlerUtilityTest, ComputePairwiseCellOffsetsC08Test_1x1x1) {
  constexpr double interactionLength{1.0};
  constexpr std::array<unsigned long, 14> expectedPairOffsetDifferences{
      0, 1, 11, 12, 13, 131, 132, 133, 143, 144, 145, 155, 156, 157,
  };

  const auto actualOffsetPairs = computePairwiseCellOffsetsC08<C08OffsetMode::c08CellPairsSorting>(
      CELLS_PER_DIMENSION, CELL_LENGTH, interactionLength);
  // Ensure the correct amount of interaction pairs
  ASSERT_EQ(actualOffsetPairs.size(), 14);

  // Transform the offset pairs to offset differences and sort them in-order
  // This way, the test is agnostic towards sorting order, offset-pair-order, the concret pairs
  // E.g. (0, 1) or (1, 0) would both valid. Here it is just tested as 1
  // E.g. (0, 1) or (12, 13) would both be valid (if the pattern is applied everywhere the same, it'ls like applying
  // the same interaction, but always shifted in y += 1). Here it is just tested as 1
  std::vector<unsigned long> actualPairOffsetsDiffercnes;
  std::transform(actualOffsetPairs.begin(), actualOffsetPairs.end(), std::back_inserter(actualPairOffsetsDiffercnes),
                 [](const auto &tuple) {
                   const auto &[offset1, offset2, sort] = tuple;
                   return std::abs(static_cast<int>(offset1) - static_cast<int>(offset2));
                 });
  std::sort(actualPairOffsetsDiffercnes.begin(), actualPairOffsetsDiffercnes.end());
  ASSERT_THAT(actualPairOffsetsDiffercnes, Pointwise(Eq(), expectedPairOffsetDifferences));
}

/*
 * The given cell length and interaction length lead to an overlap of two (3x3x3).
 * Ergo, we have 3x3x3 cells and shall have 63 interaction pairs in total between the cells.
 */
TEST_F(LCC08CellHandlerUtilityTest, ComputePairwiseCellOffsetsC08Test_2x2x2) {
  constexpr double interactionLength{2.0};
  constexpr std::array<unsigned long, 63> expectedPairOffsetDifferences{
      0,   1,   2,   10,  11,  12,  13,  14,  22,  23,  24,  25,  26,  118, 119, 120, 121, 122, 130, 131, 132,
      133, 134, 142, 143, 144, 145, 146, 154, 155, 156, 157, 158, 166, 167, 168, 169, 170, 262, 263, 264, 265,
      266, 274, 275, 276, 277, 278, 286, 287, 288, 289, 290, 298, 299, 300, 301, 302, 310, 311, 312, 313, 314,
  };

  const auto actualOffsetPairs = computePairwiseCellOffsetsC08<C08OffsetMode::c08CellPairsSorting>(
      CELLS_PER_DIMENSION, CELL_LENGTH, interactionLength);
  // Ensure the correct amount of interaction pairs
  ASSERT_EQ(actualOffsetPairs.size(), 63);

  // Flatten to offset differences (explaination, see ComputePairwiseCellOffsetsC08Test_1x1x1 test case)
  std::vector<unsigned long> actualPairOffsetsDiffercnes =
      transformAndSortOffsetPairs<C08OffsetMode::c08CellPairsSorting>(actualOffsetPairs);
  ASSERT_THAT(actualPairOffsetsDiffercnes, Pointwise(Eq(), expectedPairOffsetDifferences));
}

/*
 * The given cell length and interaction length lead to an overlap of three (4x4x4).
 * Ergo, we have 4x4x4 cells and shall have 168 interaction pairs in total between the cells.
 */
TEST_F(LCC08CellHandlerUtilityTest, ComputePairwiseCellOffsetsC08Test_3x3x3) {
  constexpr double interactionLength{3.0};
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

  const auto actualOffsetPairs = computePairwiseCellOffsetsC08<C08OffsetMode::c08CellPairsSorting>(
      CELLS_PER_DIMENSION, CELL_LENGTH, interactionLength);
  // Ensure the correct amount of interaction pairs
  EXPECT_EQ(actualOffsetPairs.size(), 168);

  // Flatten to offset differences (explaination, see ComputePairwiseCellOffsetsC08Test_1x1x1 test case)
  std::vector<unsigned long> actualPairOffsetsDiffercnes =
      transformAndSortOffsetPairs<C08OffsetMode::c08CellPairsSorting>(actualOffsetPairs);
  ASSERT_THAT(actualPairOffsetsDiffercnes, Pointwise(Eq(), expectedPairOffsetDifferences));
}