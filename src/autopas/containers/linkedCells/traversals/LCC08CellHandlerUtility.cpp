/**
 * @file LCC08CellHandlerUtility.cpp
 * @author J. Schuhmacher
 * @date 11.07.2024
 */

#include "LCC08CellHandlerUtility.h"

namespace autopas::LCC08CellHandlerUtility {

namespace internal {

constexpr int calculateCell2Index(int overlap1, const C08CellDirection &direction, int z) {
  switch (direction) {
    case C08CellDirection::frontLeft:
      return z;
    case C08CellDirection::backLeft:
      return overlap1 * (overlap1 - 1) + z;
    case C08CellDirection::frontRight:
      return overlap1 * overlap1 * (overlap1 - 1) + z;
    case C08CellDirection::backRight:
      return overlap1 * overlap1 * overlap1 - overlap1 + z;
  }
  throw std::runtime_error(ENUM_EXTENSION_EXCEPTION);
}

constexpr std::pair<int, int> toMask(const C08CellDirection &direction) {
  switch (direction) {
    case C08CellDirection::frontLeft:
      return std::make_pair(0, 0);
    case C08CellDirection::backLeft:
      return std::make_pair(0, 1);
    case C08CellDirection::frontRight:
      return std::make_pair(1, 0);
    case C08CellDirection::backRight:
      return std::make_pair(1, 1);
  }
  throw std::runtime_error(ENUM_EXTENSION_EXCEPTION);
}

constexpr bool includeCellPair(const C08CellDirection &direction, const std::array<int, 3> &overlap, int x, int y,
                               int z) {
  switch (direction) {
    case C08CellDirection::backLeft:
      return y != overlap[1] and z != 0;
    case C08CellDirection::frontRight:
      return x != overlap[0] and (y != 0 or z != 0);
    case C08CellDirection::backRight:
      return y != overlap[1] and x != overlap[0] and z != 0;
    case C08CellDirection::frontLeft:
      return true;
  }
  throw std::runtime_error(ENUM_EXTENSION_EXCEPTION);
}

std::vector<int> computeRelativeCellOffsets(const std::array<unsigned long, 3> &cellsPerDimension,
                                            const std::array<int, 3> &overlap) {
  using utils::ArrayMath::addScalar;
  using utils::ArrayMath::prod;
  using utils::ArrayUtils::static_cast_copy_array;
  using utils::ThreeDimensionalMapping::threeToOneD;
  std::vector<int> cellOffsets{};
  cellOffsets.reserve(prod(addScalar(overlap, 1)));

  // Iterate over the three spatial dimension and calculate the possible offsets for cell-inetractions in
  // the C08 base step (filtering & pairing happens later)
  for (int x = 0ul; x <= overlap[0]; ++x) {
    for (int y = 0ul; y <= overlap[1]; ++y) {
      for (int z = 0ul; z <= overlap[2]; ++z) {
        cellOffsets.push_back(threeToOneD(x, y, z, static_cast_copy_array<int>(cellsPerDimension)));
      }
    }
  }
  return cellOffsets;
}

std::array<double, 3> computeSortingDirection(const std::array<int, 3> &cellsPerDimension,
                                              const std::array<double, 3> &distVec1, int offset2) {
  using namespace autopas::utils::ArrayMath::literals;
  // First cell represented through distVec1, second cell represented through BaseCellVec (and offset2)
  const auto baseCellVec = utils::ArrayUtils::static_cast_copy_array<double>(
      utils::ThreeDimensionalMapping::oneToThreeD(offset2, cellsPerDimension));

  // In case the sorting direction is 0, 0, 0 ==> fix to 1, 1, 1
  std::array<double, 3> sortDir = distVec1 - baseCellVec;
  if (std::all_of(sortDir.begin(), sortDir.end(), [](const auto &val) { return val == 0; })) {
    sortDir = {1., 1., 1.};
  }

  // Normalize and return
  return utils::ArrayMath::normalize(sortDir);
}

}  // namespace internal

template <C08OffsetMode Mode>
OffsetPairType<Mode> computePairwiseCellOffsetsC08(const std::array<unsigned long, 3> &cellsPerDimension,
                                                   const std::array<double, 3> &cellLength, double interactionLength) {
  using namespace autopas::utils::ArrayMath::literals;
  using namespace internal;
  using utils::ArrayMath::ceilToInt;
  using utils::ArrayUtils::static_cast_copy_array;
  using utils::ThreeDimensionalMapping::threeToOneD;

  // Output of the function: Vector of pairwise cell indices & projection axis
  OffsetPairType<Mode> resultOffsetsC08{};

  // The overlap with interacting cells (see autopas::CBasedTraversal constructor)
  const std::array<int, 3> overlap{ceilToInt(interactionLength / cellLength)};
  // The offsets to the cells in 1D coordinates (in case overlap = [1, 1, 1] ==> cellOffsets contains 2x2x2 = 8 values)
  const std::vector<int> cellOffsets = computeRelativeCellOffsets(cellsPerDimension, overlap);
  const std::array<int, 3> &cellsPerDimIntegral = static_cast_copy_array<int>(cellsPerDimension);

  // Small constants used multiple times in the code below
  // Note, Legacy: With ov1 --> No support for assymetric overlaps!
  const int ov1 = overlap[0] + 1;
  const double interactionLengthSquare{interactionLength * interactionLength};

  // Due to 2D output, we need to first resize the outer vector to the x dimension
  if constexpr (Mode == C08OffsetMode::c04CellPairs) {
    resultOffsetsC08.resize(ov1);
  }

  // Iteration to build the cell pairs required for the C08 base step, x, y, z reprent the spatial dimension
  for (int x = 0; x <= overlap[0]; ++x) {
    for (int y = 0; y <= overlap[1]; ++y) {
      // First offset of the cell-pair, computed withhout z index
      const int offset1 = threeToOneD(x, y, 0, cellsPerDimIntegral);
      // The distance vector for offset1
      const std::array<double, 3> distVec1{
          {static_cast<double>(x) * cellLength[0], static_cast<double>(y) * cellLength[1], 0.0}};
      for (int z = 0; z <= overlap[2]; ++z) {
        // The z component is always calculated in the same way
        const double distVecZComponent = std::max(0.0, z - 1.0) * cellLength[2];
        // Iteration over the four directions for building the cell pairs with offset1 as first partner
        for (const C08CellDirection &direction : ALL_DIRECTIONS) {
          // Depending on the overlap and the direction from the base cell, we skip certain pairs
          // 1. Exclusion Criteria to skip computations by not including certain cell-combinations
          if (includeCellPair(direction, overlap, x, y, z)) {
            // Calculate the offset2 of the second partaking cell for interactions
            const int offset2 = cellOffsets[calculateCell2Index(ov1, direction, z)];

            // Get the masking values, these depend on the direction and are used to calculate the distance vector for
            // the cell with offset2
            const auto &[maskX, maskY] = toMask(direction);
            const std::array<double, 3> distVec2{
                std::max(0.0, (overlap[0] - x) * maskX + x * (1 - maskX) - 1.0) * cellLength[0],
                std::max(0.0, (overlap[1] - y) * maskY + y * (1 - maskY) - 1.0) * cellLength[1],
                distVecZComponent,
            };

            const auto distSquare = utils::ArrayMath::dot(distVec2, distVec2);
            // 2. Exclusion Criteria to skip computations by not including certain cell-combinations
            //      if the distance between cell centers is actually in interactionLength
            if (distSquare <= interactionLengthSquare) {
              // Calculate the sorting direction if sorting is enabled
              if constexpr (Mode == C08OffsetMode::c08CellPairsSorting) {
                const auto sortDir = computeSortingDirection(cellsPerDimIntegral, distVec1, offset2);
                // Append the cell pairs & their sorting direction
                resultOffsetsC08.emplace_back(offset2, offset1, sortDir);
              } else if constexpr (Mode == C08OffsetMode::c04CellPairs) {
                // The c04CellPairs requires a vector of vector of pairs where the first vector is resolved in the
                // x-axis
                resultOffsetsC08[x].emplace_back(offset2, offset1);
              } else {
                // Just append the pair without anything else
                resultOffsetsC08.emplace_back(offset2, offset1);
              }
            }
          }
        }
      }
    }
  }
  return resultOffsetsC08;
}

/*
 * Explicit Template Instantation - Required since the definition of computePairwiseCellOffsetsC08(..)
 * is not in header file. However, the Mode variable is finite, and all instances can be created before
 * being used by explicit template instantiation
 */

//! @cond Doxygen_Suppress
template std::vector<OffsetPairSorting> computePairwiseCellOffsetsC08<C08OffsetMode::c08CellPairsSorting>(
    const std::array<unsigned long, 3> &cellsPerDimension, const std::array<double, 3> &cellLength,
    double interactionLength);

/** Template Sepcialization to return C08 cell paris without sorting */
template std::vector<OffsetPair> computePairwiseCellOffsetsC08<C08OffsetMode::c08CellPairs>(
    const std::array<unsigned long, 3> &cellsPerDimension, const std::array<double, 3> &cellLength,
    double interactionLength);

/** Template Sepcialization to return C04 cell paris, i.e. C08 resolved on x axis */
template std::vector<OffsetPairVector> computePairwiseCellOffsetsC08<C08OffsetMode::c04CellPairs>(
    const std::array<unsigned long, 3> &cellsPerDimension, const std::array<double, 3> &cellLength,
    double interactionLength);
//! @endcond

}  // namespace autopas::LCC08CellHandlerUtility