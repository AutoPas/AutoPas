/**
 * @file LCC08CellHandlerUtility.cpp
 * @author J. Schuhmacher
 * @date 11.07.2024
 */

#include "LCC08CellHandlerUtility.h"

namespace autopas::LCC08CellHandlerUtility {

namespace internal {

constexpr std::pair<int, int> toMaskXY(const C08CellDirection &direction) {
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

std::array<double, 3> computeSortingDirection(const std::array<double, 3> &offset1Vector,
                                              const std::array<double, 3> &offset2Vector,
                                              const std::array<double, 3> &cellLength) {
  using namespace autopas::utils::ArrayMath::literals;
  // In case the sorting direction is 0, 0, 0 ==> fix to 1, 1, 1
  std::array<double, 3> sortDir = offset1Vector - offset2Vector;
  if (std::all_of(sortDir.begin(), sortDir.end(), [](const auto &val) { return val == 0; })) {
    sortDir = {1., 1., 1.};
  }

  // Normalize and return
  return utils::ArrayMath::normalize(sortDir * cellLength);
}

}  // namespace internal

template <C08OffsetMode Mode>
OffsetPairType<Mode> computePairwiseCellOffsetsC08(const std::array<unsigned long, 3> &cellsPerDimension,
                                                   const std::array<double, 3> &cellLength, double interactionLength) {
  using namespace autopas::utils::ArrayMath::literals;
  using namespace internal;
  using utils::ArrayMath::ceilAndCast;
  using utils::ArrayUtils::static_cast_copy_array;
  using utils::ThreeDimensionalMapping::threeToOneD;

  // Output of the function: Vector of pairwise cell indices & projection axis
  OffsetPairType<Mode> resultOffsetsC08{};

  // The overlap with interacting cells (see autopas::CBasedTraversal constructor)
  const std::array<int, 3> overlap{ceilAndCast(interactionLength / cellLength)};
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
      // Calculate the first partaking cell's offset relative to base cell. The first offset never has a z component
      // These vertical interactions (along the z-axis) are all incoperated by the second cell's index
      const int offset1 = threeToOneD(x, y, 0, cellsPerDimIntegral);
      for (int z = 0; z <= overlap[2]; ++z) {
        // The z component is always calculated in the same way and does not depend on the direction
        const double distVecZComponent = std::max(0.0, z - 1.0) * cellLength[2];
        // Iteration over the four directions for building the cell pairs with offset1 as first partner
        // The direction can be interpreted as x_2 and y_2 values to calculate interaction to x_1 (= x) and y_1 (= y)
        for (const C08CellDirection &direction : ALL_DIRECTIONS) {
          // Depending on the overlap and the direction from the base cell, we skip certain pairs
          // 1. Exclusion Criteria to skip computations by not including certain cell-combinations
          if (includeCellPair(direction, overlap, x, y, z)) {
            // Get the masking values for a direction, (e.g. base cell/ front left ==> (0, 0), backLeft ==> (0, 1))
            const auto &[maskX, maskY] = toMaskXY(direction);
            // Calculate the offset2 using the mask. Here, the z-component can have a value (up/down interaction)
            // The x and y values are dependent on the direction (non-zero-values), and the overlap (step-width)
            const int offset2 = threeToOneD(maskX * overlap[0], maskY * overlap[1], z, cellsPerDimIntegral);

            // The "distVec" is the direction between borders of cells.
            // Using masking for x, y: Either max(0, overlap - x - 1) or max(0, x - 1) --> the same holds for y
            const std::array<double, 3> distVec{
                std::max(0.0, (overlap[0] - x) * maskX + x * (1 - maskX) - 1.0) * cellLength[0],
                std::max(0.0, (overlap[1] - y) * maskY + y * (1 - maskY) - 1.0) * cellLength[1],
                distVecZComponent,
            };

            // 2. Exclusion Criteria to skip computations by not including certain cell-combinations
            //      if the distance between cell centers is actually in interactionLength
            if (utils::ArrayMath::dot(distVec, distVec) <= interactionLengthSquare) {
              // Calculate the sorting direction if sorting is enabled
              if constexpr (Mode == C08OffsetMode::c08CellPairsSorting) {
                // These are respectivley the 3D coordinates of the offsets of cell1 and cell2, as double elements
                // The cellLength is utuilized to modfiy the direction of the sortingVector in case the cells
                // are less squarish, but more lengthy
                const auto sortDir = computeSortingDirection(
                    {
                        static_cast<double>(x),
                        static_cast<double>(y),
                        0.0,
                    },
                    {
                        static_cast<double>(maskX * overlap[0]),
                        static_cast<double>(maskY * overlap[1]),
                        static_cast<double>(z),
                    },
                    cellLength);
                // Append the cell pairs & their sorting direction
                resultOffsetsC08.emplace_back(offset2, offset1, sortDir);
              } else if constexpr (Mode == C08OffsetMode::c04CellPairs) {
                // The c04CellPairs requires a vector of vector of pairs where the first vector is resolved in x-axis
                resultOffsetsC08[x].emplace_back(offset2, offset1);
              } else {
                // Like c08CellPairsSorting, but without sorting
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