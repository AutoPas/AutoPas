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
  if constexpr (Mode == C08OffsetMode::c04NoSorting) {
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
              if constexpr (Mode == C08OffsetMode::sorting) {
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
              } else if constexpr (Mode == C08OffsetMode::c04NoSorting) {
                // The c04NoSorting requires a vector of vector of pairs where the first vector is resolved in x-axis
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

template <C08OffsetMode Mode>
OffsetTripletType<Mode> computeTriwiseCellOffsetsC08(const std::array<unsigned long, 3> &cellsPerDimension,
                                                     const std::array<double, 3> &cellLength,
                                                     double interactionLength) {
  using namespace utils::ArrayMath::literals;
  using namespace internal;
  using utils::ArrayMath::ceilAndCast;

  // Output of the function: Vector of triwise cell indices & projection axis
  OffsetTripletType<Mode> resultOffsetsC08{};

  // Helper function to get minimal distance between two cells
  auto cellDistance = [&](long x1, long y1, long z1, long x2, long y2, long z2) {
    return std::array<double, 3>{std::max(0l, (std::abs(x1 - x2) - 1l)) * cellLength[0],
                                 std::max(0l, (std::abs(y1 - y2) - 1l)) * cellLength[1],
                                 std::max(0l, (std::abs(z1 - z2) - 1l)) * cellLength[2]};
  };

  auto is_valid_distance = [&](long x1, long y1, long z1, long x2, long y2, long z2, auto interactionlengthsq) {
    auto dist = cellDistance(x1, y1, z1, x2, y2, z2);
    return utils::ArrayMath::dot(dist, dist) > interactionlengthsq;
  };

  const std::array<int, 3> overlap{ceilAndCast(interactionLength / cellLength)};

  const auto interactionLengthSquare(interactionLength * interactionLength);

  if constexpr (Mode == C08OffsetMode::sorting) {
    resultOffsetsC08.emplace_back(0, 0, 0, std::array<double, 3>{1., 1., 1.});
  }  // Due to 2D output, we need to first resize the outer vector to the x dimension
  else if constexpr (Mode == C08OffsetMode::c04NoSorting) {
    resultOffsetsC08.resize(overlap[0] + 1);
  } else {
    resultOffsetsC08.emplace_back(0, 0, 0);
  }

  // offsets for the first cell
  for (long x1 = 0; x1 <= static_cast<long>(overlap[0]); ++x1) {
    for (long y1 = 0; y1 <= static_cast<long>(overlap[1]); ++y1) {
      for (long z1 = 0; z1 <= static_cast<long>(overlap[2]); ++z1) {
        // offsets for the second cell
        for (long x2 = 0; x2 <= static_cast<long>(overlap[0]); ++x2) {
          for (long y2 = 0; y2 <= static_cast<long>(overlap[1]); ++y2) {
            for (long z2 = 0; z2 <= static_cast<long>(overlap[2]); ++z2) {
              // check distance between cell 1 and cell 2
              const auto dist12 = cellDistance(x1, y1, z1, x2, y2, z2);
              const double dist12Square = utils::ArrayMath::dot(dist12, dist12);
              if (dist12Square >= interactionLengthSquare) continue;

              // offsets for the third cell
              for (long x3 = 0; x3 <= static_cast<long>(overlap[0]); ++x3) {
                for (long y3 = 0; y3 <= static_cast<long>(overlap[1]); ++y3) {
                  for (long z3 = 0; z3 <= static_cast<long>(overlap[2]); ++z3) {
                    // check distance between cell 1 and cell 3
                    const auto dist13 = cellDistance(x1, y1, z1, x3, y3, z3);
                    const double dist13Square = utils::ArrayMath::dot(dist13, dist13);
                    if (dist13Square >= interactionLengthSquare) continue;

                    // check distance between cell 2 and cell 3
                    const auto dist23 = cellDistance(x2, y2, z2, x3, y3, z3);
                    const double dist23Squared = utils::ArrayMath::dot(dist23, dist23);
                    if (dist23Squared >= interactionLengthSquare) continue;
                    const long offset1 = utils::ThreeDimensionalMapping::threeToOneD(
                        x1, y1, z1, utils::ArrayUtils::static_cast_copy_array<long>(cellsPerDimension));

                    const long offset2 = utils::ThreeDimensionalMapping::threeToOneD(
                        x2, y2, z2, utils::ArrayUtils::static_cast_copy_array<long>(cellsPerDimension));

                    const long offset3 = utils::ThreeDimensionalMapping::threeToOneD(
                        x3, y3, z3, utils::ArrayUtils::static_cast_copy_array<long>(cellsPerDimension));

                    if (offset1 > offset2 or offset2 >= offset3) continue;

                    if ((x1 == 0 or x2 == 0 or x3 == 0) and (y1 == 0 or y2 == 0 or y3 == 0) and
                        (z1 == 0 or z2 == 0 or z3 == 0)) {
                      if constexpr (Mode == C08OffsetMode::sorting) {
                        const std::array<double, 3> sortDirection = {
                            (x1 + x2) * cellLength[0], (y1 + y2) * cellLength[1], (z1 + z2) * cellLength[2]};
                        resultOffsetsC08.emplace_back(offset1, offset2, offset3,
                                                      utils::ArrayMath::normalize(sortDirection));
                      } else if constexpr (Mode == C08OffsetMode::c04NoSorting) {
                        resultOffsetsC08[x1].emplace_back(offset1, offset2, offset3);
                      } else {
                        resultOffsetsC08.emplace_back(offset1, offset2, offset3);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return resultOffsetsC08;
}

template <C08OffsetMode Mode>
OffsetTripletType<Mode> computeTriwiseCellOffsetsC08Optimized(const std::array<unsigned long, 3> &cellsPerDimension,
                                                              const std::array<double, 3> &cellLength,
                                                              double interactionLength) {
  using namespace utils::ArrayMath::literals;
  using namespace internal;
  using utils::ArrayMath::ceilAndCast;

  // Output of the function: Vector of triwise cell indices & projection axis
  OffsetTripletType<Mode> resultOffsetsC08{};

  const std::array<int, 3> overlap{ceilAndCast(interactionLength / cellLength)};

  long ovX = static_cast<long>(overlap[0]) + 1l;
  long ovY = static_cast<long>(overlap[1]) + 1l;
  long ovZ = static_cast<long>(overlap[2]) + 1l;

  const auto interactionLengthSquare(interactionLength * interactionLength);
  const auto cellLength1 = cellLength[0];
  const auto cellLength2 = cellLength[1];
  const auto cellLength3 = cellLength[2];
  std::vector<std::vector<std::vector<long>>> cells(ovX, std::vector<std::vector<long>>(ovY, std::vector<long>(ovZ)));
  auto is_valid_distance = [&](long x1, long y1, long z1, long x2, long y2, long z2) {
    auto dist = std::array<double, 3>{std::max(0l, (std::abs(x1 - x2) - 1l)) * cellLength1,
                                      std::max(0l, (std::abs(y1 - y2) - 1l)) * cellLength2,
                                      std::max(0l, (std::abs(z1 - z2) - 1l)) * cellLength3};
    return utils::ArrayMath::dot(dist, dist) <= interactionLengthSquare;
  };

  auto emplaceOffset = [&](long x1, long y1, long z1, long x2, long y2, long z2, long x3, long y3, long z3) {
    if (is_valid_distance(x2, y2, z2, x3, y3, z3)) {
      if (is_valid_distance(x1, y1, z1, x3, y3, z3)) {
        std::array<double, 3> sortingDirection = {(x1 + x2) * cellLength1, (y1 + y2) * cellLength2,
                                                  (z1 + z2) * cellLength3};
        resultOffsetsC08.emplace_back(cells[x1][y1][z1], cells[x2][y2][z2], cells[x3][y3][z3],
                                      utils::ArrayMath::normalize(sortingDirection));
      }
    }
  };

  for (long x = 0; x < ovX; ++x) {
    for (long y = 0; y < ovY; ++y) {
      for (long z = 0; z < ovZ; ++z) {
        cells[x][y][z] = utils::ThreeDimensionalMapping::threeToOneD(
            x, y, z, utils::ArrayUtils::static_cast_copy_array<long>(cellsPerDimension));
      }
    }
  }
  // base cell offsets
  resultOffsetsC08.emplace_back(0, 0, 0, std::array<double, 3>{1., 1., 1.});
  for (long x = 0; x < ovX; ++x) {
    for (long y = 0; y < ovY; ++y) {
      for (long z = 0; z < ovZ; ++z) {
        if (!is_valid_distance(0, 0, 0, x, y, z)) {
          continue;
        }
        for (long z2 = z + 1; z2 < ovZ; ++z2) {
          emplaceOffset(0, 0, 0, x, y, z, x, y, z2);
        }
        for (long y2 = y + 1; y2 < ovZ; ++y2) {
          for (long z2 = 0; z2 < ovZ; ++z2) {
            emplaceOffset(0, 0, 0, x, y, z, x, y2, z2);
          }
        }
        for (long x2 = x + 1; x2 < ovX; ++x2) {
          for (long y2 = 0; y2 < ovY; ++y2) {
            for (long z2 = 0; z2 < ovZ; ++z2) {
              emplaceOffset(0, 0, 0, x, y, z, x2, y2, z2);
            }
          }
        }
      }
    }
  }
  // 3 planes
  for (long x = 1; x < ovX; ++x) {
    for (long y = 1; y < ovY; ++y) {
      for (long x2 = 1; x2 < ovX; ++x2) {
        for (long z = 1; z < ovZ; ++z) {
          if (!is_valid_distance(x, y, 0, x2, 0, z)) {
            continue;
          }
          for (long y2 = 1; y2 < ovY; ++y2) {
            for (long z2 = 1; z2 < ovZ; ++z2) {
              emplaceOffset(x, y, 0, x2, 0, z, 0, y2, z2);
            }
          }
        }
      }
    }
  }
  // 3 edges
  for (long x = 1; x < ovX; ++x) {
    for (long y = 1; y < ovY; ++y) {
      for (long z = 1; z < ovZ; ++z) {
        emplaceOffset(x, 0, 0, 0, y, 0, 0, 0, z);
      }
    }
  }
  // 2 edges, 1 Wildcard
  // edge x, edge y
  for (long x = 1; x < ovX; ++x) {
    for (long y = 1; y < ovY; ++y) {
      if (!is_valid_distance(x, 0, 0, 0, y, 0)) {
        continue;
      }
      // 2 edges (same cell on one edge counted twice)
      if (cells[x][0][0] < cells[0][y][0]) {
        emplaceOffset(x, 0, 0, x, 0, 0, 0, y, 0);
      } else {
        emplaceOffset(0, y, 0, 0, y, 0, x, 0, 0);
      }
      // 2 edges (2 different cells on one edge)
      for (long x2 = x + 1; x2 < ovX; ++x2) {
        emplaceOffset(x, 0, 0, 0, y, 0, x2, 0, 0);
      }
      // 2 edges, 1 plane/inner
      for (long x2 = 1; x2 < ovX; ++x2) {
        for (long z = 1; z < ovZ; ++z) {
          emplaceOffset(x, 0, 0, 0, y, 0, x2, 0, z);
        }
        for (long y2 = 1; y2 < ovY; ++y2) {
          emplaceOffset(x, 0, 0, 0, y, 0, x2, y2, 0);
          // 2 edges, 1 inner
          for (long z = 1; z < ovZ; ++z) {
            emplaceOffset(x, 0, 0, 0, y, 0, x2, y2, z);
          }
        }
      }
      for (long y2 = 1; y2 < ovY; ++y2) {
        for (long z = 1; z < ovZ; ++z) {
          emplaceOffset(x, 0, 0, 0, y, 0, 0, y2, z);
        }
      }
    }
  }
  // 2 edges (2 different cells on one edge)
  for (long y = 1; y < ovY; ++y) {
    for (long x = 1; x < ovX; ++x) {
      for (long y2 = y + 1; y2 < ovY; ++y2) {
        emplaceOffset(0, y, 0, x, 0, 0, 0, y2, 0);
      }
    }
  }

  // edge x, edge z
  for (long x = 1; x < ovX; ++x) {
    for (long z = 1; z < ovZ; ++z) {
      if (!is_valid_distance(x, 0, 0, 0, 0, z)) {
        continue;
      }
      // 2 edges (same cell on one edge counted twice)
      if (cells[x][0][0] < cells[0][0][z]) {
        emplaceOffset(x, 0, 0, x, 0, 0, 0, 0, z);
      } else {
        emplaceOffset(0, 0, z, 0, 0, z, x, 0, 0);
      }
      // 2 edges (2 different cells on one edge)
      for (long x2 = x + 1; x2 < ovX; ++x2) {
        emplaceOffset(x, 0, 0, 0, 0, z, x2, 0, 0);
      }
      // 2 edges, 1 plane/inner
      for (long x2 = 1; x2 < ovX; ++x2) {
        for (long z2 = 1; z2 < ovZ; ++z2) {
          emplaceOffset(x, 0, 0, 0, 0, z, x2, 0, z2);
        }
        for (long y = 1; y < ovY; ++y) {
          emplaceOffset(x, 0, 0, 0, 0, z, x2, y, 0);
          // 2 edges, 1 inner
          for (long z2 = 1; z2 < ovZ; ++z2) {
            emplaceOffset(x, 0, 0, 0, 0, z, x2, y, z2);
          }
        }
      }
      for (long y = 1; y < ovY; ++y) {
        for (long z2 = 1; z2 < ovZ; ++z2) {
          emplaceOffset(x, 0, 0, 0, 0, z, 0, y, z2);
        }
      }
    }
  }
  // 2 edges (2 different cells on one edge)
  for (long z = 1; z < ovY; ++z) {
    for (long x = 1; x < ovX; ++x) {
      for (long z2 = z + 1; z2 < ovY; ++z2) {
        emplaceOffset(0, 0, z, x, 0, 0, 0, 0, z2);
      }
    }
  }

  // edge y, edge z
  for (long y = 1; y < ovY; ++y) {
    for (long z = 1; z < ovZ; ++z) {
      if (!is_valid_distance(0, y, 0, 0, 0, z)) {
        continue;
      }
      // 2 edges (same cell on one edge counted twice)
      if (cells[0][y][0] < cells[0][0][z]) {
        emplaceOffset(0, y, 0, 0, y, 0, 0, 0, z);
      } else {
        emplaceOffset(0, 0, z, 0, 0, z, 0, y, 0);
      }
      // 2 edges (2 different cells on one edge)
      for (long y2 = y + 1; y2 < ovX; ++y2) {
        emplaceOffset(0, y, 0, 0, 0, z, 0, y2, 0);
      }
      // 2 edges, 1 plane/inner
      for (long x = 1; x < ovX; ++x) {
        for (long z2 = 1; z2 < ovZ; ++z2) {
          emplaceOffset(0, y, 0, 0, 0, z, x, 0, z2);
        }
        for (long y2 = 1; y2 < ovY; ++y2) {
          emplaceOffset(0, y, 0, 0, 0, z, x, y2, 0);
          // 2 edges, 1 inner
          for (long z2 = 1; z2 < ovZ; ++z2) {
            emplaceOffset(0, y, 0, 0, 0, z, x, y2, z2);
          }
        }
      }
      for (long y2 = 1; y2 < ovY; ++y2) {
        for (long z2 = 1; z2 < ovZ; ++z2) {
          emplaceOffset(0, y, 0, 0, 0, z, 0, y2, z2);
        }
      }
    }
  }
  // 2 edges (2 different cells on one edge)
  for (long z = 1; z < ovY; ++z) {
    for (long y = 1; y < ovY; ++y) {
      for (long z2 = z + 1; z2 < ovY; ++z2) {
        emplaceOffset(0, 0, z, 0, y, 0, 0, 0, z2);
      }
    }
  }

  // edge x, plane yz
  for (long x = 1; x < ovX; ++x) {
    for (long y = 1; y < ovY; ++y) {
      for (long z = 1; z < ovZ; ++z) {
        if (!is_valid_distance(x, 0, 0, 0, y, z)) {
          continue;
        }
        // 1 edge, 1 plane (same cell on one edge counted twice)
        if (cells[x][0][0] < cells[0][y][z]) {
          emplaceOffset(x, 0, 0, x, 0, 0, 0, y, z);
        } else {
          emplaceOffset(0, y, z, 0, y, z, x, 0, 0);
        }
        // 1 edge, 1 plane (2 different cells same edge)
        for (long x2 = x + 1; x2 < ovX; ++x2) {
          emplaceOffset(x, 0, 0, x2, 0, 0, 0, y, z);
        }
        // 1 edge, 1 plane, 1 plane/inner
        for (long x2 = 1; x2 < ovX; ++x2) {
          for (long y2 = 1; y2 < ovY; ++y2) {
            emplaceOffset(x, 0, 0, 0, y, z, x2, y2, 0);
            // inner
            for (long z2 = 1; z2 < ovZ; ++z2) {
              emplaceOffset(x, 0, 0, 0, y, z, x2, y2, z2);
            }
          }
        }
        for (long x2 = 1; x2 < ovX; ++x2) {
          for (long z2 = 1; z2 < ovZ; ++z2) {
            emplaceOffset(x, 0, 0, 0, y, z, x2, 0, z2);
          }
        }
        // 1 edge 1 plane (2 different cells same plane)
        for (long y2 = 1; y2 < ovY; ++y2) {
          for (long z2 = 1; z2 < ovZ; ++z2) {
            if (cells[0][y][z] < cells[0][y2][z2]) {
              emplaceOffset(x, 0, 0, 0, y, z, 0, y2, z2);
            }
          }
        }
      }
    }
  }
  // edge y, plane xz
  for (long y = 1; y < ovY; ++y) {
    for (long x = 1; x < ovX; ++x) {
      for (long z = 1; z < ovZ; ++z) {
        if (!is_valid_distance(0, y, 0, x, 0, z)) {
          continue;
        }
        // 1 edge, 1 plane (same cell on one edge counted twice)
        if (cells[0][y][0] < cells[x][0][z]) {
          emplaceOffset(0, y, 0, 0, y, 0, x, 0, z);
        } else {
          emplaceOffset(x, 0, z, x, 0, z, 0, y, 0);
        }
        // 1 edge, 1 plane (2 different cells same edge)
        for (long y2 = y + 1; y2 < ovX; ++y2) {
          emplaceOffset(0, y, 0, 0, y2, 0, x, 0, z);
        }
        // 1 edge, 1 plane, 1 plane/inner
        for (long x2 = 1; x2 < ovX; ++x2) {
          for (long y2 = 1; y2 < ovY; ++y2) {
            emplaceOffset(0, y, 0, x, 0, z, x2, y2, 0);
            // inner
            for (long z2 = 1; z2 < ovZ; ++z2) {
              emplaceOffset(0, y, 0, x, 0, z, x2, y2, z2);
            }
          }
        }
        for (long y2 = 1; y2 < ovY; ++y2) {
          for (long z2 = 1; z2 < ovZ; ++z2) {
            emplaceOffset(0, y, 0, x, 0, z, 0, y2, z2);
          }
        }
        // 1 edge 1 plane (2 different cells same plane)
        for (long x2 = 1; x2 < ovX; ++x2) {
          for (long z2 = 1; z2 < ovZ; ++z2) {
            if (cells[x][0][z] < cells[x2][0][z2]) {
              emplaceOffset(0, y, 0, x, 0, z, x2, 0, z2);
            }
          }
        }
      }
    }
  }
  // edge z, plane xy
  for (long z = 1; z < ovZ; ++z) {
    for (long x = 1; x < ovX; ++x) {
      for (long y = 1; y < ovY; ++y) {
        if (!is_valid_distance(0, 0, z, x, y, 0)) {
          continue;
        }
        // 1 edge, 1 plane (same cell on one edge counted twice)
        if (cells[0][0][z] < cells[x][y][0]) {
          emplaceOffset(0, 0, z, 0, 0, z, x, y, 0);
        } else {
          emplaceOffset(x, y, 0, x, y, 0, 0, 0, z);
        }
        // 1 edge, 1 plane (2 different cells same edge)
        for (long z2 = z + 1; z2 < ovX; ++z2) {
          emplaceOffset(0, 0, z, 0, 0, z2, x, y, 0);
        }
        // 1 edge, 1 plane, 1 plane/inner
        for (long x2 = 1; x2 < ovX; ++x2) {
          for (long y2 = 1; y2 < ovY; ++y2) {
            // 1 edge 1 plane (2 different cells same plane)
            if (cells[x][y][0] < cells[x2][y2][0]) {
              emplaceOffset(0, 0, z, x, y, 0, x2, y2, 0);
            }
            // inner
            for (long z2 = 1; z2 < ovZ; ++z2) {
              emplaceOffset(0, 0, z, x, y, 0, x2, y2, z2);
            }
          }
        }
        for (long y2 = 1; y2 < ovY; ++y2) {
          for (long z2 = 1; z2 < ovZ; ++z2) {
            emplaceOffset(0, 0, z, x, y, 0, 0, y2, z2);
          }
        }
        for (long x2 = 1; x2 < ovX; ++x2) {
          for (long z2 = 1; z2 < ovZ; ++z2) {
            emplaceOffset(0, 0, z, x, y, 0, x2, 0, z2);
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
template std::vector<OffsetPairSorting> computePairwiseCellOffsetsC08<C08OffsetMode::sorting>(
    const std::array<unsigned long, 3> &cellsPerDimension, const std::array<double, 3> &cellLength,
    double interactionLength);

/** Template Sepcialization to return C08 cell paris without sorting */
template std::vector<OffsetPair> computePairwiseCellOffsetsC08<C08OffsetMode::noSorting>(
    const std::array<unsigned long, 3> &cellsPerDimension, const std::array<double, 3> &cellLength,
    double interactionLength);

/** Template Sepcialization to return C04 cell paris, i.e. C08 resolved on x axis */
template std::vector<OffsetPairVector> computePairwiseCellOffsetsC08<C08OffsetMode::c04NoSorting>(
    const std::array<unsigned long, 3> &cellsPerDimension, const std::array<double, 3> &cellLength,
    double interactionLength);

template std::vector<OffsetTripletSorting> computeTriwiseCellOffsetsC08<C08OffsetMode::sorting>(
    const std::array<unsigned long, 3> &cellsPerDimension, const std::array<double, 3> &cellLength,
    double interactionLength);

/** Template Sepcialization to return C08 cell paris without sorting */
template std::vector<OffsetTriplet> computeTriwiseCellOffsetsC08<C08OffsetMode::noSorting>(
    const std::array<unsigned long, 3> &cellsPerDimension, const std::array<double, 3> &cellLength,
    double interactionLength);

/** Template Sepcialization to return C04 cell paris, i.e. C08 resolved on x axis */
template std::vector<OffsetTripletVector> computeTriwiseCellOffsetsC08<C08OffsetMode::c04NoSorting>(
    const std::array<unsigned long, 3> &cellsPerDimension, const std::array<double, 3> &cellLength,
    double interactionLength);
//! @endcond

}  // namespace autopas::LCC08CellHandlerUtility