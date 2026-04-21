/**
 * @file CellBlock3DTest.cpp
 * @author tchipev
 * @date 19.01.18
 */

#include "CellBlock3DTest.h"

#include "autopas/utils/ArrayUtils.h"
#include "autopasTools/generators/GridGenerator.h"

void CellBlock3DTest::SetUp() {
  // Initialize different cell blocks
  _cells_1x1x1 = std::make_unique<autopas::internal::CellBlock3D<FMCell>>(
      _vec1, std::array<double, 3>({0.0, 0.0, 0.0}), std::array<double, 3>({10.0, 10.0, 10.0}), 10.0);
  _cells_1x1x1_cs2 = std::make_unique<autopas::internal::CellBlock3D<FMCell>>(
      _vec1_cs, std::array<double, 3>({0.0, 0.0, 0.0}), std::array<double, 3>({10.0, 10.0, 10.0}), 5.0, 2.0);
  _cells_2x2x2 = std::make_unique<autopas::internal::CellBlock3D<FMCell>>(
      _vec2, std::array<double, 3>({0.0, 0.0, 0.0}), std::array<double, 3>({10.0, 10.0, 10.0}), 5.0);
  _cells_2x2x2_cs05 = std::make_unique<autopas::internal::CellBlock3D<FMCell>>(
      _vec2_cs, std::array<double, 3>({0.0, 0.0, 0.0}), std::array<double, 3>({10.0, 10.0, 10.0}), 10.0, 0.5);
  _cells_3x3x3 = std::make_unique<autopas::internal::CellBlock3D<FMCell>>(
      _vec3, std::array<double, 3>({0.0, 0.0, 0.0}), std::array<double, 3>({10.0, 10.0, 10.0}), 3.0);
  _cells_11x4x4_nonZeroBoxMin = std::make_unique<autopas::internal::CellBlock3D<FMCell>>(
      _vec4, std::array<double, 3>({2. / 3., 0.0, 0.0}), std::array<double, 3>({1., .125, .125}), 3. / 100.);
  _cells_19x19x19 = std::make_unique<autopas::internal::CellBlock3D<FMCell>>(
      _vec19, std::array<double, 3>({0.0, 0.0, 0.0}), std::array<double, 3>({58.5, 58.5, 58.5}), 3.0);
}

void testIndex(std::unique_ptr<autopas::internal::CellBlock3D<FMCell>> &cellBlock, std::array<double, 3> &start,
               std::array<double, 3> &dr, std::array<int, 3> &numParts) {
  auto mesh = CellBlock3DTest::getMesh(start, dr, numParts);

  unsigned long counter = 0ul;
  for (auto &m : mesh) {
    unsigned long index = cellBlock->get1DIndexOfPosition(m);
    ASSERT_EQ(index, counter) << "Pos: [" << autopas::utils::ArrayUtils::to_string(m) << "]";
    ++counter;
  }
}

TEST_F(CellBlock3DTest, test1x1x1) {
  std::array<double, 3> start = {-5., -5., -5.}, dr = {10.0, 10.0, 10.0};
  std::array<int, 3> numParts = {3, 3, 3};
  testIndex(_cells_1x1x1, start, dr, numParts);
}

TEST_F(CellBlock3DTest, test1x1x1_cs2) {
  std::array<double, 3> start = {-5., -5., -5.}, dr = {10.0, 10.0, 10.0};
  std::array<int, 3> numParts = {3, 3, 3};
  testIndex(_cells_1x1x1_cs2, start, dr, numParts);
}

TEST_F(CellBlock3DTest, test2x2x2) {
  std::array<double, 3> start = {-2.5, -2.5, -2.5}, dr = {5.0, 5.0, 5.0};
  std::array<int, 3> numParts = {4, 4, 4};
  testIndex(_cells_2x2x2, start, dr, numParts);
}

TEST_F(CellBlock3DTest, test2x2x2_cs05) {
  std::array<double, 3> start = {-7.5, -7.5, -7.5}, dr = {5.0, 5.0, 5.0};
  std::array<int, 3> numParts = {6, 6, 6};
  testIndex(_cells_2x2x2_cs05, start, dr, numParts);
}

TEST_F(CellBlock3DTest, test3x3x3) {
  std::array<double, 3> start = {-1.6, -1.6, -1.6}, dr = {3.3, 3.3, 3.3};
  std::array<int, 3> numParts = {5, 5, 5};
  testIndex(_cells_3x3x3, start, dr, numParts);
}

void testBoundary(std::unique_ptr<autopas::internal::CellBlock3D<FMCell>> &cellBlock, std::array<double, 3> boxMin,
                  std::array<double, 3> boxMax) {
  std::array<std::array<double, 4>, 3> possibleShifts = {};
  for (unsigned short dim = 0; dim < 3; ++dim) {
    possibleShifts[dim] = {std::nextafter(boxMin[dim], -1.),  // slightly below boxmin (outside)
                           boxMin[dim],                       // at boxmin (inside)
                           std::nextafter(boxMax[dim], -1.),  // slightly below boxmax (inside)
                           boxMax[dim]};                      // boxmax (outside)
  }
  std::array<size_t, 3> cellsPerDimWithHalo = cellBlock->getCellsPerDimensionWithHalo();
  auto haloThickness(cellBlock->getCellsPerInteractionLength());
  std::array<size_t, 3> ind = {};
  for (ind[0] = 0; ind[0] < 4; ++ind[0]) {
    for (ind[1] = 0; ind[1] < 4; ++ind[1]) {
      for (ind[2] = 0; ind[2] < 4; ++ind[2]) {
        std::array<double, 3> position = {possibleShifts[0][ind[0]], possibleShifts[1][ind[1]],
                                          possibleShifts[2][ind[2]]};
        auto pos = cellBlock->get3DIndexOfPosition(position);

        for (int d = 0; d < 3; ++d) {
          switch (ind[d]) {
            case 0:
              // slightly below boxmin (outside)
              EXPECT_EQ(pos[d], haloThickness - 1)
                  << " for d = " << d << ", ind[d] = " << ind[d] << ", position[d] = " << position[d]
                  << ", cellsPerDimWithHalo[d]: " << cellsPerDimWithHalo[d];
              break;
            case 1:
              // at boxmin (inside)
              EXPECT_EQ(pos[d], haloThickness)
                  << " for d = " << d << ", ind[d] = " << ind[d] << ", position[d] = " << position[d]
                  << ", cellsPerDimWithHalo[d]: " << cellsPerDimWithHalo[d];
              break;
            case 2:
              // slightly below boxmax (inside)
              EXPECT_EQ(pos[d], cellsPerDimWithHalo[d] - haloThickness - 1)
                  << " for d = " << d << ", ind[d] = " << ind[d] << ", position[d] = " << position[d]
                  << ", cellsPerDimWithHalo[d]: " << cellsPerDimWithHalo[d];
              break;
            case 3:
              // boxmax (outside)
              EXPECT_EQ(pos[d], cellsPerDimWithHalo[d] - haloThickness)
                  << " for d = " << d << ", ind[d] = " << ind[d] << ", position[d] = " << position[d]
                  << ", cellsPerDimWithHalo[d]: " << cellsPerDimWithHalo[d];
              break;
            default: {
              FAIL();
            }
          }
        }
      }
    }
  }
}

TEST_F(CellBlock3DTest, testBoundaries_cs_geq1) {
  testBoundary(_cells_1x1x1, {0., 0., 0.}, {10., 10., 10.});

  testBoundary(_cells_2x2x2, {0., 0., 0.}, {10., 10., 10.});

  testBoundary(_cells_1x1x1_cs2, {0., 0., 0.}, {10., 10., 10.});

  testBoundary(_cells_3x3x3, {0., 0., 0.}, {10., 10., 10.});

  testBoundary(_cells_11x4x4_nonZeroBoxMin, {2. / 3., 0., 0.}, {1., .125, .125});

  testBoundary(_cells_19x19x19, {0., 0., 0.}, {58.5, 58.5, 58.5});
}

TEST_F(CellBlock3DTest, testBoundaries_cs_leq1) { testBoundary(_cells_2x2x2_cs05, {0., 0., 0.}, {10., 10., 10.}); }

std::vector<std::array<double, 3>> CellBlock3DTest::getMesh(std::array<double, 3> start, std::array<double, 3> dr,
                                                            std::array<int, 3> numParts) {
  std::vector<std::array<double, 3>> ret;
  ret.reserve(numParts[0] * numParts[1] * numParts[2]);

  for (int z = 0; z < numParts[2]; ++z) {
    for (int y = 0; y < numParts[1]; ++y) {
      for (int x = 0; x < numParts[0]; ++x) {
        std::array<double, 3> pos{};
        pos[0] = start[0] + x * dr[0];
        pos[1] = start[1] + y * dr[1];
        pos[2] = start[2] + z * dr[2];
        ret.push_back(pos);
      }
    }
  }
  return ret;
}

size_t getNumberOfParticlesInBox(std::unique_ptr<autopas::internal::CellBlock3D<FMCell>> &cellBlock,
                                 std::vector<FMCell> &vec) {
  const Molecule defaultParticle;
  autopasTools::generators::GridGenerator::fillWithParticles(
      vec, cellBlock->getCellsPerDimensionWithHalo(), cellBlock->getCellsPerDimensionWithHalo(), defaultParticle);
  cellBlock->clearHaloCells();
  return std::accumulate(vec.begin(), vec.end(), 0, [](auto acc, auto &e) { return acc + e.size(); });
}

TEST_F(CellBlock3DTest, testClearHaloParticles) {
  EXPECT_EQ(getNumberOfParticlesInBox(_cells_1x1x1, _vec1), 1);
  EXPECT_EQ(getNumberOfParticlesInBox(_cells_1x1x1_cs2, _vec1_cs), 1);
  EXPECT_EQ(getNumberOfParticlesInBox(_cells_2x2x2, _vec2), 2 * 2 * 2);
  EXPECT_EQ(getNumberOfParticlesInBox(_cells_2x2x2_cs05, _vec2_cs), 2 * 2 * 2);
  EXPECT_EQ(getNumberOfParticlesInBox(_cells_3x3x3, _vec3), 3 * 3 * 3);
  EXPECT_EQ(getNumberOfParticlesInBox(_cells_11x4x4_nonZeroBoxMin, _vec4), 11 * 4 * 4);
  EXPECT_EQ(getNumberOfParticlesInBox(_cells_19x19x19, _vec19), 19 * 19 * 19);
}

/**
 * Checks if the OwnershipState is set correctly, depending on whether it is a cell that can contain only owned
 * particles, only halo particles, owned and halo particles, or if it is an empty cell.
 */
TEST_F(CellBlock3DTest, testCellOwnership) {
  std::size_t numCells = _cells_1x1x1->getNumCells();
  for (int i = 0; i < numCells; i++) {
    if (_cells_1x1x1->cellCanContainHaloParticles(i)) {
      EXPECT_TRUE(toInt64(_cells_1x1x1->getCell(i).getPossibleParticleOwnerships() & autopas::OwnershipState::halo));
    } else if (_cells_1x1x1->cellCanContainOwnedParticles(i)) {
      EXPECT_TRUE(toInt64(_cells_1x1x1->getCell(i).getPossibleParticleOwnerships() & autopas::OwnershipState::owned));
    } else {
      EXPECT_TRUE(_cells_1x1x1->getCell(i).getPossibleParticleOwnerships() == autopas::OwnershipState::dummy);
    }
  }
}

/**
 * Tests if the cell bounding box is calculated correctly by calling getCellBoundingBox with 1D and 3D indices
 * and comparing it to the precalculated values.
 */
TEST_F(CellBlock3DTest, testGetCellBoundingBox) {
  using namespace autopas::utils::ArrayMath::literals;

  constexpr double desiredCellLength = 0.2;
  constexpr double cellSizeFactor = 1.0;
  constexpr size_t cellsPerDimensionWithoutHalo = 3;
  constexpr size_t cellsPerDimensionWithHalo = cellsPerDimensionWithoutHalo + 2;
  constexpr std::array<double, 3> boxMin = {0.1, 0.1, 0.1};
  constexpr std::array<double, 3> boxMax = boxMin + cellsPerDimensionWithoutHalo * desiredCellLength;
  autopas::internal::CellBlock3D<FMCell> cellBlock{_vec1, boxMin, boxMax, desiredCellLength, cellSizeFactor};

  // helper function to avoid code duplication
  // calls getCellBoundingBox with 1D and 3D indices and compares it to the expected values
  const auto testCellBoundingBox = [&](const std::array<size_t, 3> &cellIndex3D,
                                       const std::array<double, 3> &expectedMin,
                                       const std::array<double, 3> &expectedMax) {
    // inner helper function
    const auto actualTest = [&](const auto cellIndex) {
      using autopas::utils::ArrayUtils::operator<<;
      auto [actualMin, actualMax] = cellBlock.getCellBoundingBox(cellIndex);
      EXPECT_THAT(actualMin, testing::Pointwise(testing::DoubleNear(1e-15), expectedMin));
      EXPECT_THAT(actualMax, testing::Pointwise(testing::DoubleNear(1e-15), expectedMax));
      if (HasFailure()) {
        std::cout << "Failed for cellIndex: " << cellIndex << std::endl;
      }
    };
    const auto cellIndex1D = autopas::utils::ThreeDimensionalMapping::threeToOneD(
        cellIndex3D, {cellsPerDimensionWithHalo, cellsPerDimensionWithHalo, cellsPerDimensionWithHalo});
    actualTest(cellIndex1D);
    actualTest(cellIndex3D);
  };

  // first cell in the halo box
  testCellBoundingBox({0, 0, 0}, boxMin - desiredCellLength, boxMin);
  // first cell in the actual box
  testCellBoundingBox({1, 1, 1}, boxMin, boxMin + desiredCellLength);
  // last cell in the actual box
  testCellBoundingBox({cellsPerDimensionWithHalo - 2, cellsPerDimensionWithHalo - 2, cellsPerDimensionWithHalo - 2},
                      boxMax - desiredCellLength, boxMax);
  // last cell in the halo box
  testCellBoundingBox({cellsPerDimensionWithHalo - 1, cellsPerDimensionWithHalo - 1, cellsPerDimensionWithHalo - 1},
                      boxMax, boxMax + desiredCellLength);
}