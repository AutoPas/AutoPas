/**
 * @file CellBlock3DTest.cpp
 * @author tchipev
 * @date 19.01.18
 */

#include "CellBlock3DTest.h"

#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/generators/GridGenerator.h"

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
  autopas::generators::GridGenerator::fillWithParticles(vec, cellBlock->getCellsPerDimensionWithHalo(),
                                                        cellBlock->getCellsPerDimensionWithHalo(), defaultParticle);
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

/**
 * Tests whether updateHaloParticle() finds and updates matching halo particles in containing and neighboring halo
 * cells (including follow-up updates after a particle has already left its cell, at both lower and upper boundaries)
 * while ignoring owned cells, mismatched IDs, and particles outside allowedDistance.
 */
TEST_F(CellBlock3DTest, testUpdateHaloParticle) {
  constexpr double cellLength = 1.0;
  constexpr std::array<double, 3> boxMin = {0.0, 0.0, 0.0};
  constexpr std::array<double, 3> boxMax = {3.0, 3.0, 3.0};
  std::vector<FMCell> cells;
  autopas::internal::CellBlock3D<FMCell> cellBlock{cells, boxMin, boxMax, cellLength, 1.0};

  // Case 1 (lower boundary): Matching halo particle in its containing halo cell {0, 2, 2} is updated.
  {
    Molecule original({-0.5, 1.5, 1.5}, {0.0, 0.0, 0.0}, 1);
    original.setOwnershipState(autopas::OwnershipState::halo);
    cellBlock.getContainingCell(original.getR()).addParticle(original);

    Molecule updated({-0.45, 1.5, 1.5}, {1.0, 2.0, 3.0}, 1);
    updated.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_TRUE(cellBlock.updateHaloParticle(updated, 0.1));

    const auto &stored = cellBlock.getCell({0, 2, 2})[0];
    EXPECT_EQ(stored.getR(), updated.getR());
    EXPECT_EQ(stored.getV(), updated.getV());
    cellBlock.clearHaloCells();
  }

  // Case 2 (lower boundary): Halo particle moved across a cell boundary into an adjacent halo cell within
  // allowedDistance, and is updated a second time while already outside its stored cell.
  {
    // Original stored in halo cell {0, 0, 1}, updated positions are in adjacent halo cell {0, 1, 1}.
    Molecule original({-0.05, -0.02, 0.05}, {0.0, 0.0, 0.0}, 2);
    original.setOwnershipState(autopas::OwnershipState::halo);
    cellBlock.getContainingCell(original.getR()).addParticle(original);

    Molecule updated({-0.05, 0.05, 0.05}, {4.0, 5.0, 6.0}, 2);
    updated.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_TRUE(cellBlock.updateHaloParticle(updated, 0.1));

    const auto &storedAfterFirst = cellBlock.getCell({0, 0, 1})[0];
    EXPECT_EQ(storedAfterFirst.getR(), updated.getR());
    EXPECT_EQ(storedAfterFirst.getV(), updated.getV());

    Molecule updatedAgain({-0.05, 0.07, 0.05}, {7.0, 8.0, 9.0}, 2);
    updatedAgain.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_TRUE(cellBlock.updateHaloParticle(updatedAgain, 0.1));

    const auto &storedAfterSecond = cellBlock.getCell({0, 0, 1})[0];
    EXPECT_EQ(storedAfterSecond.getR(), updatedAgain.getR());
    EXPECT_EQ(storedAfterSecond.getV(), updatedAgain.getV());
    cellBlock.clearHaloCells();
  }

  // Case 3 (lower boundary): Halo particle moved slightly into an owned cell -> still finds and updates it in the
  // neighboring halo cell, including on a subsequent update while already inside the owned cell.
  {
    Molecule original({-0.02, 0.5, 0.5}, {0.0, 0.0, 0.0}, 3);
    original.setOwnershipState(autopas::OwnershipState::halo);
    cellBlock.getContainingCell(original.getR()).addParticle(original);

    Molecule updated({0.02, 0.5, 0.5}, {2.0, 3.0, 4.0}, 3);
    updated.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_TRUE(cellBlock.updateHaloParticle(updated, 0.1));

    const auto &storedAfterFirst = cellBlock.getCell({0, 1, 1})[0];
    EXPECT_EQ(storedAfterFirst.getR(), updated.getR());
    EXPECT_EQ(storedAfterFirst.getV(), updated.getV());

    Molecule updatedAgain({0.06, 0.5, 0.5}, {5.0, 6.0, 7.0}, 3);
    updatedAgain.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_TRUE(cellBlock.updateHaloParticle(updatedAgain, 0.1));

    const auto &storedAfterSecond = cellBlock.getCell({0, 1, 1})[0];
    EXPECT_EQ(storedAfterSecond.getR(), updatedAgain.getR());
    EXPECT_EQ(storedAfterSecond.getV(), updatedAgain.getV());
    cellBlock.clearHaloCells();
  }

  // Case 4: Matching particle in a nearby owned cell must NOT be updated.
  {
    Molecule owned({0.02, 0.02, 0.02}, {0.0, 0.0, 0.0}, 4);
    owned.setOwnershipState(autopas::OwnershipState::owned);
    cellBlock.getContainingCell(owned.getR()).addParticle(owned);

    Molecule haloUpdate({-0.02, 0.02, 0.02}, {9.0, 9.0, 9.0}, 4);
    haloUpdate.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_FALSE(cellBlock.updateHaloParticle(haloUpdate, 0.1));

    const auto &stored = cellBlock.getCell({1, 1, 1})[0];
    EXPECT_EQ(stored.getR(), owned.getR());
    EXPECT_EQ(stored.getV(), owned.getV());
    cellBlock.getCell({1, 1, 1}).clear();
  }

  // Case 5: Particle outside allowedDistance or with a mismatched ID is not updated.
  {
    Molecule original({-0.5, 1.5, 1.5}, {0.0, 0.0, 0.0}, 5);
    original.setOwnershipState(autopas::OwnershipState::halo);
    cellBlock.getContainingCell(original.getR()).addParticle(original);

    Molecule wrongId({-0.48, 1.5, 1.5}, {1.0, 1.0, 1.0}, 99);
    wrongId.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_FALSE(cellBlock.updateHaloParticle(wrongId, 0.1));

    Molecule tooFar({-0.3, 1.5, 1.5}, {1.0, 1.0, 1.0}, 5);
    tooFar.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_FALSE(cellBlock.updateHaloParticle(tooFar, 0.1));

    const auto &stored = cellBlock.getCell({0, 2, 2})[0];
    EXPECT_EQ(stored.getR(), original.getR());
    EXPECT_EQ(stored.getV(), original.getV());
    cellBlock.clearHaloCells();
  }

  // Case 6: Position deep inside owned domain does not update owned particles and returns false.
  {
    Molecule owned({1.5, 1.5, 1.5}, {0.0, 0.0, 0.0}, 6);
    owned.setOwnershipState(autopas::OwnershipState::owned);
    cellBlock.getContainingCell(owned.getR()).addParticle(owned);

    Molecule haloUpdate({1.51, 1.5, 1.5}, {7.0, 8.0, 9.0}, 6);
    haloUpdate.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_FALSE(cellBlock.updateHaloParticle(haloUpdate, 0.1));

    const auto &stored = cellBlock.getCell({2, 2, 2})[0];
    EXPECT_EQ(stored.getR(), owned.getR());
    EXPECT_EQ(stored.getV(), owned.getV());
    cellBlock.getCell({2, 2, 2}).clear();
  }

  // Case 7 (upper boundary, Case 1): Matching halo particle in its containing upper halo cell {4, 2, 2} is updated.
  {
    Molecule original({3.5, 1.5, 1.5}, {0.0, 0.0, 0.0}, 7);
    original.setOwnershipState(autopas::OwnershipState::halo);
    cellBlock.getContainingCell(original.getR()).addParticle(original);

    Molecule updated({3.45, 1.5, 1.5}, {1.0, 2.0, 3.0}, 7);
    updated.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_TRUE(cellBlock.updateHaloParticle(updated, 0.1));

    const auto &stored = cellBlock.getCell({4, 2, 2})[0];
    EXPECT_EQ(stored.getR(), updated.getR());
    EXPECT_EQ(stored.getV(), updated.getV());
    cellBlock.clearHaloCells();
  }

  // Case 8 (upper boundary, Case 2): Halo particle moved into an adjacent upper halo cell and updated a second time.
  {
    // Original stored in upper halo cell {4, 4, 3}, updated positions are in adjacent halo cell {4, 3, 3}.
    Molecule original({3.05, 3.02, 2.95}, {0.0, 0.0, 0.0}, 8);
    original.setOwnershipState(autopas::OwnershipState::halo);
    cellBlock.getContainingCell(original.getR()).addParticle(original);

    Molecule updated({3.05, 2.95, 2.95}, {4.0, 5.0, 6.0}, 8);
    updated.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_TRUE(cellBlock.updateHaloParticle(updated, 0.1));

    const auto &storedAfterFirst = cellBlock.getCell({4, 4, 3})[0];
    EXPECT_EQ(storedAfterFirst.getR(), updated.getR());
    EXPECT_EQ(storedAfterFirst.getV(), updated.getV());

    Molecule updatedAgain({3.05, 2.93, 2.95}, {7.0, 8.0, 9.0}, 8);
    updatedAgain.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_TRUE(cellBlock.updateHaloParticle(updatedAgain, 0.1));

    const auto &storedAfterSecond = cellBlock.getCell({4, 4, 3})[0];
    EXPECT_EQ(storedAfterSecond.getR(), updatedAgain.getR());
    EXPECT_EQ(storedAfterSecond.getV(), updatedAgain.getV());
    cellBlock.clearHaloCells();
  }

  // Case 9 (upper boundary, Case 3): Upper halo particle moved into an owned cell and updated a second time.
  {
    // Original stored in upper halo cell {4, 3, 3}, updated positions are in owned cell {3, 3, 3}.
    Molecule original({3.02, 2.5, 2.5}, {0.0, 0.0, 0.0}, 9);
    original.setOwnershipState(autopas::OwnershipState::halo);
    cellBlock.getContainingCell(original.getR()).addParticle(original);

    Molecule updated({2.98, 2.5, 2.5}, {2.0, 3.0, 4.0}, 9);
    updated.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_TRUE(cellBlock.updateHaloParticle(updated, 0.1));

    const auto &storedAfterFirst = cellBlock.getCell({4, 3, 3})[0];
    EXPECT_EQ(storedAfterFirst.getR(), updated.getR());
    EXPECT_EQ(storedAfterFirst.getV(), updated.getV());

    Molecule updatedAgain({2.94, 2.5, 2.5}, {5.0, 6.0, 7.0}, 9);
    updatedAgain.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_TRUE(cellBlock.updateHaloParticle(updatedAgain, 0.1));

    const auto &storedAfterSecond = cellBlock.getCell({4, 3, 3})[0];
    EXPECT_EQ(storedAfterSecond.getR(), updatedAgain.getR());
    EXPECT_EQ(storedAfterSecond.getV(), updatedAgain.getV());
    cellBlock.clearHaloCells();
  }
}

/**
 * Tests updateHaloParticle() with cellSizeFactor = 0.5 (two halo cell layers per side) for same-cell,
 * adjacent-halo-cell, and halo-into-owned-cell updates (with follow-up updates) at both lower and upper boundaries.
 */
TEST_F(CellBlock3DTest, testUpdateHaloParticleCSF05) {
  constexpr double interactionLength = 1.0;
  constexpr double cellSizeFactor = 0.5;
  constexpr std::array<double, 3> boxMin = {0.0, 0.0, 0.0};
  constexpr std::array<double, 3> boxMax = {3.0, 3.0, 3.0};
  std::vector<FMCell> cells;
  autopas::internal::CellBlock3D<FMCell> cellBlock{cells, boxMin, boxMax, interactionLength, cellSizeFactor};
  ASSERT_EQ(cellBlock.getCellsPerInteractionLength(), 2ul);

  // Lower boundary, Case 1: Matching halo particle in outer halo cell {0, 4, 4} is updated in place.
  {
    Molecule original({-0.75, 1.25, 1.25}, {0.0, 0.0, 0.0}, 1);
    original.setOwnershipState(autopas::OwnershipState::halo);
    cellBlock.getContainingCell(original.getR()).addParticle(original);

    Molecule updated({-0.70, 1.25, 1.25}, {1.0, 2.0, 3.0}, 1);
    updated.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_TRUE(cellBlock.updateHaloParticle(updated, 0.1));

    const auto &stored = cellBlock.getCell({0, 4, 4})[0];
    EXPECT_EQ(stored.getR(), updated.getR());
    EXPECT_EQ(stored.getV(), updated.getV());
    cellBlock.clearHaloCells();
  }

  // Lower boundary, Case 2: Halo particle crosses from outer halo layer {0, 4, 4} into inner halo layer {1, 4, 4}
  // and is updated again while already outside its stored cell.
  {
    Molecule original({-0.52, 1.25, 1.25}, {0.0, 0.0, 0.0}, 2);
    original.setOwnershipState(autopas::OwnershipState::halo);
    cellBlock.getContainingCell(original.getR()).addParticle(original);

    Molecule updated({-0.48, 1.25, 1.25}, {4.0, 5.0, 6.0}, 2);
    updated.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_TRUE(cellBlock.updateHaloParticle(updated, 0.1));

    const auto &storedAfterFirst = cellBlock.getCell({0, 4, 4})[0];
    EXPECT_EQ(storedAfterFirst.getR(), updated.getR());
    EXPECT_EQ(storedAfterFirst.getV(), updated.getV());

    Molecule updatedAgain({-0.44, 1.25, 1.25}, {7.0, 8.0, 9.0}, 2);
    updatedAgain.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_TRUE(cellBlock.updateHaloParticle(updatedAgain, 0.1));

    const auto &storedAfterSecond = cellBlock.getCell({0, 4, 4})[0];
    EXPECT_EQ(storedAfterSecond.getR(), updatedAgain.getR());
    EXPECT_EQ(storedAfterSecond.getV(), updatedAgain.getV());
    cellBlock.clearHaloCells();
  }

  // Lower boundary, Case 3: Halo particle moves from inner halo layer {1, 2, 2} into owned cell {2, 2, 2}
  // and is updated again while already inside the owned cell.
  {
    Molecule original({-0.02, 0.25, 0.25}, {0.0, 0.0, 0.0}, 3);
    original.setOwnershipState(autopas::OwnershipState::halo);
    cellBlock.getContainingCell(original.getR()).addParticle(original);

    Molecule updated({0.02, 0.25, 0.25}, {2.0, 3.0, 4.0}, 3);
    updated.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_TRUE(cellBlock.updateHaloParticle(updated, 0.1));

    const auto &storedAfterFirst = cellBlock.getCell({1, 2, 2})[0];
    EXPECT_EQ(storedAfterFirst.getR(), updated.getR());
    EXPECT_EQ(storedAfterFirst.getV(), updated.getV());

    Molecule updatedAgain({0.06, 0.25, 0.25}, {5.0, 6.0, 7.0}, 3);
    updatedAgain.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_TRUE(cellBlock.updateHaloParticle(updatedAgain, 0.1));

    const auto &storedAfterSecond = cellBlock.getCell({1, 2, 2})[0];
    EXPECT_EQ(storedAfterSecond.getR(), updatedAgain.getR());
    EXPECT_EQ(storedAfterSecond.getV(), updatedAgain.getV());
    cellBlock.clearHaloCells();
  }

  // Upper boundary, Case 1: Matching halo particle in outer upper halo cell {9, 4, 4} is updated in place.
  {
    Molecule original({3.75, 1.25, 1.25}, {0.0, 0.0, 0.0}, 4);
    original.setOwnershipState(autopas::OwnershipState::halo);
    cellBlock.getContainingCell(original.getR()).addParticle(original);

    Molecule updated({3.70, 1.25, 1.25}, {1.0, 2.0, 3.0}, 4);
    updated.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_TRUE(cellBlock.updateHaloParticle(updated, 0.1));

    const auto &stored = cellBlock.getCell({9, 4, 4})[0];
    EXPECT_EQ(stored.getR(), updated.getR());
    EXPECT_EQ(stored.getV(), updated.getV());
    cellBlock.clearHaloCells();
  }

  // Upper boundary, Case 2: Halo particle crosses from outer upper halo layer {9, 4, 4} into inner upper halo layer
  // {8, 4, 4} and is updated again while already outside its stored cell.
  {
    Molecule original({3.52, 1.25, 1.25}, {0.0, 0.0, 0.0}, 5);
    original.setOwnershipState(autopas::OwnershipState::halo);
    cellBlock.getContainingCell(original.getR()).addParticle(original);

    Molecule updated({3.48, 1.25, 1.25}, {4.0, 5.0, 6.0}, 5);
    updated.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_TRUE(cellBlock.updateHaloParticle(updated, 0.1));

    const auto &storedAfterFirst = cellBlock.getCell({9, 4, 4})[0];
    EXPECT_EQ(storedAfterFirst.getR(), updated.getR());
    EXPECT_EQ(storedAfterFirst.getV(), updated.getV());

    Molecule updatedAgain({3.44, 1.25, 1.25}, {7.0, 8.0, 9.0}, 5);
    updatedAgain.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_TRUE(cellBlock.updateHaloParticle(updatedAgain, 0.1));

    const auto &storedAfterSecond = cellBlock.getCell({9, 4, 4})[0];
    EXPECT_EQ(storedAfterSecond.getR(), updatedAgain.getR());
    EXPECT_EQ(storedAfterSecond.getV(), updatedAgain.getV());
    cellBlock.clearHaloCells();
  }

  // Upper boundary, Case 3: Halo particle moves from inner upper halo layer {8, 7, 7} into owned cell {7, 7, 7}
  // and is updated again while already inside the owned cell.
  {
    Molecule original({3.02, 2.75, 2.75}, {0.0, 0.0, 0.0}, 6);
    original.setOwnershipState(autopas::OwnershipState::halo);
    cellBlock.getContainingCell(original.getR()).addParticle(original);

    Molecule updated({2.98, 2.75, 2.75}, {2.0, 3.0, 4.0}, 6);
    updated.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_TRUE(cellBlock.updateHaloParticle(updated, 0.1));

    const auto &storedAfterFirst = cellBlock.getCell({8, 7, 7})[0];
    EXPECT_EQ(storedAfterFirst.getR(), updated.getR());
    EXPECT_EQ(storedAfterFirst.getV(), updated.getV());

    Molecule updatedAgain({2.94, 2.75, 2.75}, {5.0, 6.0, 7.0}, 6);
    updatedAgain.setOwnershipState(autopas::OwnershipState::halo);
    EXPECT_TRUE(cellBlock.updateHaloParticle(updatedAgain, 0.1));

    const auto &storedAfterSecond = cellBlock.getCell({8, 7, 7})[0];
    EXPECT_EQ(storedAfterSecond.getR(), updatedAgain.getR());
    EXPECT_EQ(storedAfterSecond.getV(), updatedAgain.getV());
    cellBlock.clearHaloCells();
  }
}
