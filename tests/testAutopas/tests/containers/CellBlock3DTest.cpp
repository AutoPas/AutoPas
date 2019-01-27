/**
 * @file CellBlock3DTest.cpp
 * @author tchipev
 * @date 19.01.18
 */

#include "CellBlock3DTest.h"

TEST_F(CellBlock3DTest, test1x1x1) {
  std::array<double, 3> start = {-5., -5., -5.}, dr = {10.0, 10.0, 10.0};
  std::array<int, 3> numParts = {3, 3, 3};
  auto mesh = getMesh(start, dr, numParts);

  int counter = 0;
  for (auto &m : mesh) {
    unsigned long index = _cells_1x1x1.get1DIndexOfPosition(m);
    ASSERT_EQ(index, counter);
    ++counter;
  }
}

TEST_F(CellBlock3DTest, test2x2x2) {
  std::array<double, 3> start = {-2.5, -2.5, -2.5}, dr = {5.0, 5.0, 5.0};
  std::array<int, 3> numParts = {4, 4, 4};
  auto mesh = getMesh(start, dr, numParts);

  int counter = 0;
  for (auto &m : mesh) {
    unsigned long index = _cells_2x2x2.get1DIndexOfPosition(m);
    ASSERT_EQ(index, counter);
    ++counter;
  }
}

TEST_F(CellBlock3DTest, test3x3x3) {
  std::array<double, 3> start = {-1.6, -1.6, -1.6}, dr = {3.3, 3.3, 3.3};
  std::array<int, 3> numParts = {5, 5, 5};
  auto mesh = getMesh(start, dr, numParts);

  int counter = 0;
  for (auto &m : mesh) {
    unsigned long index = _cells_3x3x3.get1DIndexOfPosition(m);
    ASSERT_EQ(index, counter);
    ++counter;
  }
}

void testBoundary(autopas::CellBlock3D<autopas::FullParticleCell<autopas::MoleculeLJ>> &cellBlock,
                  std::array<double, 3> boxMin, std::array<double, 3> boxMax) {
  std::array<std::array<double, 4>, 3> possibleShifts = {};
  for (unsigned short dim = 0; dim < 3; ++dim) {
    possibleShifts[dim] = {std::nextafter(boxMin[dim], -1.),  // slightly below boxmin (outside)
                           boxMin[dim],                       // at boxmin (inside)
                           std::nextafter(boxMax[dim], -1.),  // slightly below boxmax (inside)
                           boxMax[dim]};                      // boxmax (outside)
  }
  std::array<size_t, 3> cellsPerDimWithHalo = cellBlock.getCellsPerDimensionWithHalo();
  std::array<size_t, 3> ind = {};
  for (ind[0] = 0; ind[0] < 4; ++ind[0]) {
    for (ind[1] = 0; ind[1] < 4; ++ind[1]) {
      for (ind[2] = 0; ind[2] < 4; ++ind[2]) {
        std::array<double, 3> position = {possibleShifts[0][ind[0]], possibleShifts[1][ind[1]],
                                          possibleShifts[2][ind[2]]};
        auto pos = cellBlock.get3DIndexOfPosition(position);

        for (int d = 0; d < 3; ++d) {
          switch (ind[d]) {
            case 0:
              // slightly below boxmin (outside)
              EXPECT_EQ(pos[d], 0) << " for d = " << d << ", ind[d] = " << ind[d] << ", position[d] = " << position[d]
                                   << ", cellsPerDimWithHalo[d]: " << cellsPerDimWithHalo[d];
              break;
            case 1:
              // at boxmin (inside)
              EXPECT_EQ(pos[d], 1) << " for d = " << d << ", ind[d] = " << ind[d] << ", position[d] = " << position[d]
                                   << ", cellsPerDimWithHalo[d]: " << cellsPerDimWithHalo[d];
              break;
            case 2:
              // slightly below boxmax (inside)
              EXPECT_EQ(pos[d], cellsPerDimWithHalo[d] - 2)
                  << " for d = " << d << ", ind[d] = " << ind[d] << ", position[d] = " << position[d]
                  << ", cellsPerDimWithHalo[d]: " << cellsPerDimWithHalo[d];
              break;
            case 3:
              // boxmax (outside)
              EXPECT_EQ(pos[d], cellsPerDimWithHalo[d] - 1)
                  << " for d = " << d << ", ind[d] = " << ind[d] << ", position[d] = " << position[d]
                  << ", cellsPerDimWithHalo[d]: " << cellsPerDimWithHalo[d];
              break;
            default: { FAIL(); }
          }
        }
      }
    }
  }
}

TEST_F(CellBlock3DTest, testBoundaries) {
  testBoundary(_cells_1x1x1, {0., 0., 0.}, {10., 10., 10.});

  testBoundary(_cells_2x2x2, {0., 0., 0.}, {10., 10., 10.});

  testBoundary(_cells_3x3x3, {0., 0., 0.}, {10., 10., 10.});

  testBoundary(_cells_11x4x4_nonZeroBoxMin, {2. / 3., 0., 0.}, {1., .125, .125});

  testBoundary(_cells_19x19x19, {0., 0., 0.}, {58.5, 58.5, 58.5});
}

std::vector<std::array<double, 3>> CellBlock3DTest::getMesh(std::array<double, 3> start, std::array<double, 3> dr,
                                                            std::array<int, 3> numParts) const {
  std::vector<std::array<double, 3>> ret;

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
