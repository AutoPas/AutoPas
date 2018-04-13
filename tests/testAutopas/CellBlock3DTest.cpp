/*
 * CellBlock3DTest.cpp
 *
 *  Created on: 19 Jan 2018
 *      Author: tchipevn
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

TEST_F(CellBlock3DTest, testBoundaries) {
  std::array<double, 3> position{0, 0, 0};
  for (int x_ind : {0, 1}) {
    for (int y_ind : {0, 1}) {
      for (int z_ind : {0, 1}) {
        position = {x_ind * 10., y_ind * 10., z_ind * 10.};

        auto pos = _cells_1x1x1.get3DIndexOfPosition(position);

        ASSERT_EQ(pos[0], 1 * x_ind + 1);
        ASSERT_EQ(pos[1], 1 * y_ind + 1);
        ASSERT_EQ(pos[2], 1 * z_ind + 1);

        pos = _cells_2x2x2.get3DIndexOfPosition(position);

        ASSERT_EQ(pos[0], 2 * x_ind + 1);
        ASSERT_EQ(pos[1], 2 * y_ind + 1);
        ASSERT_EQ(pos[2], 2 * z_ind + 1);

        pos = _cells_3x3x3.get3DIndexOfPosition(position);

        ASSERT_EQ(pos[0], 3 * x_ind + 1);
        ASSERT_EQ(pos[1], 3 * y_ind + 1);
        ASSERT_EQ(pos[2], 3 * z_ind + 1);

        position = {x_ind ? 1. : 2. / 3, y_ind * .125, z_ind * .125};

        pos = _cells_11x4x4_nonZeroBoxMin.get3DIndexOfPosition(position);

        ASSERT_EQ(pos[0], 11 * x_ind + 1);
        ASSERT_EQ(pos[1], 4 * y_ind + 1);
        ASSERT_EQ(pos[2], 4 * z_ind + 1);

      }
    }
  }
}

std::vector<std::array<double, 3>> CellBlock3DTest::getMesh(
    std::array<double, 3> start, std::array<double, 3> dr,
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
