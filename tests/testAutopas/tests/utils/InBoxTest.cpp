/**
 * @file InBoxTest.cpp
 * @author seckler
 * @date 07.06.18
 */

#include <gtest/gtest.h>

#include "autopas/utils/inBox.h"

namespace InBoxTest {

TEST(InBoxTest, testInBox) {
  std::array<double, 3> lowCorner{1. / 3., 0., 2. / 3.};
  std::array<double, 3> highCorner{2. / 3., 2. / 3., 1.};

  for (double x : {0., 1. / 3., 2. / 3., 1.}) {
    for (double y : {0., 1. / 3., 2. / 3.}) {
      for (double z : {1. / 3., 2. / 3., 1.}) {
        bool inside = x >= 1. / 3. and x < 2. / 3.;
        inside &= y >= 0. and y < 2. / 3.;
        inside &= z >= 2. / 3. and z < 1.;
        const std::array<double, 3> position{x, y, z};
        EXPECT_EQ(inside, autopas::utils::inBox(position, lowCorner, highCorner));
      }
    }
  }
}

TEST(InBoxTest, testNotInBox) {
  std::array<double, 3> lowCorner{1. / 3., 0., 2. / 3.};
  std::array<double, 3> highCorner{2. / 3., 2. / 3., 1.};

  for (double x : {0., 1. / 3., 2. / 3., 1.}) {
    for (double y : {0., 1. / 3., 2. / 3.}) {
      for (double z : {1. / 3., 2. / 3., 1.}) {
        bool inside = x >= 1. / 3. and x < 2. / 3.;
        inside &= y >= 0. and y < 2. / 3.;
        inside &= z >= 2. / 3. and z < 1.;
        const std::array<double, 3> position{x, y, z};
        EXPECT_EQ(not inside, autopas::utils::notInBox(position, lowCorner, highCorner));
      }
    }
  }
}
} // end namespace InBoxTest
