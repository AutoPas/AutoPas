/**
 * @file DirectSumTest.cpp
 * @author seckler
 * @date 27.04.18
 */

#include "DirectSumTest.h"

TEST_F(DirectSumTest, testParticleAdding) {
  autopas::DirectSum<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> directSum({0., 0., 0.},
                                                                                                {10., 10., 10.}, 1.);
  int id = 1;
  for (double x : {-.5, 0., 5., 9.999, 10., 10.5}) {
    for (double y : {-.5, 0., 5., 9.999, 10., 10.5}) {
      for (double z : {-.5, 0., 5., 9.999, 10., 10.5}) {
        autopas::Particle p({x, y, z}, {0., 0., 0.}, id++);
        if (x == 10. or y == 10. or z == 10. or x == -.5 or y == -.5 or z == -.5 or x == 10.5 or y == 10.5 or
            z == 10.5) {
          EXPECT_ANY_THROW(directSum.addParticle(p));     // outside, therefore not ok!
          EXPECT_NO_THROW(directSum.addHaloParticle(p));  // outside, therefore ok!
        } else {
          EXPECT_NO_THROW(directSum.addParticle(p));       // inside, therefore ok!
          EXPECT_ANY_THROW(directSum.addHaloParticle(p));  // inside, therefore not ok!
        }
      }
    }
  }
}