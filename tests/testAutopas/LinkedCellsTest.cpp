/**
 * @file LinkedCellsTest.cpp
 * @author seckler
 * @date 27.04.18
 */

#include "LinkedCellsTest.h"
#include <cells/FullParticleCell.h>
#include <containers/LinkedCells.h>
#include <particles/Particle.h>

TEST_F(LinkedCellsTest, testParticleAdding) {
  autopas::LinkedCells<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> linkedCells(
      {0., 0., 0.}, {10., 10., 10.}, 1.);
  int id = 1;
  for (double x : {-1.5, -.5, 0., 5., 9.999, 10., 10.5, 11.5}) {
    for (double y : {-1.5, -.5, 0., 5., 9.999, 10., 10.5, 11.5}) {
      for (double z : {-1.5, -.5, 0., 5., 9.999, 10., 10.5, 11.5}) {
        autopas::Particle p({x, y, z}, {0., 0., 0.}, id++);
        if (x == -1.5 or y == -1.5 or z == -1.5 or x == 11.5 or y == 11.5 or z == 11.5) {
          EXPECT_ANY_THROW(linkedCells.addParticle(p));      // outside, therefore not ok!
          EXPECT_ANY_THROW(linkedCells.addHaloParticle(p));  // much outside, therefore not ok!
        } else if (x == 10. or y == 10. or z == 10. or x == -.5 or y == -.5 or z == -.5 or x == 10.5 or y == 10.5 or
                   z == 10.5) {
          EXPECT_ANY_THROW(linkedCells.addParticle(p));     // outside, therefore not ok!
          EXPECT_NO_THROW(linkedCells.addHaloParticle(p));  // outside, therefore ok!
        } else {
          EXPECT_NO_THROW(linkedCells.addParticle(p));       // inside, therefore ok!
          EXPECT_ANY_THROW(linkedCells.addHaloParticle(p));  // inside, therefore not ok!
        }
      }
    }
  }
}

TEST_F(LinkedCellsTest, testCheckUpdateContainerNeededNoMove) {
  {
    autopas::LinkedCells<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> linkedCells(
        {0., 0., 0.}, {10., 10., 10.}, 1.);
    int id = 1;
    for (double x : {-.5, 0., 5., 9.999, 10., 10.5}) {
      for (double y : {-.5, 0., 5., 9.999, 10., 10.5}) {
        for (double z : {-.5, 0., 5., 9.999, 10., 10.5}) {
          autopas::Particle p({x, y, z}, {0., 0., 0.}, id++);
          bool halo = false;
          for (int d = 0; d < 3; d++) {
            if (p.getR()[d] < 0. or p.getR()[d] >= 10.) {
              halo = true;
            }
          }
          if (halo) {
            linkedCells.addHaloParticle(p);
          } else {
            linkedCells.addParticle(p);
          }
          EXPECT_FALSE(linkedCells.isContainerUpdateNeeded());
        }
      }
    }
  }
  {
    autopas::LinkedCells<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> linkedCells(
        {0., 0., 0.}, {10., 10., 10.}, 3.);
    int id = 1;
    for (double x : {-1.5, -.5, 0., 1. / 3, 2. / 3, 10., 10.5, 11.5}) {
      for (double y : {-1.5, -.5, 0., 1. / 3, 2. / 3, 10., 10.5, 11.5}) {
        for (double z : {-1.5, -.5, 0., 1. / 3, 2. / 3, 10., 10.5, 11.5}) {
          autopas::Particle p({x, y, z}, {0., 0., 0.}, id++);
          bool halo = false;
          for (int d = 0; d < 3; d++) {
            if (p.getR()[d] < 0. or p.getR()[d] >= 10.) {
              halo = true;
            }
          }
          if (halo) {
            linkedCells.addHaloParticle(p);
          } else {
            linkedCells.addParticle(p);
          }
          EXPECT_FALSE(linkedCells.isContainerUpdateNeeded());
        }
      }
    }
  }
}