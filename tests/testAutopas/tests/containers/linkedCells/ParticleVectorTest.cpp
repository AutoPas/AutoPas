/**
 * @file ParticleVectorTest.cpp
 *
 * @author lunaticcoding
 * @date 08.06.2020
 */

#include "ParticleVectorTest.h"

#include <gmock/gmock-generated-matchers.h>

#include "autopas/containers/linkedCells/ParticleVector.h"

ParticleVectorTest::ParticleVectorTest() = default;

TEST_F(ParticleVectorTest, testdirtySizeAfterMarkAsClean) {
  auto particleVector = ParticleVector<autopas::Particle>();

  autopas::Particle p1({0.5, 0.5, 0.5}, {0, 0, 0}, 0);
  autopas::Particle p2({1.5, 1.5, 1.5}, {0, 0, 0}, 1);
  autopas::Particle p3({1.6, 1.5, 1.5}, {0, 0, 0}, 2);
  autopas::Particle p4({2.5, 1.5, 1.5}, {0, 0, 0}, 3);
  autopas::Particle p5({2.5, 2.5, 2.5}, {0, 0, 0}, 4);
  autopas::Particle p6({2.5, 1.5, 1.5}, {0, 0, 0}, 3);
  autopas::Particle p7({2.5, 2.5, 2.5}, {0, 0, 0}, 4);

  particleVector.push_back(p1);
  particleVector.push_back(p2);
  particleVector.push_back(p3);
  particleVector.push_back(p4);
  particleVector.push_back(p5);

  particleVector.markAsClean();

  particleVector.push_back(p6);
  particleVector.push_back(p7);

  EXPECT_EQ(particleVector.totalSize(), 7);
  EXPECT_EQ(particleVector.dirtySize(), 2);
}

TEST_F(ParticleVectorTest, testDirtySizeAfterImplicitResize) {
  auto particleVector = ParticleVector<autopas::Particle>();

  autopas::Particle p1({0.5, 0.5, 0.5}, {0, 0, 0}, 0);
  autopas::Particle p2({1.5, 1.5, 1.5}, {0, 0, 0}, 1);
  autopas::Particle p3({1.6, 1.5, 1.5}, {0, 0, 0}, 2);
  autopas::Particle p4({2.5, 1.5, 1.5}, {0, 0, 0}, 3);
  autopas::Particle p5({2.5, 2.5, 2.5}, {0, 0, 0}, 4);
  autopas::Particle p6({2.5, 2.5, 2.5}, {0, 0, 0}, 4);
  autopas::Particle p7({2.5, 2.5, 2.5}, {0, 0, 0}, 4);
  autopas::Particle p8({2.5, 2.5, 2.5}, {0, 0, 0}, 4);
  autopas::Particle p9({2.5, 2.5, 2.5}, {0, 0, 0}, 4);

  particleVector.push_back(p1);
  particleVector.push_back(p2);
  particleVector.push_back(p3);

  particleVector.markAsClean();

  particleVector.push_back(p4);
  particleVector.push_back(p5);
  particleVector.push_back(p6);
  particleVector.push_back(p7);
  particleVector.push_back(p8);
  particleVector.push_back(p9);

  EXPECT_EQ(particleVector.totalSize(), 9);
  EXPECT_EQ(particleVector.dirtySize(), 9);
}
