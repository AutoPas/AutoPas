/**
 * @file ParticleVectorTest.cpp
 *
 * @author lunaticcoding
 * @date 08.06.2020
 */

#include "ParticleVectorTest.h"

#include "autopas/containers/linkedCells/ParticleVector.h"

ParticleVectorTest::ParticleVectorTest() = default;

TEST_F(ParticleVectorTest, testdirtySizeAfterMarkAsClean) {
  auto particleVector = ParticleVector<autopas::ParticleBaseFP64>();

  autopas::ParticleBaseFP64 p1({0.5, 0.5, 0.5}, {0, 0, 0}, 0);
  autopas::ParticleBaseFP64 p2({1.5, 1.5, 1.5}, {0, 0, 0}, 1);
  autopas::ParticleBaseFP64 p3({1.6, 1.5, 1.5}, {0, 0, 0}, 2);
  autopas::ParticleBaseFP64 p4({2.5, 1.5, 1.5}, {0, 0, 0}, 3);
  autopas::ParticleBaseFP64 p5({2.5, 2.5, 2.5}, {0, 0, 0}, 4);
  autopas::ParticleBaseFP64 p6({2.5, 1.5, 1.5}, {0, 0, 0}, 3);
  autopas::ParticleBaseFP64 p7({2.5, 2.5, 2.5}, {0, 0, 0}, 4);

  particleVector.push_back(p1);
  particleVector.push_back(p2);
  particleVector.push_back(p3);
  particleVector.push_back(p4);
  particleVector.push_back(p5);

  particleVector.markAsClean();

  particleVector.push_back(p6);
  particleVector.push_back(p7);

  EXPECT_EQ(particleVector.totalSize(), 7);
  // The following line cannot be guaranteed (but works for gcc + clang)! MSVC, e.g., has a different growth rate!
  // see: https://tylerayoung.com/2020/08/20/default-capacity-growth-rate-of-c-stdvector/
  EXPECT_EQ(particleVector.dirtySize(), 2);
}

TEST_F(ParticleVectorTest, testDirtySizeAfterImplicitResize) {
  auto particleVector = ParticleVector<autopas::ParticleBaseFP64>();

  autopas::ParticleBaseFP64 p1({0.5, 0.5, 0.5}, {0, 0, 0}, 0);
  autopas::ParticleBaseFP64 p2({1.5, 1.5, 1.5}, {0, 0, 0}, 1);
  autopas::ParticleBaseFP64 p3({1.6, 1.5, 1.5}, {0, 0, 0}, 2);
  autopas::ParticleBaseFP64 p4({2.5, 1.5, 1.5}, {0, 0, 0}, 3);
  autopas::ParticleBaseFP64 p5({2.5, 2.5, 2.5}, {0, 0, 0}, 4);
  autopas::ParticleBaseFP64 p6({2.5, 2.5, 2.5}, {0, 0, 0}, 4);
  autopas::ParticleBaseFP64 p7({2.5, 2.5, 2.5}, {0, 0, 0}, 4);
  autopas::ParticleBaseFP64 p8({2.5, 2.5, 2.5}, {0, 0, 0}, 4);
  autopas::ParticleBaseFP64 p9({2.5, 2.5, 2.5}, {0, 0, 0}, 4);

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
