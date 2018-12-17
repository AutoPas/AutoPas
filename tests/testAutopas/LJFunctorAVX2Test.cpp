/**
 * @file LJFunctorAVX2Test.cpp
 * @author F. Gratl
 * @date 12/17/18
 */

#include "LJFunctorAVX2Test.h"
#include <autopas/cells/FullParticleCell.h>
#include <autopas/particles/Particle.h>
#include <autopas/utils/SoA.h>
#include <testingHelpers/RandomGenerator.h>
#include <testingHelpers/commonTypedefs.h>

TEST_F(LJFunctorAVX2Test, testLJFunctorVSLJFunctorAVX2) {
  //  auto soa1 = autopas::Particle::SoAArraysType();
  FPCell cell1NoAVX2;
  FPCell cell2NoAVX2;
  FPCell cell1AVX2;
  FPCell cell2AVX2;

  Particle defaultParticle;
  RandomGenerator::fillWithParticles(cell1AVX2, defaultParticle, {0, 0, 0}, {1, 1, 1}, 8);
  RandomGenerator::fillWithParticles(cell2AVX2, defaultParticle, {1, 0, 0}, {2, 1, 1}, 8);
}
