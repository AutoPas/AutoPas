/**
 * @file LJFunctorAVX2Test.cpp
 * @author F. Gratl
 * @date 12/17/18
 */

#include "LJFunctorAVX2Test.h"
#include <autopas/cells/FullParticleCell.h>
#include <autopas/pairwiseFunctors/LJFunctorAVX2.h>
#include <autopas/particles/Particle.h>
#include <autopas/utils/SoA.h>
#include <testingHelpers/RandomGenerator.h>
#include <testingHelpers/commonTypedefs.h>

TEST_F(LJFunctorAVX2Test, testLJFunctorVSLJFunctorAVX2) {
  //  auto soa1 = autopas::Particle::SoAArraysType();
  FPCell cell1AVX2;
  FPCell cell2AVX2;

  Particle defaultParticle;
  RandomGenerator::fillWithParticles(cell1AVX2, defaultParticle, _lowCorner,
                                     {_lowCorner[0] / 2, _lowCorner[1], _lowCorner[2]}, 8);
  RandomGenerator::fillWithParticles(cell2AVX2, defaultParticle, {_lowCorner[0] / 2, _lowCorner[1], _lowCorner[2]},
                                     _highCorner, 8);

  // copy cells
  FPCell cell1NoAVX2(cell1AVX2);
  FPCell cell2NoAVX2(cell2AVX2);

  autopas::LJFunctor<Particle, FPCell, false, false> ljFunctor(_cutoff, _epsilon, _sigma, 0.0, _lowCorner, _highCorner);
  autopas::LJFunctorAVX2<Particle, FPCell, false, false> ljFunctorAVX2(_cutoff, _epsilon, _sigma, 0.0, _lowCorner,
                                                                       _highCorner);

  ljFunctor.S
}
