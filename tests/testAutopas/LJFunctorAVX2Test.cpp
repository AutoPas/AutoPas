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

  ljFunctor.SoALoader(cell1NoAVX2, cell2NoAVX2._particleSoABuffer);
  ljFunctor.SoALoader(cell2NoAVX2, cell2NoAVX2._particleSoABuffer);
  ljFunctorAVX2.SoALoader(cell1AVX2, cell1AVX2._particleSoABuffer);
  ljFunctorAVX2.SoALoader(cell2AVX2, cell2AVX2._particleSoABuffer);

  ljFunctor.SoAFunctor(cell1NoAVX2._particleSoABuffer, cell2NoAVX2._particleSoABuffer, true);
  ljFunctorAVX2.SoAFunctor(cell1AVX2._particleSoABuffer, cell2AVX2._particleSoABuffer, true);

  ljFunctor.SoAExtractor(cell1NoAVX2, cell1NoAVX2._particleSoABuffer);
  ljFunctor.SoAExtractor(cell2NoAVX2, cell2NoAVX2._particleSoABuffer);
  ljFunctorAVX2.SoAExtractor(cell1AVX2, cell1AVX2._particleSoABuffer);
  ljFunctorAVX2.SoAExtractor(cell2AVX2, cell1AVX2._particleSoABuffer);

  // @TODO: Assert that at least something happened

  ASSERT_NE(cell1AVX2._particles.size(), 0);
  ASSERT_NE(cell2AVX2._particles.size(), 0);
  ASSERT_EQ(cell1NoAVX2._particles.size(), cell1AVX2._particles.size());
  ASSERT_EQ(cell2NoAVX2._particles.size(), cell2AVX2._particles.size());

  // sort all cells' particles by id
  std::sort(cell1NoAVX2._particles.begin(), cell1NoAVX2._particles.end(),
            [](const Particle &p1, const Particle &p2) { return p1.getID() < p2.getID(); });
  std::sort(cell2NoAVX2._particles.begin(), cell2NoAVX2._particles.end(),
            [](const Particle &p1, const Particle &p2) { return p1.getID() < p2.getID(); });
  std::sort(cell1AVX2._particles.begin(), cell1AVX2._particles.end(),
            [](const Particle &p1, const Particle &p2) { return p1.getID() < p2.getID(); });
  std::sort(cell2AVX2._particles.begin(), cell2AVX2._particles.end(),
            [](const Particle &p1, const Particle &p2) { return p1.getID() < p2.getID(); });

  for (auto &&pIterNoAVX2 = cell1NoAVX2._particles.begin(),
              pIterAVX2 = cell1AVX2._particles.begin();
       pIterNoAVX2 < cell1NoAVX2._particles.end(); ++pIterAVX2, ++pIterNoAVX2) {
    EXPECT_EQ(pIterAVX2->getID(), pIterNoAVX2->getID());

    EXPECT_EQ(pIterAVX2->getR()[0], pIterNoAVX2->getR()[0]);
    EXPECT_EQ(pIterAVX2->getR()[1], pIterNoAVX2->getR()[1]);
    EXPECT_EQ(pIterAVX2->getR()[2], pIterNoAVX2->getR()[2]);

    EXPECT_EQ(pIterAVX2->getV()[0], pIterNoAVX2->getV()[0]);
    EXPECT_EQ(pIterAVX2->getV()[1], pIterNoAVX2->getV()[1]);
    EXPECT_EQ(pIterAVX2->getV()[2], pIterNoAVX2->getV()[2]);

    EXPECT_EQ(pIterAVX2->getF()[0], pIterNoAVX2->getF()[0]);
    EXPECT_EQ(pIterAVX2->getF()[1], pIterNoAVX2->getF()[1]);
    EXPECT_EQ(pIterAVX2->getF()[2], pIterNoAVX2->getF()[2]);
  }
}
