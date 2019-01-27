/**
 * @file LJFunctorAVX2Test.cpp
 * @author F. Gratl
 * @date 12/17/18
 */

#include "LJFunctorAVX2Test.h"
#include <autopas/cells/FullParticleCell.h>
#include <autopas/pairwiseFunctors/LJFunctorAVX2.h>
#include <autopas/particles/Particle.h>
#include <testingHelpers/RandomGenerator.h>

template <class SoAType>
bool LJFunctorAVX2Test::SoAParticlesEqual(autopas::SoA<SoAType> &soa1, autopas::SoA<SoAType> &soa2) {
  EXPECT_GT(soa1.getNumParticles(), 0);
  EXPECT_EQ(soa1.getNumParticles(), soa2.getNumParticles());

  unsigned long *const __restrict__ idptr1 = soa1.template begin<Particle::AttributeNames::id>();
  unsigned long *const __restrict__ idptr2 = soa2.template begin<Particle::AttributeNames::id>();

  double *const __restrict__ xptr1 = soa1.template begin<Particle::AttributeNames::posX>();
  double *const __restrict__ yptr1 = soa1.template begin<Particle::AttributeNames::posY>();
  double *const __restrict__ zptr1 = soa1.template begin<Particle::AttributeNames::posZ>();
  double *const __restrict__ xptr2 = soa2.template begin<Particle::AttributeNames::posX>();
  double *const __restrict__ yptr2 = soa2.template begin<Particle::AttributeNames::posY>();
  double *const __restrict__ zptr2 = soa2.template begin<Particle::AttributeNames::posZ>();

  double *const __restrict__ fxptr1 = soa1.template begin<Particle::AttributeNames::forceX>();
  double *const __restrict__ fyptr1 = soa1.template begin<Particle::AttributeNames::forceY>();
  double *const __restrict__ fzptr1 = soa1.template begin<Particle::AttributeNames::forceZ>();
  double *const __restrict__ fxptr2 = soa2.template begin<Particle::AttributeNames::forceX>();
  double *const __restrict__ fyptr2 = soa2.template begin<Particle::AttributeNames::forceY>();
  double *const __restrict__ fzptr2 = soa2.template begin<Particle::AttributeNames::forceZ>();

  for (size_t i = 0; i < soa1.getNumParticles(); ++i) {
    EXPECT_EQ(idptr1[i], idptr2[i]);
    // TODO: assert near
    EXPECT_DOUBLE_EQ(xptr1[i], xptr2[i]) << "for particle pair " << idptr1[i];
    EXPECT_DOUBLE_EQ(yptr1[i], yptr2[i]) << "for particle pair " << idptr1[i];
    EXPECT_DOUBLE_EQ(zptr1[i], zptr2[i]) << "for particle pair " << idptr1[i];
    EXPECT_DOUBLE_EQ(fxptr1[i], fxptr2[i]) << "for particle pair " << idptr1[i];
    EXPECT_DOUBLE_EQ(fyptr1[i], fyptr2[i]) << "for particle pair " << idptr1[i];
    EXPECT_DOUBLE_EQ(fzptr1[i], fzptr2[i]) << "for particle pair " << idptr1[i];
  }
  return not::testing::Test::HasFailure();
}

bool LJFunctorAVX2Test::particleEqual(Particle &p1, Particle &p2) {
  EXPECT_EQ(p1.getID(), p2.getID());

  EXPECT_DOUBLE_EQ(p1.getR()[0], p2.getR()[0]) << "for particle pair " << p1.getID();
  EXPECT_DOUBLE_EQ(p1.getR()[1], p2.getR()[1]) << "for particle pair " << p1.getID();
  EXPECT_DOUBLE_EQ(p1.getR()[2], p2.getR()[2]) << "for particle pair " << p1.getID();
  EXPECT_DOUBLE_EQ(p1.getF()[0], p2.getF()[0]) << "for particle pair " << p1.getID();
  EXPECT_DOUBLE_EQ(p1.getF()[1], p2.getF()[1]) << "for particle pair " << p1.getID();
  EXPECT_DOUBLE_EQ(p1.getF()[2], p2.getF()[2]) << "for particle pair " << p1.getID();

  return not::testing::Test::HasFailure();
}

bool LJFunctorAVX2Test::AoSParticlesEqual(FPCell &cell1, FPCell &cell2) {
  EXPECT_GT(cell1.numParticles(), 0);
  EXPECT_EQ(cell1.numParticles(), cell2.numParticles());

  bool ret = true;
  for (size_t i = 0; i < cell1.numParticles(); ++i) {
    ret &= particleEqual(cell1._particles[i], cell2._particles[i]);
  }

  return ret;
}

void LJFunctorAVX2Test::testLJFunctorVSLJFunctorAVXTwoCells(bool newton3) {
  FPCell cell1AVX2;
  FPCell cell2AVX2;

  size_t numParticles = 7;

  Particle defaultParticle({0, 0, 0}, {0, 0, 0}, 0);
  RandomGenerator::fillWithParticles(cell1AVX2, defaultParticle, _lowCorner,
                                     {_highCorner[0] / 2, _highCorner[1], _highCorner[2]}, numParticles);
  RandomGenerator::fillWithParticles(cell2AVX2, defaultParticle, {_highCorner[0] / 2, _lowCorner[1], _lowCorner[2]},
                                     _highCorner, numParticles);

  // copy cells
  FPCell cell1NoAVX2(cell1AVX2);
  FPCell cell2NoAVX2(cell2AVX2);

  autopas::LJFunctor<Particle, FPCell, false, false> ljFunctorNoAVX2(_cutoff, _epsilon, _sigma, 0.0, _lowCorner,
                                                                     _highCorner);
  autopas::LJFunctorAVX2<Particle, FPCell, false, false> ljFunctorAVX2(_cutoff, _epsilon, _sigma, 0.0, _lowCorner,
                                                                       _highCorner);

  ASSERT_TRUE(AoSParticlesEqual(cell1AVX2, cell1NoAVX2)) << "Cells 1 not equal after copy initialization.";
  ASSERT_TRUE(AoSParticlesEqual(cell2AVX2, cell2NoAVX2)) << "Cells 2 not equal after copy initialization.";

  ljFunctorNoAVX2.SoALoader(cell1NoAVX2, cell1NoAVX2._particleSoABuffer);
  ljFunctorNoAVX2.SoALoader(cell2NoAVX2, cell2NoAVX2._particleSoABuffer);
  ljFunctorAVX2.SoALoader(cell1AVX2, cell1AVX2._particleSoABuffer);
  ljFunctorAVX2.SoALoader(cell2AVX2, cell2AVX2._particleSoABuffer);

  ASSERT_TRUE(SoAParticlesEqual(cell1AVX2._particleSoABuffer, cell1NoAVX2._particleSoABuffer))
      << "Cells 1 not equal after loading.";
  ASSERT_TRUE(SoAParticlesEqual(cell2AVX2._particleSoABuffer, cell2NoAVX2._particleSoABuffer))
      << "Cells 2 not equal after loading.";

  ljFunctorNoAVX2.SoAFunctor(cell1NoAVX2._particleSoABuffer, cell2NoAVX2._particleSoABuffer, newton3);
  ljFunctorAVX2.SoAFunctor(cell1AVX2._particleSoABuffer, cell2AVX2._particleSoABuffer, newton3);

  ASSERT_TRUE(SoAParticlesEqual(cell1AVX2._particleSoABuffer, cell1NoAVX2._particleSoABuffer))
      << "Cells 1 not equal after applying functor.";
  ASSERT_TRUE(SoAParticlesEqual(cell2AVX2._particleSoABuffer, cell2NoAVX2._particleSoABuffer))
      << "Cells 2 not equal after applying functor.";

  ljFunctorAVX2.SoAExtractor(cell1AVX2, cell1AVX2._particleSoABuffer);
  ljFunctorAVX2.SoAExtractor(cell2AVX2, cell2AVX2._particleSoABuffer);
  ljFunctorAVX2.SoAExtractor(cell1NoAVX2, cell1NoAVX2._particleSoABuffer);
  ljFunctorAVX2.SoAExtractor(cell2NoAVX2, cell2NoAVX2._particleSoABuffer);

  ASSERT_TRUE(AoSParticlesEqual(cell1AVX2, cell1NoAVX2)) << "Cells 1 not equal after extracting.";
  ASSERT_TRUE(AoSParticlesEqual(cell2AVX2, cell2NoAVX2)) << "Cells 2 not equal after extracting.";
}

void LJFunctorAVX2Test::testLJFunctorVSLJFunctorAVXOneCell(bool newton3) {
  FPCell cellAVX2;

  size_t numParticles = 7;

  Particle defaultParticle({0, 0, 0}, {0, 0, 0}, 0);
  RandomGenerator::fillWithParticles(cellAVX2, defaultParticle, _lowCorner, _highCorner, numParticles);

  // copy cells
  FPCell cellNoAVX2(cellAVX2);

  autopas::LJFunctor<Particle, FPCell, false, false> ljFunctorNoAVX2(_cutoff, _epsilon, _sigma, 0.0, _lowCorner,
                                                                     _highCorner);
  autopas::LJFunctorAVX2<Particle, FPCell, false, false> ljFunctorAVX2(_cutoff, _epsilon, _sigma, 0.0, _lowCorner,
                                                                       _highCorner);

  ASSERT_TRUE(AoSParticlesEqual(cellAVX2, cellNoAVX2)) << "Cells not equal after copy initialization.";

  ljFunctorNoAVX2.SoALoader(cellNoAVX2, cellNoAVX2._particleSoABuffer);
  ljFunctorAVX2.SoALoader(cellAVX2, cellAVX2._particleSoABuffer);

  ASSERT_TRUE(SoAParticlesEqual(cellAVX2._particleSoABuffer, cellNoAVX2._particleSoABuffer))
      << "Cells not equal after loading.";

  ljFunctorNoAVX2.SoAFunctor(cellNoAVX2._particleSoABuffer, newton3);
  ljFunctorAVX2.SoAFunctor(cellAVX2._particleSoABuffer, newton3);

  ASSERT_TRUE(SoAParticlesEqual(cellAVX2._particleSoABuffer, cellNoAVX2._particleSoABuffer))
      << "Cells not equal after applying functor.";

  ljFunctorAVX2.SoAExtractor(cellAVX2, cellAVX2._particleSoABuffer);
  ljFunctorAVX2.SoAExtractor(cellNoAVX2, cellNoAVX2._particleSoABuffer);

  ASSERT_TRUE(AoSParticlesEqual(cellAVX2, cellNoAVX2)) << "Cells 1 not equal after extracting.";
}

TEST_F(LJFunctorAVX2Test, testLJFunctorVSLJFunctorAVX2OneCellNewton3) { testLJFunctorVSLJFunctorAVXOneCell(true); }

TEST_F(LJFunctorAVX2Test, testLJFunctorVSLJFunctorAVX2OneCellNoNewton3) { testLJFunctorVSLJFunctorAVXOneCell(false); }

TEST_F(LJFunctorAVX2Test, testLJFunctorVSLJFunctorAVX2TwoCellNewton3) { testLJFunctorVSLJFunctorAVXTwoCells(true); }

TEST_F(LJFunctorAVX2Test, testLJFunctorVSLJFunctorAVX2TwoCellNoNewton3) { testLJFunctorVSLJFunctorAVXTwoCells(false); }