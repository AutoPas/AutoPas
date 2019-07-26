/**
 * @file LJFunctorAVXTest.cpp
 * @author F. Gratl
 * @date 12/17/18
 */

#ifdef __AVX__

#include "LJFunctorAVXTest.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/pairwiseFunctors/LJFunctorAVX.h"
#include "autopas/particles/Particle.h"
#include "testingHelpers/RandomGenerator.h"

template <class SoAType>
bool LJFunctorAVXTest::SoAParticlesEqual(autopas::SoA<SoAType> &soa1, autopas::SoA<SoAType> &soa2) {
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

    double tolerance = 1e-8;
    EXPECT_NEAR(xptr1[i], xptr2[i], tolerance) << "for particle pair " << idptr1[i];
    EXPECT_NEAR(yptr1[i], yptr2[i], tolerance) << "for particle pair " << idptr1[i];
    EXPECT_NEAR(zptr1[i], zptr2[i], tolerance) << "for particle pair " << idptr1[i];
    EXPECT_NEAR(fxptr1[i], fxptr2[i], tolerance) << "for particle pair " << idptr1[i];
    EXPECT_NEAR(fyptr1[i], fyptr2[i], tolerance) << "for particle pair " << idptr1[i];
    EXPECT_NEAR(fzptr1[i], fzptr2[i], tolerance) << "for particle pair " << idptr1[i];
  }
  // clang-format off
  return not ::testing::Test::HasFailure();
  // clang-format on
}

bool LJFunctorAVXTest::particleEqual(Particle &p1, Particle &p2) {
  EXPECT_EQ(p1.getID(), p2.getID());

  double tolerance = 1e-8;
  EXPECT_NEAR(p1.getR()[0], p2.getR()[0], tolerance) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getR()[1], p2.getR()[1], tolerance) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getR()[2], p2.getR()[2], tolerance) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getF()[0], p2.getF()[0], tolerance) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getF()[1], p2.getF()[1], tolerance) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getF()[2], p2.getF()[2], tolerance) << "for particle pair " << p1.getID();

  // clang-format off
  return not ::testing::Test::HasFailure();
  // clang-format on
}

bool LJFunctorAVXTest::AoSParticlesEqual(FPCell &cell1, FPCell &cell2) {
  EXPECT_GT(cell1.numParticles(), 0);
  EXPECT_EQ(cell1.numParticles(), cell2.numParticles());

  bool ret = true;
  for (size_t i = 0; i < cell1.numParticles(); ++i) {
    ret &= particleEqual(cell1._particles[i], cell2._particles[i]);
  }

  return ret;
}

void LJFunctorAVXTest::testLJFunctorVSLJFunctorAVXTwoCells(bool newton3) {
  FPCell cell1AVX;
  FPCell cell2AVX;

  size_t numParticles = 7;

  Particle defaultParticle({0, 0, 0}, {0, 0, 0}, 0);
  RandomGenerator::fillWithParticles(cell1AVX, defaultParticle, _lowCorner,
                                     {_highCorner[0] / 2, _highCorner[1], _highCorner[2]}, numParticles);
  RandomGenerator::fillWithParticles(cell2AVX, defaultParticle, {_highCorner[0] / 2, _lowCorner[1], _lowCorner[2]},
                                     _highCorner, numParticles);

  // copy cells
  FPCell cell1NoAVX(cell1AVX);
  FPCell cell2NoAVX(cell2AVX);
    double universalValue=1; //epsilon=sigma=mass=1.0
    ParticleClassLibrary PCL = ParticleClassLibrary(universalValue,universalValue,universalValue);
  autopas::LJFunctor<Particle, FPCell, autopas::FunctorN3Modes::Both, true> ljFunctorNoAVX(_cutoff, PCL, 0.0);
  autopas::LJFunctorAVX<Particle, FPCell, autopas::FunctorN3Modes::Both, true> ljFunctorAVX(_cutoff, _epsilon, _sigma,
                                                                                            0.0);

  ljFunctorAVX.initTraversal();
  ljFunctorNoAVX.initTraversal();

  ASSERT_TRUE(AoSParticlesEqual(cell1AVX, cell1NoAVX)) << "Cells 1 not equal after copy initialization.";
  ASSERT_TRUE(AoSParticlesEqual(cell2AVX, cell2NoAVX)) << "Cells 2 not equal after copy initialization.";

  ljFunctorNoAVX.SoALoader(cell1NoAVX, cell1NoAVX._particleSoABuffer);
  ljFunctorNoAVX.SoALoader(cell2NoAVX, cell2NoAVX._particleSoABuffer);
  ljFunctorAVX.SoALoader(cell1AVX, cell1AVX._particleSoABuffer);
  ljFunctorAVX.SoALoader(cell2AVX, cell2AVX._particleSoABuffer);

  ASSERT_TRUE(SoAParticlesEqual(cell1AVX._particleSoABuffer, cell1NoAVX._particleSoABuffer))
      << "Cells 1 not equal after loading.";
  ASSERT_TRUE(SoAParticlesEqual(cell2AVX._particleSoABuffer, cell2NoAVX._particleSoABuffer))
      << "Cells 2 not equal after loading.";

  ljFunctorNoAVX.SoAFunctor(cell1NoAVX._particleSoABuffer, cell2NoAVX._particleSoABuffer, newton3);
  ljFunctorAVX.SoAFunctor(cell1AVX._particleSoABuffer, cell2AVX._particleSoABuffer, newton3);

  ASSERT_TRUE(SoAParticlesEqual(cell1AVX._particleSoABuffer, cell1NoAVX._particleSoABuffer))
      << "Cells 1 not equal after applying functor.";
  ASSERT_TRUE(SoAParticlesEqual(cell2AVX._particleSoABuffer, cell2NoAVX._particleSoABuffer))
      << "Cells 2 not equal after applying functor.";

  ljFunctorAVX.SoAExtractor(cell1AVX, cell1AVX._particleSoABuffer);
  ljFunctorAVX.SoAExtractor(cell2AVX, cell2AVX._particleSoABuffer);
  ljFunctorAVX.SoAExtractor(cell1NoAVX, cell1NoAVX._particleSoABuffer);
  ljFunctorAVX.SoAExtractor(cell2NoAVX, cell2NoAVX._particleSoABuffer);

  ASSERT_TRUE(AoSParticlesEqual(cell1AVX, cell1NoAVX)) << "Cells 1 not equal after extracting.";
  ASSERT_TRUE(AoSParticlesEqual(cell2AVX, cell2NoAVX)) << "Cells 2 not equal after extracting.";

  ljFunctorAVX.endTraversal(newton3);
  ljFunctorNoAVX.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorAVX.getUpot(), ljFunctorNoAVX.getUpot(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorAVX.getVirial(), ljFunctorNoAVX.getVirial(), tolerance) << "global virial";
}

void LJFunctorAVXTest::testLJFunctorVSLJFunctorAVXOneCell(bool newton3) {
  FPCell cellAVX;

  size_t numParticles = 7;

  Particle defaultParticle({0, 0, 0}, {0, 0, 0}, 0);
  RandomGenerator::fillWithParticles(cellAVX, defaultParticle, _lowCorner, _highCorner, numParticles);

  // copy cells
  FPCell cellNoAVX(cellAVX);
    double universalValue=1; //epsilon=sigma=mass=1.0
    ParticleClassLibrary PCL = ParticleClassLibrary(universalValue,universalValue,universalValue);
  autopas::LJFunctor<Particle, FPCell, autopas::FunctorN3Modes::Both, true> ljFunctorNoAVX(_cutoff, PCL, 0.0);
  autopas::LJFunctorAVX<Particle, FPCell, autopas::FunctorN3Modes::Both, true> ljFunctorAVX(_cutoff, _epsilon, _sigma,
                                                                                            0.0);

  ASSERT_TRUE(AoSParticlesEqual(cellAVX, cellNoAVX)) << "Cells not equal after copy initialization.";

  ljFunctorAVX.initTraversal();
  ljFunctorNoAVX.initTraversal();

  ljFunctorNoAVX.SoALoader(cellNoAVX, cellNoAVX._particleSoABuffer);
  ljFunctorAVX.SoALoader(cellAVX, cellAVX._particleSoABuffer);

  ASSERT_TRUE(SoAParticlesEqual(cellAVX._particleSoABuffer, cellNoAVX._particleSoABuffer))
      << "Cells not equal after loading.";

  ljFunctorNoAVX.SoAFunctor(cellNoAVX._particleSoABuffer, newton3);
  ljFunctorAVX.SoAFunctor(cellAVX._particleSoABuffer, newton3);

  ASSERT_TRUE(SoAParticlesEqual(cellAVX._particleSoABuffer, cellNoAVX._particleSoABuffer))
      << "Cells not equal after applying functor.";

  ljFunctorAVX.SoAExtractor(cellAVX, cellAVX._particleSoABuffer);
  ljFunctorAVX.SoAExtractor(cellNoAVX, cellNoAVX._particleSoABuffer);

  ASSERT_TRUE(AoSParticlesEqual(cellAVX, cellNoAVX)) << "Cells 1 not equal after extracting.";

  ljFunctorAVX.endTraversal(newton3);
  ljFunctorNoAVX.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorAVX.getUpot(), ljFunctorNoAVX.getUpot(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorAVX.getVirial(), ljFunctorNoAVX.getVirial(), tolerance) << "global virial";
}

TEST_F(LJFunctorAVXTest, testLJFunctorVSLJFunctorAVXOneCellNewton3) { testLJFunctorVSLJFunctorAVXOneCell(true); }

TEST_F(LJFunctorAVXTest, testLJFunctorVSLJFunctorAVXOneCellNoNewton3) { testLJFunctorVSLJFunctorAVXOneCell(false); }

TEST_F(LJFunctorAVXTest, testLJFunctorVSLJFunctorAVXTwoCellNewton3) { testLJFunctorVSLJFunctorAVXTwoCells(true); }

TEST_F(LJFunctorAVXTest, testLJFunctorVSLJFunctorAVXTwoCellNoNewton3) { testLJFunctorVSLJFunctorAVXTwoCells(false); }

#endif  // __AVX__
