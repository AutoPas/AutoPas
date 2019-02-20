/**
 * @file LJFunctorCudaTest.cpp
 * @author jspahl
 * @date 2/11/18
 */

#if defined(AUTOPAS_CUDA)

#include "LJFunctorCudaTest.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/pairwiseFunctors/LJFunctor.h"
#include "autopas/particles/Particle.h"
#include "testingHelpers/RandomGenerator.h"

template <class SoAType>
bool LJFunctorCudaTest::SoAParticlesEqual(autopas::SoA<SoAType> &soa1, autopas::SoA<SoAType> &soa2) {
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

    EXPECT_NEAR(xptr1[i], xptr2[i], _maxError) << "for particle pair " << idptr1[i];
    EXPECT_NEAR(yptr1[i], yptr2[i], _maxError) << "for particle pair " << idptr1[i];
    EXPECT_NEAR(zptr1[i], zptr2[i], _maxError) << "for particle pair " << idptr1[i];
    EXPECT_NEAR(fxptr1[i], fxptr2[i], _maxError) << "for particle pair " << idptr1[i];
    EXPECT_NEAR(fyptr1[i], fyptr2[i], _maxError) << "for particle pair " << idptr1[i];
    EXPECT_NEAR(fzptr1[i], fzptr2[i], _maxError) << "for particle pair " << idptr1[i];
  }
  // clang-format off
	return not ::testing::Test::HasFailure();
  // clang-format on
}

bool LJFunctorCudaTest::particleEqual(Particle &p1, Particle &p2) {
  EXPECT_EQ(p1.getID(), p2.getID());

  EXPECT_NEAR(p1.getR()[0], p2.getR()[0], _maxError) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getR()[1], p2.getR()[1], _maxError) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getR()[2], p2.getR()[2], _maxError) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getF()[0], p2.getF()[0], _maxError) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getF()[1], p2.getF()[1], _maxError) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getF()[2], p2.getF()[2], _maxError) << "for particle pair " << p1.getID();

  // clang-format off
	return not ::testing::Test::HasFailure();
  // clang-format on
}

bool LJFunctorCudaTest::AoSParticlesEqual(FPCell &cell1, FPCell &cell2) {
  EXPECT_GT(cell1.numParticles(), 0);
  EXPECT_EQ(cell1.numParticles(), cell2.numParticles());

  bool ret = true;
  for (size_t i = 0; i < cell1.numParticles(); ++i) {
    ret &= particleEqual(cell1._particles[i], cell2._particles[i]);
  }

  return ret;
}

void LJFunctorCudaTest::testLJFunctorVSLJFunctorCudaTwoCells(size_t numParticles, size_t numParticles2, bool newton3) {
  FPCell cell1Cuda;
  FPCell cell2Cuda;

  Particle defaultParticle({0, 0, 0}, {0, 0, 0}, 0);
  RandomGenerator::fillWithParticles(cell1Cuda, defaultParticle, _lowCorner,
                                     {_highCorner[0] / 2, _highCorner[1], _highCorner[2]}, numParticles);
  RandomGenerator::fillWithParticles(cell2Cuda, defaultParticle, {_highCorner[0] / 2, _lowCorner[1], _lowCorner[2]},
                                     _highCorner, numParticles2);

  // copy cells
  FPCell cell1NoCuda(cell1Cuda);
  FPCell cell2NoCuda(cell2Cuda);

  autopas::LJFunctor<Particle, FPCell, false, false> ljFunctorNoCuda(_cutoff, _epsilon, _sigma, 0.0, _lowCorner,
                                                                     _highCorner);
  autopas::LJFunctor<Particle, FPCell, false, false> ljFunctorCuda(_cutoff, _epsilon, _sigma, 0.0, _lowCorner,
                                                                   _highCorner);

  ASSERT_TRUE(AoSParticlesEqual(cell1Cuda, cell1NoCuda)) << "Cells 1 not equal after copy initialization.";
  ASSERT_TRUE(AoSParticlesEqual(cell2Cuda, cell2NoCuda)) << "Cells 2 not equal after copy initialization.";

  ljFunctorNoCuda.SoALoader(cell1NoCuda, cell1NoCuda._particleSoABuffer);
  ljFunctorNoCuda.SoALoader(cell2NoCuda, cell2NoCuda._particleSoABuffer);
  ljFunctorCuda.SoALoader(cell1Cuda, cell1Cuda._particleSoABuffer);
  ljFunctorCuda.SoALoader(cell2Cuda, cell2Cuda._particleSoABuffer);

  ljFunctorCuda.deviceSoALoader(cell1Cuda._particleSoABuffer, cell1Cuda._particleSoABufferDevice);
  ljFunctorCuda.deviceSoALoader(cell2Cuda._particleSoABuffer, cell2Cuda._particleSoABufferDevice);

  // Functor calls
  ljFunctorNoCuda.SoAFunctor(cell1NoCuda._particleSoABuffer, cell2NoCuda._particleSoABuffer, newton3);
  ljFunctorCuda.CudaFunctor(cell1Cuda._particleSoABufferDevice, cell2Cuda._particleSoABufferDevice, newton3);

  ljFunctorCuda.deviceSoAExtractor(cell1Cuda._particleSoABuffer, cell1Cuda._particleSoABufferDevice);
  ljFunctorCuda.deviceSoAExtractor(cell2Cuda._particleSoABuffer, cell2Cuda._particleSoABufferDevice);

  ASSERT_TRUE(SoAParticlesEqual(cell1Cuda._particleSoABuffer, cell1NoCuda._particleSoABuffer))
      << "Cells 1 not equal after applying functor and extracting to SoA.";
  ASSERT_TRUE(SoAParticlesEqual(cell2Cuda._particleSoABuffer, cell2NoCuda._particleSoABuffer))
      << "Cells 2 not equal after applying functor and extracting to SoA.";

  ljFunctorCuda.SoAExtractor(cell1Cuda, cell1Cuda._particleSoABuffer);
  ljFunctorCuda.SoAExtractor(cell2Cuda, cell2Cuda._particleSoABuffer);
  ljFunctorCuda.SoAExtractor(cell1NoCuda, cell1NoCuda._particleSoABuffer);
  ljFunctorCuda.SoAExtractor(cell2NoCuda, cell2NoCuda._particleSoABuffer);

  ASSERT_TRUE(AoSParticlesEqual(cell1Cuda, cell1NoCuda)) << "Cells 1 not equal after extracting.";
  ASSERT_TRUE(AoSParticlesEqual(cell2Cuda, cell2NoCuda)) << "Cells 2 not equal after extracting.";
}

void LJFunctorCudaTest::testLJFunctorVSLJFunctorCudaOneCell(size_t numParticles, bool newton3) {
  FPCell cellCuda;

  Particle defaultParticle({0, 0, 0}, {0, 0, 0}, 0);
  RandomGenerator::fillWithParticles(cellCuda, defaultParticle, _lowCorner, _highCorner, numParticles);

  // copy cells
  FPCell cellNoCuda(cellCuda);

  autopas::LJFunctor<Particle, FPCell, false, false> ljFunctorNoCuda(_cutoff, _epsilon, _sigma, 0.0, _lowCorner,
                                                                     _highCorner);
  autopas::LJFunctor<Particle, FPCell, false, false> ljFunctorCuda(_cutoff, _epsilon, _sigma, 0.0, _lowCorner,
                                                                   _highCorner);

  ASSERT_TRUE(AoSParticlesEqual(cellCuda, cellNoCuda)) << "Cells not equal after copy initialization.";

  ljFunctorNoCuda.SoALoader(cellNoCuda, cellNoCuda._particleSoABuffer);
  ljFunctorCuda.SoALoader(cellCuda, cellCuda._particleSoABuffer);
  ljFunctorCuda.deviceSoALoader(cellCuda._particleSoABuffer, cellCuda._particleSoABufferDevice);

  // functor calls
  ljFunctorNoCuda.SoAFunctor(cellNoCuda._particleSoABuffer, newton3);
  ljFunctorCuda.CudaFunctor(cellCuda._particleSoABufferDevice, newton3);

  ljFunctorCuda.deviceSoAExtractor(cellCuda._particleSoABuffer, cellCuda._particleSoABufferDevice);

  ASSERT_TRUE(SoAParticlesEqual(cellCuda._particleSoABuffer, cellNoCuda._particleSoABuffer))
      << "Cells not equal after applying functor and extracting to SoA.";

  ljFunctorCuda.SoAExtractor(cellCuda, cellCuda._particleSoABuffer);
  ljFunctorCuda.SoAExtractor(cellNoCuda, cellNoCuda._particleSoABuffer);

  ASSERT_TRUE(AoSParticlesEqual(cellCuda, cellNoCuda)) << "Cells not equal after extracting.";
}

TEST_F(LJFunctorCudaTest, testLJFunctorVSLJFunctorCudaOneCellNoNewton3_7Particles) {
  testLJFunctorVSLJFunctorCudaOneCell(7);
}
TEST_F(LJFunctorCudaTest, testLJFunctorVSLJFunctorCudaOneCellNewton3_7Particles) {
  testLJFunctorVSLJFunctorCudaOneCell(7, true);
}
TEST_F(LJFunctorCudaTest, testLJFunctorVSLJFunctorCudaTwoCellNoNewton3_7_7Particles) {
  testLJFunctorVSLJFunctorCudaTwoCells(7, 7);
}
TEST_F(LJFunctorCudaTest, testLJFunctorVSLJFunctorCudaTwoCellNewton3_7_21Particles) {
  testLJFunctorVSLJFunctorCudaTwoCells(7, 21, true);
}

TEST_F(LJFunctorCudaTest, testLJFunctorVSLJFunctorCudaTwoCellNoNewton3_7_34Particles) {
  testLJFunctorVSLJFunctorCudaTwoCells(7, 34);
}

TEST_F(LJFunctorCudaTest, testLJFunctorVSLJFunctorCudaOneCellNoNewton3_32Particles) {
  testLJFunctorVSLJFunctorCudaOneCell(32);
}
TEST_F(LJFunctorCudaTest, testLJFunctorVSLJFunctorCudaOneCellNewton3_32Particles) {
  testLJFunctorVSLJFunctorCudaOneCell(32, true);
}
TEST_F(LJFunctorCudaTest, testLJFunctorVSLJFunctorCudaTwoCellNoNewton3_32_32Particles) {
  testLJFunctorVSLJFunctorCudaTwoCells(32, 32);
}
TEST_F(LJFunctorCudaTest, testLJFunctorVSLJFunctorCudaTwoCellNewton3_32_32Particles) {
  testLJFunctorVSLJFunctorCudaTwoCells(32, 32, true);
}

TEST_F(LJFunctorCudaTest, testLJFunctorVSLJFunctorCudaOneCellNoNewton3_34Particles) {
  testLJFunctorVSLJFunctorCudaOneCell(34);
}
TEST_F(LJFunctorCudaTest, testLJFunctorVSLJFunctorCudaOneCellNewton3_34Particles) {
  testLJFunctorVSLJFunctorCudaOneCell(34, true);
}
TEST_F(LJFunctorCudaTest, testLJFunctorVSLJFunctorCudaOneCellNewton3_66Particles) {
  testLJFunctorVSLJFunctorCudaOneCell(66, true);
}
TEST_F(LJFunctorCudaTest, testLJFunctorVSLJFunctorCudaTwoCellNoNewton3_34_34Particles) {
  testLJFunctorVSLJFunctorCudaTwoCells(34, 34);
}
TEST_F(LJFunctorCudaTest, testLJFunctorVSLJFunctorCudaTwoCellNewton3_35_34Particles) {
  testLJFunctorVSLJFunctorCudaTwoCells(35, 34, true);
}
TEST_F(LJFunctorCudaTest, testLJFunctorVSLJFunctorCudaTwoCellNewton3_21_35Particles) {
  testLJFunctorVSLJFunctorCudaTwoCells(21, 35, true);
}
TEST_F(LJFunctorCudaTest, testLJFunctorVSLJFunctorCudaTwoCellNewton3_35_21Particles) {
  testLJFunctorVSLJFunctorCudaTwoCells(35, 21, true);
}
TEST_F(LJFunctorCudaTest, testLJFunctorVSLJFunctorCudaTwoCellNewton3_35_64Particles) {
  testLJFunctorVSLJFunctorCudaTwoCells(35, 64, true);
}
TEST_F(LJFunctorCudaTest, testLJFunctorVSLJFunctorCudaTwoCellNewton3_64_36Particles) {
  testLJFunctorVSLJFunctorCudaTwoCells(64, 36, true);
}
TEST_F(LJFunctorCudaTest, testLJFunctorVSLJFunctorCudaTwoCellNewton3_128_124Particles) {
  testLJFunctorVSLJFunctorCudaTwoCells(128, 124, true);
}
TEST_F(LJFunctorCudaTest, testLJFunctorVSLJFunctorCudaTwoCellNewton3_128_128Particles) {
  testLJFunctorVSLJFunctorCudaTwoCells(128, 128, true);
}
TEST_F(LJFunctorCudaTest, testLJFunctorVSLJFunctorCudaTwoCellNoNewton3_34_7Particles) {
  testLJFunctorVSLJFunctorCudaTwoCells(34, 7);
}

#endif  // AUTOPAS_CUDA
