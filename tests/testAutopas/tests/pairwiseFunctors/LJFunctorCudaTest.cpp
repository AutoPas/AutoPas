/**
 * @file LJFunctorCudaTest.cpp
 * @author jspahl
 * @date 2/11/18
 */

#if defined(AUTOPAS_CUDA)

#include "LJFunctorCudaTest.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/pairwiseFunctors/LJFunctor.h"
#include "testingHelpers/RandomGenerator.h"

template <class SoAType>
bool LJFunctorCudaTest::SoAParticlesEqual(autopas::SoA<SoAType> &soa1, autopas::SoA<SoAType> &soa2) {
  EXPECT_GT(soa1.getNumParticles(), 0);
  EXPECT_EQ(soa1.getNumParticles(), soa2.getNumParticles());

  unsigned long *const __restrict__ idptr1 = soa1.template begin<Particle::AttributeNames::id>();
  unsigned long *const __restrict__ idptr2 = soa2.template begin<Particle::AttributeNames::id>();

  auto *const __restrict__ xptr1 = soa1.template begin<Particle::AttributeNames::posX>();
  auto *const __restrict__ yptr1 = soa1.template begin<Particle::AttributeNames::posY>();
  auto *const __restrict__ zptr1 = soa1.template begin<Particle::AttributeNames::posZ>();
  auto *const __restrict__ xptr2 = soa2.template begin<Particle::AttributeNames::posX>();
  auto *const __restrict__ yptr2 = soa2.template begin<Particle::AttributeNames::posY>();
  auto *const __restrict__ zptr2 = soa2.template begin<Particle::AttributeNames::posZ>();

  auto *const __restrict__ fxptr1 = soa1.template begin<Particle::AttributeNames::forceX>();
  auto *const __restrict__ fyptr1 = soa1.template begin<Particle::AttributeNames::forceY>();
  auto *const __restrict__ fzptr1 = soa1.template begin<Particle::AttributeNames::forceZ>();
  auto *const __restrict__ fxptr2 = soa2.template begin<Particle::AttributeNames::forceX>();
  auto *const __restrict__ fyptr2 = soa2.template begin<Particle::AttributeNames::forceY>();
  auto *const __restrict__ fzptr2 = soa2.template begin<Particle::AttributeNames::forceZ>();

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

template <typename ParticleType, bool useNewton3, bool calculateGlobals>
void LJFunctorCudaTest::testLJFunctorVSLJFunctorCudaTwoCells(size_t numParticles, size_t numParticles2) {
  autopas::FullParticleCell<ParticleType> cell1Cuda;
  autopas::FullParticleCell<ParticleType> cell2Cuda;

  ParticleType defaultParticle({0, 0, 0}, {0, 0, 0}, 0);
  RandomGenerator::fillWithParticles(cell1Cuda, defaultParticle, _lowCorner,
                                     {_highCorner[0] / 2, _highCorner[1], _highCorner[2]}, numParticles);
  RandomGenerator::fillWithParticles(cell2Cuda, defaultParticle, {_highCorner[0] / 2, _lowCorner[1], _lowCorner[2]},
                                     _highCorner, numParticles2);

  // copy cells
  autopas::FullParticleCell<ParticleType> cell1NoCuda(cell1Cuda);
  autopas::FullParticleCell<ParticleType> cell2NoCuda(cell2Cuda);

  autopas::LJFunctor<Particle, autopas::FullParticleCell<ParticleType>, autopas::FunctorN3Modes::Both, calculateGlobals>
      ljFunctorNoCuda(_cutoff, _epsilon, _sigma, 0.0);
  autopas::LJFunctor<Particle, autopas::FullParticleCell<ParticleType>, autopas::FunctorN3Modes::Both, calculateGlobals>
      ljFunctorCuda(_cutoff, _epsilon, _sigma, 0.0);

  ljFunctorCuda.getCudaWrapper()->setNumThreads(32);

  ljFunctorNoCuda.initTraversal();
  ljFunctorCuda.initTraversal();

  ASSERT_TRUE(AoSParticlesEqual(cell1Cuda, cell1NoCuda)) << "Cells 1 not equal after copy initialization.";
  ASSERT_TRUE(AoSParticlesEqual(cell2Cuda, cell2NoCuda)) << "Cells 2 not equal after copy initialization.";

  ljFunctorNoCuda.SoALoader(cell1NoCuda, cell1NoCuda._particleSoABuffer);
  ljFunctorNoCuda.SoALoader(cell2NoCuda, cell2NoCuda._particleSoABuffer);
  ljFunctorCuda.SoALoader(cell1Cuda, cell1Cuda._particleSoABuffer);
  ljFunctorCuda.SoALoader(cell2Cuda, cell2Cuda._particleSoABuffer);

  ljFunctorCuda.deviceSoALoader(cell1Cuda._particleSoABuffer, cell1Cuda._particleSoABufferDevice);
  ljFunctorCuda.deviceSoALoader(cell2Cuda._particleSoABuffer, cell2Cuda._particleSoABufferDevice);

  // Functor calls
  ljFunctorNoCuda.SoAFunctor(cell1NoCuda._particleSoABuffer, cell2NoCuda._particleSoABuffer, useNewton3);
  ljFunctorCuda.CudaFunctor(cell1Cuda._particleSoABufferDevice, cell2Cuda._particleSoABufferDevice, useNewton3);

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

  ljFunctorNoCuda.endTraversal(useNewton3);
  ljFunctorCuda.endTraversal(useNewton3);

  ASSERT_TRUE(AoSParticlesEqual(cell1Cuda, cell1NoCuda)) << "Cells 1 not equal after extracting.";
  ASSERT_TRUE(AoSParticlesEqual(cell2Cuda, cell2NoCuda)) << "Cells 2 not equal after extracting.";

  if (calculateGlobals) {
    ASSERT_NEAR(ljFunctorNoCuda.getUpot(), ljFunctorCuda.getUpot(), 1.0e-13);
    ASSERT_NEAR(ljFunctorNoCuda.getVirial(), ljFunctorCuda.getVirial(), 1.0e-13);
  }
}

template <typename ParticleType, bool useNewton3, bool calculateGlobals>
void LJFunctorCudaTest::testLJFunctorVSLJFunctorCudaOneCell(size_t numParticles) {
  autopas::FullParticleCell<ParticleType> cellCuda;

  ParticleType defaultParticle({0, 0, 0}, {0, 0, 0}, 0);
  RandomGenerator::fillWithParticles(cellCuda, defaultParticle, _lowCorner, _highCorner, numParticles);

  // copy cells
  autopas::FullParticleCell<ParticleType> cellNoCuda(cellCuda);

  autopas::LJFunctor<Particle, autopas::FullParticleCell<ParticleType>, autopas::FunctorN3Modes::Both, calculateGlobals>
      ljFunctorNoCuda(_cutoff, _epsilon, _sigma, 0.0);
  autopas::LJFunctor<Particle, autopas::FullParticleCell<ParticleType>, autopas::FunctorN3Modes::Both, calculateGlobals>
      ljFunctorCuda(_cutoff, _epsilon, _sigma, 0.0);
  ljFunctorCuda.getCudaWrapper()->setNumThreads(32);

  ljFunctorNoCuda.initTraversal();
  ljFunctorCuda.initTraversal();

  ASSERT_TRUE(AoSParticlesEqual(cellCuda, cellNoCuda)) << "Cells not equal after copy initialization.";

  ljFunctorNoCuda.SoALoader(cellNoCuda, cellNoCuda._particleSoABuffer);
  ljFunctorCuda.SoALoader(cellCuda, cellCuda._particleSoABuffer);
  ljFunctorCuda.deviceSoALoader(cellCuda._particleSoABuffer, cellCuda._particleSoABufferDevice);

  // functor calls
  ljFunctorNoCuda.SoAFunctor(cellNoCuda._particleSoABuffer, useNewton3);
  ljFunctorCuda.CudaFunctor(cellCuda._particleSoABufferDevice, useNewton3);

  ljFunctorCuda.deviceSoAExtractor(cellCuda._particleSoABuffer, cellCuda._particleSoABufferDevice);

  ASSERT_TRUE(SoAParticlesEqual(cellCuda._particleSoABuffer, cellNoCuda._particleSoABuffer))
      << "Cells not equal after applying functor and extracting to SoA.";

  ljFunctorCuda.SoAExtractor(cellCuda, cellCuda._particleSoABuffer);
  ljFunctorCuda.SoAExtractor(cellNoCuda, cellNoCuda._particleSoABuffer);

  ljFunctorNoCuda.endTraversal(useNewton3);
  ljFunctorCuda.endTraversal(useNewton3);

  ASSERT_TRUE(AoSParticlesEqual(cellCuda, cellNoCuda)) << "Cells not equal after extracting.";

  if (calculateGlobals) {
    ASSERT_NEAR(ljFunctorNoCuda.getUpot(), ljFunctorCuda.getUpot(), 1.0e-13);
    ASSERT_NEAR(ljFunctorNoCuda.getVirial(), ljFunctorCuda.getVirial(), 1.0e-13);
  }
}

TEST_P(LJFunctorCudaTest, testLJFunctorVSLJFunctorCuda) {
  auto options = GetParam();
  if (std::get<3>(options) < 0) {
    if (std::get<0>(options)) {
      if (std::get<1>(options)) {
        testLJFunctorVSLJFunctorCudaOneCell<Particle, true, true>(std::get<2>(options));
      } else {
        testLJFunctorVSLJFunctorCudaOneCell<Particle, true, false>(std::get<2>(options));
      }
    } else {
      if (std::get<1>(options)) {
        testLJFunctorVSLJFunctorCudaOneCell<Particle, false, true>(std::get<2>(options));
      } else {
        testLJFunctorVSLJFunctorCudaOneCell<Particle, false, false>(std::get<2>(options));
      }
    }
  } else {
    if (std::get<0>(options)) {
      if (std::get<1>(options)) {
        testLJFunctorVSLJFunctorCudaTwoCells<Particle, true, true>(std::get<2>(options), std::get<3>(options));
      } else {
        testLJFunctorVSLJFunctorCudaTwoCells<Particle, true, false>(std::get<2>(options), std::get<3>(options));
      }
    } else {
      if (std::get<1>(options)) {
        testLJFunctorVSLJFunctorCudaTwoCells<Particle, false, true>(std::get<2>(options), std::get<3>(options));
      } else {
        testLJFunctorVSLJFunctorCudaTwoCells<Particle, false, false>(std::get<2>(options), std::get<3>(options));
      }
    }
  }
}

INSTANTIATE_TEST_CASE_P(Generated, LJFunctorCudaTest,
                        ::testing::Combine(::testing::Bool(), ::testing::Bool(),
                                           ::testing::ValuesIn({1, 2, 4, 16, 31, 32, 33, 55, 64, 65}),
                                           ::testing::ValuesIn({-1, 1, 4, 16, 31, 32, 33, 55, 64, 65})));

#endif  // AUTOPAS_CUDA
