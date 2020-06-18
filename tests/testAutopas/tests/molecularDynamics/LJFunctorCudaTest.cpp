/**
 * @file LJFunctorCudaTest.cpp
 * @author jspahl
 * @date 2/11/18
 */

#if defined(AUTOPAS_CUDA)

#include "LJFunctorCudaTest.h"

#include "autopas/cells/FullParticleCell.h"
#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopas/utils/StaticSelectorMacros.h"
#include "autopasTools/generators/RandomGenerator.h"

template <class Particle>
bool LJFunctorCudaTest::SoAParticlesEqual(autopas::SoA<typename Particle::SoAArraysType> &soa1,
                                          autopas::SoA<typename Particle::SoAArraysType> &soa2) {
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

bool LJFunctorCudaTest::AoSParticlesEqual(FMCell &cell1, FMCell &cell2) {
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
  FMCell cell1Cuda;
  FMCell cell2Cuda;

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::RandomGenerator::fillWithParticles(
      cell1Cuda, defaultParticle, _lowCorner, {_highCorner[0] / 2, _highCorner[1], _highCorner[2]}, numParticles);
  autopasTools::generators::RandomGenerator::fillWithParticles(
      cell2Cuda, defaultParticle, {_highCorner[0] / 2, _lowCorner[1], _lowCorner[2]}, _highCorner, numParticles2);

  // copy cells
  FMCell cell1NoCuda(cell1Cuda);
  FMCell cell2NoCuda(cell2Cuda);

  constexpr bool shifting = true;
  constexpr bool mixing = false;
  autopas::LJFunctor<Molecule, FMCell, shifting, mixing, autopas::FunctorN3Modes::Both, calculateGlobals>
      ljFunctorNoCuda(_cutoff);
  ljFunctorNoCuda.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  autopas::LJFunctor<Molecule, FMCell, shifting, mixing, autopas::FunctorN3Modes::Both, calculateGlobals> ljFunctorCuda(
      _cutoff);
  ljFunctorCuda.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

  ljFunctorCuda.getCudaWrapper()->setNumThreads(32);

  ljFunctorNoCuda.initTraversal();
  ljFunctorCuda.initTraversal();

  ASSERT_TRUE(AoSParticlesEqual(cell1Cuda, cell1NoCuda)) << "Cells 1 not equal after copy initialization.";
  ASSERT_TRUE(AoSParticlesEqual(cell2Cuda, cell2NoCuda)) << "Cells 2 not equal after copy initialization.";

  ljFunctorNoCuda.SoALoader(cell1NoCuda, cell1NoCuda._particleSoABuffer, 0);
  ljFunctorNoCuda.SoALoader(cell2NoCuda, cell2NoCuda._particleSoABuffer, 0);
  ljFunctorCuda.SoALoader(cell1Cuda, cell1Cuda._particleSoABuffer, 0);
  ljFunctorCuda.SoALoader(cell2Cuda, cell2Cuda._particleSoABuffer, 0);

  ljFunctorCuda.deviceSoALoader(cell1Cuda._particleSoABuffer, cell1Cuda._particleSoABufferDevice);
  ljFunctorCuda.deviceSoALoader(cell2Cuda._particleSoABuffer, cell2Cuda._particleSoABufferDevice);

  // Functor calls
  ljFunctorNoCuda.SoAFunctorPair(cell1NoCuda._particleSoABuffer, cell2NoCuda._particleSoABuffer, useNewton3);
  ljFunctorCuda.CudaFunctor(cell1Cuda._particleSoABufferDevice, cell2Cuda._particleSoABufferDevice, useNewton3);

  ljFunctorCuda.deviceSoAExtractor(cell1Cuda._particleSoABuffer, cell1Cuda._particleSoABufferDevice);
  ljFunctorCuda.deviceSoAExtractor(cell2Cuda._particleSoABuffer, cell2Cuda._particleSoABufferDevice);

  ASSERT_TRUE(SoAParticlesEqual<Molecule>(cell1Cuda._particleSoABuffer, cell1NoCuda._particleSoABuffer))
      << "Cells 1 not equal after applying functor and extracting to SoA.";
  ASSERT_TRUE(SoAParticlesEqual<Molecule>(cell2Cuda._particleSoABuffer, cell2NoCuda._particleSoABuffer))
      << "Cells 2 not equal after applying functor and extracting to SoA.";

  ljFunctorCuda.SoAExtractor(cell1Cuda, cell1Cuda._particleSoABuffer, 0);
  ljFunctorCuda.SoAExtractor(cell2Cuda, cell2Cuda._particleSoABuffer, 0);
  ljFunctorCuda.SoAExtractor(cell1NoCuda, cell1NoCuda._particleSoABuffer, 0);
  ljFunctorCuda.SoAExtractor(cell2NoCuda, cell2NoCuda._particleSoABuffer, 0);

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
  FMCell cellCuda;

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::RandomGenerator::fillWithParticles(cellCuda, defaultParticle, _lowCorner, _highCorner,
                                                               numParticles);

  // copy cells
  FMCell cellNoCuda(cellCuda);
  constexpr bool shifting = true;
  constexpr bool mixing = false;
  autopas::LJFunctor<Molecule, FMCell, shifting, mixing, autopas::FunctorN3Modes::Both, calculateGlobals>
      ljFunctorNoCuda(_cutoff);
  ljFunctorNoCuda.setParticleProperties(sqrt(_epsilon * _epsilon) * 24.0,
                                        ((_sigma + _sigma) / 2) * (_sigma + _sigma) / 2);
  autopas::LJFunctor<Molecule, FMCell, shifting, mixing, autopas::FunctorN3Modes::Both, calculateGlobals> ljFunctorCuda(
      _cutoff);
  ljFunctorCuda.setParticleProperties(sqrt(_epsilon * _epsilon) * 24.0,
                                      ((_sigma + _sigma) / 2) * (_sigma + _sigma) / 2);
  ljFunctorCuda.getCudaWrapper()->setNumThreads(32);

  ljFunctorNoCuda.initTraversal();
  ljFunctorCuda.initTraversal();

  ASSERT_TRUE(AoSParticlesEqual(cellCuda, cellNoCuda)) << "Cells not equal after copy initialization.";

  ljFunctorNoCuda.SoALoader(cellNoCuda, cellNoCuda._particleSoABuffer, 0);
  ljFunctorCuda.SoALoader(cellCuda, cellCuda._particleSoABuffer, 0);
  ljFunctorCuda.deviceSoALoader(cellCuda._particleSoABuffer, cellCuda._particleSoABufferDevice);

  // functor calls
  ljFunctorNoCuda.SoAFunctorSingle(cellNoCuda._particleSoABuffer, useNewton3);
  ljFunctorCuda.CudaFunctor(cellCuda._particleSoABufferDevice, useNewton3);

  ljFunctorCuda.deviceSoAExtractor(cellCuda._particleSoABuffer, cellCuda._particleSoABufferDevice);

  ASSERT_TRUE(SoAParticlesEqual<Molecule>(cellCuda._particleSoABuffer, cellNoCuda._particleSoABuffer))
      << "Cells not equal after applying functor and extracting to SoA.";

  ljFunctorCuda.SoAExtractor(cellCuda, cellCuda._particleSoABuffer, 0);
  ljFunctorCuda.SoAExtractor(cellNoCuda, cellNoCuda._particleSoABuffer, 0);

  ljFunctorNoCuda.endTraversal(useNewton3);
  ljFunctorCuda.endTraversal(useNewton3);

  ASSERT_TRUE(AoSParticlesEqual(cellCuda, cellNoCuda)) << "Cells not equal after extracting.";

  if (calculateGlobals) {
    ASSERT_NEAR(ljFunctorNoCuda.getUpot(), ljFunctorCuda.getUpot(), 1.0e-13);
    ASSERT_NEAR(ljFunctorNoCuda.getVirial(), ljFunctorCuda.getVirial(), 1.0e-13);
  }
}

TEST_P(LJFunctorCudaTest, testLJFunctorVSLJFunctorCuda) {
  /// @todo c++20: replace below std::get<> with structured bindings. c++20 allows captured values in lambdas.
  // auto [newton3, calculateGlobals, numParticlesFirstCell, numParticlesSecondCell] = GetParam();
  auto options = GetParam();
  auto newton3 = std::get<0>(options);
  auto calculateGlobals = std::get<1>(options);
  auto numParticlesFirstCell = std::get<2>(options);
  auto numParticlesSecondCell = std::get<3>(options);

  // using nested withStaticBool is not possible because of bug in gcc < 9 (and the intel compiler)
  /// @todo c++20: gcc < 9 can probably be dropped, replace with nested lambdas.
  if (newton3) {
    autopas::utils::withStaticBool(calculateGlobals, [&](auto calculateGlobalsC) {
      if (numParticlesSecondCell == 0) {
        testLJFunctorVSLJFunctorCudaOneCell<Particle, true /*newton3*/, calculateGlobalsC>(numParticlesFirstCell);
      } else {
        testLJFunctorVSLJFunctorCudaTwoCells<Particle, true /*newton3*/, calculateGlobalsC>(numParticlesFirstCell,
                                                                                            numParticlesSecondCell);
      }
    });
  } else {
    autopas::utils::withStaticBool(calculateGlobals, [&](auto calculateGlobalsC) {
      if (numParticlesSecondCell == 0) {
        testLJFunctorVSLJFunctorCudaOneCell<Particle, false /*newton3*/, calculateGlobalsC>(numParticlesFirstCell);
      } else {
        testLJFunctorVSLJFunctorCudaTwoCells<Particle, false /*newton3*/, calculateGlobalsC>(numParticlesFirstCell,
                                                                                             numParticlesSecondCell);
      }
    });
  }
}

static auto toString = [](const auto &info) {
  auto [newton3, calculateGlobals, numParticlesFirstCell, numParticlesSecondCell] = info.param;
  std::stringstream resStream;
  resStream << (newton3 ? "N3" : "noN3") << (calculateGlobals ? "globals" : "noGlobals") << numParticlesFirstCell << "x"
            << numParticlesSecondCell;
  std::string res = resStream.str();
  std::replace(res.begin(), res.end(), '-', '_');
  std::replace(res.begin(), res.end(), '.', '_');
  return res;
};

INSTANTIATE_TEST_SUITE_P(
    Generated, LJFunctorCudaTest,
    ::testing::Combine(::testing::Bool(), ::testing::Bool(),
                       ::testing::ValuesIn({1, 2, 4, 16, 31, 32, 33, 55, 64, 65}) /* numParticlesFirstCell */,
                       ::testing::ValuesIn({0, 1, 4, 16, 31, 32, 33, 55, 64, 65}) /* numParticlesSecondCell */),
    toString);

#endif  // AUTOPAS_CUDA
