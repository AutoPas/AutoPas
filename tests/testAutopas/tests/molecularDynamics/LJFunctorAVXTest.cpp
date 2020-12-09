/**
 * @file LJFunctorAVXTest.cpp
 * @author F. Gratl
 * @date 12/17/18
 */

#ifdef __AVX__

#include "LJFunctorAVXTest.h"

#include "autopas/cells/FullParticleCell.h"
#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopas/molecularDynamics/LJFunctorAVX.h"
#include "autopas/particles/Particle.h"
#include "autopasTools/generators/RandomGenerator.h"

namespace LJFunctorAVXTest {

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
    EXPECT_NEAR(xptr1[i], xptr2[i], tolerance) << "for particle pair " << idptr1[i] << "and i=" << i;
    EXPECT_NEAR(yptr1[i], yptr2[i], tolerance) << "for particle pair " << idptr1[i] << "and i=" << i;
    EXPECT_NEAR(zptr1[i], zptr2[i], tolerance) << "for particle pair " << idptr1[i] << "and i=" << i;
    EXPECT_NEAR(fxptr1[i], fxptr2[i], tolerance) << "for particle pair " << idptr1[i] << "and i=" << i;
    EXPECT_NEAR(fyptr1[i], fyptr2[i], tolerance) << "for particle pair " << idptr1[i] << "and i=" << i;
    EXPECT_NEAR(fzptr1[i], fzptr2[i], tolerance) << "for particle pair " << idptr1[i] << "and i=" << i;
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

bool LJFunctorAVXTest::AoSParticlesEqual(FMCell &cell1, FMCell &cell2) {
  EXPECT_GT(cell1.numParticles(), 0);
  EXPECT_EQ(cell1.numParticles(), cell2.numParticles());

  bool ret = true;
  for (size_t i = 0; i < cell1.numParticles(); ++i) {
    ret &= particleEqual(cell1._particles[i], cell2._particles[i]);
  }

  return ret;
}

void LJFunctorAVXTest::testLJFunctorVSLJFunctorAVXTwoCells(bool newton3, bool doDeleteSomeParticles,
                                                           bool useUnalignedViews) {
  FMCell cell1AVX;
  FMCell cell2AVX;

  size_t numParticles = 7;

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::RandomGenerator::fillWithParticles(
      cell1AVX, defaultParticle, _lowCorner, {_highCorner[0] / 2, _highCorner[1], _highCorner[2]}, numParticles);
  autopasTools::generators::RandomGenerator::fillWithParticles(
      cell2AVX, defaultParticle, {_highCorner[0] / 2, _lowCorner[1], _lowCorner[2]}, _highCorner, numParticles);

  if (doDeleteSomeParticles) {
    for (auto &particle : cell1AVX) {
      if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
    }
    for (auto &particle : cell2AVX) {
      if (particle.getID() == 4) autopas::internal::markParticleAsDeleted(particle);
    }
  }

  // copy cells
  FMCell cell1NoAVX(cell1AVX);
  FMCell cell2NoAVX(cell2AVX);

  constexpr bool shifting = true;
  constexpr bool mixing = false;
  autopas::LJFunctor<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorNoAVX(_cutoff);
  ljFunctorNoAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  autopas::LJFunctorAVX<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorAVX(_cutoff);
  ljFunctorAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

  ljFunctorAVX.initTraversal();
  ljFunctorNoAVX.initTraversal();

  ASSERT_TRUE(AoSParticlesEqual(cell1AVX, cell1NoAVX)) << "Cells 1 not equal after copy initialization.";
  ASSERT_TRUE(AoSParticlesEqual(cell2AVX, cell2NoAVX)) << "Cells 2 not equal after copy initialization.";

  ljFunctorNoAVX.SoALoader(cell1NoAVX, cell1NoAVX._particleSoABuffer, 0);
  ljFunctorNoAVX.SoALoader(cell2NoAVX, cell2NoAVX._particleSoABuffer, 0);
  ljFunctorAVX.SoALoader(cell1AVX, cell1AVX._particleSoABuffer, 0);
  ljFunctorAVX.SoALoader(cell2AVX, cell2AVX._particleSoABuffer, 0);

  ASSERT_TRUE(SoAParticlesEqual(cell1AVX._particleSoABuffer, cell1NoAVX._particleSoABuffer))
      << "Cells 1 not equal after loading.";
  ASSERT_TRUE(SoAParticlesEqual(cell2AVX._particleSoABuffer, cell2NoAVX._particleSoABuffer))
      << "Cells 2 not equal after loading.";

  if (useUnalignedViews) {
    ljFunctorNoAVX.SoAFunctorPair(cell1NoAVX._particleSoABuffer.constructView(1, cell1NoAVX.numParticles()),
                                  cell2NoAVX._particleSoABuffer.constructView(1, cell2NoAVX.numParticles()), newton3);
    ljFunctorAVX.SoAFunctorPair(cell1AVX._particleSoABuffer.constructView(1, cell1AVX.numParticles()),
                                cell2AVX._particleSoABuffer.constructView(1, cell2AVX.numParticles()), newton3);
  } else {
    ljFunctorNoAVX.SoAFunctorPair(cell1NoAVX._particleSoABuffer, cell2NoAVX._particleSoABuffer, newton3);
    ljFunctorAVX.SoAFunctorPair(cell1AVX._particleSoABuffer, cell2AVX._particleSoABuffer, newton3);
  }
  ASSERT_TRUE(SoAParticlesEqual(cell1AVX._particleSoABuffer, cell1NoAVX._particleSoABuffer))
      << "Cells 1 not equal after applying functor.";
  ASSERT_TRUE(SoAParticlesEqual(cell2AVX._particleSoABuffer, cell2NoAVX._particleSoABuffer))
      << "Cells 2 not equal after applying functor.";

  ljFunctorAVX.SoAExtractor(cell1AVX, cell1AVX._particleSoABuffer, 0);
  ljFunctorAVX.SoAExtractor(cell2AVX, cell2AVX._particleSoABuffer, 0);
  ljFunctorAVX.SoAExtractor(cell1NoAVX, cell1NoAVX._particleSoABuffer, 0);
  ljFunctorAVX.SoAExtractor(cell2NoAVX, cell2NoAVX._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cell1AVX, cell1NoAVX)) << "Cells 1 not equal after extracting.";
  ASSERT_TRUE(AoSParticlesEqual(cell2AVX, cell2NoAVX)) << "Cells 2 not equal after extracting.";

  ljFunctorAVX.endTraversal(newton3);
  ljFunctorNoAVX.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorAVX.getUpot(), ljFunctorNoAVX.getUpot(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorAVX.getVirial(), ljFunctorNoAVX.getVirial(), tolerance) << "global virial";
}

void LJFunctorAVXTest::testLJFunctorVSLJFunctorAVXOneCell(bool newton3, bool doDeleteSomeParticles,
                                                          bool useUnalignedViews) {
  FMCell cellAVX;

  size_t numParticles = 7;

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::RandomGenerator::fillWithParticles(cellAVX, defaultParticle, _lowCorner, _highCorner,
                                                               numParticles);

  if (doDeleteSomeParticles) {
    for (auto &particle : cellAVX) {
      if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
    }
  }

  // copy cells
  FMCell cellNoAVX(cellAVX);
  constexpr bool shifting = true;
  constexpr bool mixing = false;
  autopas::LJFunctor<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorNoAVX(_cutoff);
  ljFunctorNoAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  autopas::LJFunctorAVX<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorAVX(_cutoff);
  ljFunctorAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

  ASSERT_TRUE(AoSParticlesEqual(cellAVX, cellNoAVX)) << "Cells not equal after copy initialization.";

  ljFunctorAVX.initTraversal();
  ljFunctorNoAVX.initTraversal();

  ljFunctorNoAVX.SoALoader(cellNoAVX, cellNoAVX._particleSoABuffer, 0);
  ljFunctorAVX.SoALoader(cellAVX, cellAVX._particleSoABuffer, 0);

  ASSERT_TRUE(SoAParticlesEqual(cellAVX._particleSoABuffer, cellNoAVX._particleSoABuffer))
      << "Cells not equal after loading.";

  if (useUnalignedViews) {
    ljFunctorNoAVX.SoAFunctorSingle(cellNoAVX._particleSoABuffer.constructView(1, cellNoAVX.numParticles()), newton3);
    ljFunctorAVX.SoAFunctorSingle(cellAVX._particleSoABuffer.constructView(1, cellAVX.numParticles()), newton3);
  } else {
    ljFunctorNoAVX.SoAFunctorSingle(cellNoAVX._particleSoABuffer, newton3);
    ljFunctorAVX.SoAFunctorSingle(cellAVX._particleSoABuffer, newton3);
  }
  ASSERT_TRUE(SoAParticlesEqual(cellAVX._particleSoABuffer, cellNoAVX._particleSoABuffer))
      << "Cells not equal after applying functor.";

  ljFunctorAVX.SoAExtractor(cellAVX, cellAVX._particleSoABuffer, 0);
  ljFunctorAVX.SoAExtractor(cellNoAVX, cellNoAVX._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cellAVX, cellNoAVX)) << "Cells 1 not equal after extracting.";

  ljFunctorAVX.endTraversal(newton3);
  ljFunctorNoAVX.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorAVX.getUpot(), ljFunctorNoAVX.getUpot(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorAVX.getVirial(), ljFunctorNoAVX.getVirial(), tolerance) << "global virial";
}

void LJFunctorAVXTest::testLJFunctorVSLJFunctorAVXVerlet(bool newton3, bool doDeleteSomeParticles) {
  FMCell cellAVX;

  constexpr size_t numParticles = 7;

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::RandomGenerator::fillWithParticles(cellAVX, defaultParticle, _lowCorner, _highCorner,
                                                               numParticles);

  if (doDeleteSomeParticles) {
    // mark some particles as deleted to test if the functor handles them correctly
    for (auto &particle : cellAVX) {
      if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
    }
  }

  // generate neighbor lists
  std::array<std::vector<size_t, autopas::AlignedAllocator<size_t>>, numParticles> neighborLists;
  for (size_t i = 0; i < numParticles; ++i) {
    for (size_t j = newton3 ? i + 1 : 0; j < numParticles; ++j) {
      if (i == j) {
        continue;
      }
      auto dr = autopas::utils::ArrayMath::sub(cellAVX[i].getR(), cellAVX[j].getR());
      double dr2 = autopas::utils::ArrayMath::dot(dr, dr);
      if (dr2 <= _interactionLengthSquare) {
        neighborLists[i].push_back(j);
      }
    }
  }

  // copy cells
  FMCell cellNoAVX(cellAVX);
  constexpr bool shifting = true;
  constexpr bool mixing = false;
  constexpr bool calculateGlobals = true;
  autopas::LJFunctor<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, calculateGlobals> ljFunctorNoAVX(
      _cutoff);
  ljFunctorNoAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  autopas::LJFunctorAVX<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, calculateGlobals> ljFunctorAVX(
      _cutoff);
  ljFunctorAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

  ASSERT_TRUE(AoSParticlesEqual(cellAVX, cellNoAVX)) << "Cells not equal after copy initialization.";

  ljFunctorAVX.initTraversal();
  ljFunctorNoAVX.initTraversal();

  ljFunctorNoAVX.SoALoader(cellNoAVX, cellNoAVX._particleSoABuffer, 0);
  ljFunctorAVX.SoALoader(cellAVX, cellAVX._particleSoABuffer, 0);

  ASSERT_TRUE(SoAParticlesEqual(cellAVX._particleSoABuffer, cellNoAVX._particleSoABuffer))
      << "Cells not equal after loading.";

  for (size_t i = 0; i < numParticles; ++i) {
    ljFunctorNoAVX.SoAFunctorVerlet(cellNoAVX._particleSoABuffer, i, neighborLists[i], newton3);
    ljFunctorAVX.SoAFunctorVerlet(cellAVX._particleSoABuffer, i, neighborLists[i], newton3);
  }

  ASSERT_TRUE(SoAParticlesEqual(cellAVX._particleSoABuffer, cellNoAVX._particleSoABuffer))
      << "Cells not equal after applying functor.";

  ljFunctorAVX.SoAExtractor(cellAVX, cellAVX._particleSoABuffer, 0);
  ljFunctorAVX.SoAExtractor(cellNoAVX, cellNoAVX._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cellAVX, cellNoAVX)) << "Cells not equal after extracting.";

  ljFunctorAVX.endTraversal(newton3);
  ljFunctorNoAVX.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorAVX.getUpot(), ljFunctorNoAVX.getUpot(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorAVX.getVirial(), ljFunctorNoAVX.getVirial(), tolerance) << "global virial";
}

void LJFunctorAVXTest::testLJFunctorVSLJFunctorAVXAoS(bool newton3, bool doDeleteSomeParticles) {
  FMCell cellAVX;

  constexpr size_t numParticles = 7;

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::RandomGenerator::fillWithParticles(cellAVX, defaultParticle, _lowCorner, _highCorner,
                                                               numParticles);

  if (doDeleteSomeParticles) {
    // mark some particles as deleted to test if the functor handles them correctly
    for (auto &particle : cellAVX) {
      if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
    }
  }

  // copy cells
  FMCell cellNoAVX(cellAVX);
  constexpr bool shifting = true;
  constexpr bool mixing = false;
  autopas::LJFunctor<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorNoAVX(_cutoff);
  ljFunctorNoAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  autopas::LJFunctorAVX<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorAVX(_cutoff);
  ljFunctorAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

  ASSERT_TRUE(AoSParticlesEqual(cellAVX, cellNoAVX)) << "Cells not equal after copy initialization.";

  ljFunctorAVX.initTraversal();
  ljFunctorNoAVX.initTraversal();

  for (size_t i = 0; i < numParticles; ++i) {
    for (size_t j = newton3 ? i + 1 : 0; j < numParticles; ++j) {
      if (i == j) {
        continue;
      }
      ljFunctorNoAVX.AoSFunctor(cellNoAVX[i], cellNoAVX[j], newton3);
      ljFunctorAVX.AoSFunctor(cellAVX[i], cellAVX[j], newton3);
    }
  }

  ASSERT_TRUE(AoSParticlesEqual(cellAVX, cellNoAVX)) << "Cells not equal after applying AoSfunctor.";

  ljFunctorAVX.endTraversal(newton3);
  ljFunctorNoAVX.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorAVX.getUpot(), ljFunctorNoAVX.getUpot(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorAVX.getVirial(), ljFunctorNoAVX.getVirial(), tolerance) << "global virial";
}

TEST_P(LJFunctorAVXTest, testLJFunctorVSLJFunctorAVXAoS) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testLJFunctorVSLJFunctorAVXAoS(newton3, doDeleteSomeParticle);
}

TEST_P(LJFunctorAVXTest, testLJFunctorVSLJFunctorAVXVerlet) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testLJFunctorVSLJFunctorAVXVerlet(newton3, doDeleteSomeParticle);
}

TEST_P(LJFunctorAVXTest, testLJFunctorVSLJFunctorAVXOneCellAlignedAccess) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testLJFunctorVSLJFunctorAVXOneCell(newton3, doDeleteSomeParticle, false);
}

TEST_P(LJFunctorAVXTest, testLJFunctorVSLJFunctorAVXOneCellUseUnalignedViews) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testLJFunctorVSLJFunctorAVXOneCell(newton3, doDeleteSomeParticle, true);
}

TEST_P(LJFunctorAVXTest, testLJFunctorVSLJFunctorAVXTwoCellsAlignedAccess) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testLJFunctorVSLJFunctorAVXTwoCells(newton3, doDeleteSomeParticle, false);
}

TEST_P(LJFunctorAVXTest, testLJFunctorVSLJFunctorAVXTwoCellsUseUnalignedViews) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testLJFunctorVSLJFunctorAVXTwoCells(newton3, doDeleteSomeParticle, true);
}

/**
 * Lambda to generate a readable string out of the parameters of this test.
 */
static auto toString = [](const auto &info) {
  auto [newton3, doDeleteSomeParticle] = info.param;
  std::stringstream resStream;
  resStream << (newton3 ? "N3" : "noN3") << "_" << (doDeleteSomeParticle ? "withDeletions" : "noDeletions");
  std::string res = resStream.str();
  std::replace(res.begin(), res.end(), '-', '_');
  std::replace(res.begin(), res.end(), '.', '_');
  return res;
};

INSTANTIATE_TEST_SUITE_P(Generated, LJFunctorAVXTest, ::testing::Combine(::testing::Bool(), ::testing::Bool()),
                         toString);

#endif  // __AVX__

}  // end namespace LJFunctorAVXTest
