

#include "LJFunctorMIPPTest.h"

#include "autopas/cells/FullParticleCell.h"
#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopas/molecularDynamics/LJFunctorMIPP.h"
#include "autopas/particles/Particle.h"
#include "autopasTools/generators/RandomGenerator.h"

template <class SoAType>
bool LJFunctorMIPPTest::SoAParticlesEqual(autopas::SoA<SoAType> &soa1, autopas::SoA<SoAType> &soa2) {
  EXPECT_GT(soa1.getNumberOfParticles(), 0);
  EXPECT_EQ(soa1.getNumberOfParticles(), soa2.getNumberOfParticles());

  unsigned long *const __restrict idptr1 = soa1.template begin<Particle::AttributeNames::id>();
  unsigned long *const __restrict idptr2 = soa2.template begin<Particle::AttributeNames::id>();

  double *const __restrict xptr1 = soa1.template begin<Particle::AttributeNames::posX>();
  double *const __restrict yptr1 = soa1.template begin<Particle::AttributeNames::posY>();
  double *const __restrict zptr1 = soa1.template begin<Particle::AttributeNames::posZ>();
  double *const __restrict xptr2 = soa2.template begin<Particle::AttributeNames::posX>();
  double *const __restrict yptr2 = soa2.template begin<Particle::AttributeNames::posY>();
  double *const __restrict zptr2 = soa2.template begin<Particle::AttributeNames::posZ>();

  double *const __restrict fxptr1 = soa1.template begin<Particle::AttributeNames::forceX>();
  double *const __restrict fyptr1 = soa1.template begin<Particle::AttributeNames::forceY>();
  double *const __restrict fzptr1 = soa1.template begin<Particle::AttributeNames::forceZ>();
  double *const __restrict fxptr2 = soa2.template begin<Particle::AttributeNames::forceX>();
  double *const __restrict fyptr2 = soa2.template begin<Particle::AttributeNames::forceY>();
  double *const __restrict fzptr2 = soa2.template begin<Particle::AttributeNames::forceZ>();

  for (size_t i = 0; i < soa1.getNumberOfParticles(); ++i) {
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

bool LJFunctorMIPPTest::particleEqual(Particle &p1, Particle &p2) {
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

bool LJFunctorMIPPTest::AoSParticlesEqual(FMCell &cell1, FMCell &cell2) {
  EXPECT_GT(cell1.numParticles(), 0);
  EXPECT_EQ(cell1.numParticles(), cell2.numParticles());

  bool ret = true;
  for (size_t i = 0; i < cell1.numParticles(); ++i) {
    ret = ret and particleEqual(cell1._particles[i], cell2._particles[i]);
  }

  return ret;
}

void LJFunctorMIPPTest::testLJFunctorVSLJFunctorMIPPTwoCells(bool newton3, bool doDeleteSomeParticles,
                                                           bool useUnalignedViews) {
  FMCell cell1MIPP;
  FMCell cell2MIPP;

  size_t numParticles = 7;

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::RandomGenerator::fillWithParticles(
      cell1MIPP, defaultParticle, _lowCorner, {_highCorner[0] / 2, _highCorner[1], _highCorner[2]}, numParticles);
  autopasTools::generators::RandomGenerator::fillWithParticles(
      cell2MIPP, defaultParticle, {_highCorner[0] / 2, _lowCorner[1], _lowCorner[2]}, _highCorner, numParticles);

  if (doDeleteSomeParticles) {
    for (auto &particle : cell1MIPP) {
      if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
    }
    for (auto &particle : cell2MIPP) {
      if (particle.getID() == 4) autopas::internal::markParticleAsDeleted(particle);
    }
  }

  // copy cells
  FMCell cell1NoMIPP(cell1MIPP);
  FMCell cell2NoMIPP(cell2MIPP);

  constexpr bool shifting = true;
  constexpr bool mixing = false;
  autopas::LJFunctor<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorNoMIPP(_cutoff);
  ljFunctorNoMIPP.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  autopas::LJFunctorMIPP<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorMIPP(_cutoff);
  ljFunctorMIPP.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

  ljFunctorMIPP.initTraversal();
  ljFunctorNoMIPP.initTraversal();

  ASSERT_TRUE(AoSParticlesEqual(cell1MIPP, cell1NoMIPP)) << "Cells 1 not equal after copy initialization.";
  ASSERT_TRUE(AoSParticlesEqual(cell2MIPP, cell2NoMIPP)) << "Cells 2 not equal after copy initialization.";

  ljFunctorNoMIPP.SoALoader(cell1NoMIPP, cell1NoMIPP._particleSoABuffer, 0);
  ljFunctorNoMIPP.SoALoader(cell2NoMIPP, cell2NoMIPP._particleSoABuffer, 0);
  ljFunctorMIPP.SoALoader(cell1MIPP, cell1MIPP._particleSoABuffer, 0);
  ljFunctorMIPP.SoALoader(cell2MIPP, cell2MIPP._particleSoABuffer, 0);

  ASSERT_TRUE(SoAParticlesEqual(cell1MIPP._particleSoABuffer, cell1NoMIPP._particleSoABuffer))
      << "Cells 1 not equal after loading.";
  ASSERT_TRUE(SoAParticlesEqual(cell2MIPP._particleSoABuffer, cell2NoMIPP._particleSoABuffer))
      << "Cells 2 not equal after loading.";

  if (useUnalignedViews) {
    ljFunctorNoMIPP.SoAFunctorPair(cell1NoMIPP._particleSoABuffer.constructView(1, cell1NoMIPP.numParticles()),
                                  cell2NoMIPP._particleSoABuffer.constructView(1, cell2NoMIPP.numParticles()), newton3);
    ljFunctorMIPP.SoAFunctorPair(cell1MIPP._particleSoABuffer.constructView(1, cell1MIPP.numParticles()),
                                cell2MIPP._particleSoABuffer.constructView(1, cell2MIPP.numParticles()), newton3);
  } else {
    ljFunctorNoMIPP.SoAFunctorPair(cell1NoMIPP._particleSoABuffer, cell2NoMIPP._particleSoABuffer, newton3);
    ljFunctorMIPP.SoAFunctorPair(cell1MIPP._particleSoABuffer, cell2MIPP._particleSoABuffer, newton3);
  }
  ASSERT_TRUE(SoAParticlesEqual(cell1MIPP._particleSoABuffer, cell1NoMIPP._particleSoABuffer))
      << "Cells 1 not equal after applying functor.";
  ASSERT_TRUE(SoAParticlesEqual(cell2MIPP._particleSoABuffer, cell2NoMIPP._particleSoABuffer))
      << "Cells 2 not equal after applying functor.";

  ljFunctorMIPP.SoAExtractor(cell1MIPP, cell1MIPP._particleSoABuffer, 0);
  ljFunctorMIPP.SoAExtractor(cell2MIPP, cell2MIPP._particleSoABuffer, 0);
  ljFunctorMIPP.SoAExtractor(cell1NoMIPP, cell1NoMIPP._particleSoABuffer, 0);
  ljFunctorMIPP.SoAExtractor(cell2NoMIPP, cell2NoMIPP._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cell1MIPP, cell1NoMIPP)) << "Cells 1 not equal after extracting.";
  ASSERT_TRUE(AoSParticlesEqual(cell2MIPP, cell2NoMIPP)) << "Cells 2 not equal after extracting.";

  ljFunctorMIPP.endTraversal(newton3);
  ljFunctorNoMIPP.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorMIPP.getUpot(), ljFunctorNoMIPP.getUpot(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorMIPP.getVirial(), ljFunctorNoMIPP.getVirial(), tolerance) << "global virial";
}

void LJFunctorMIPPTest::testLJFunctorVSLJFunctorMIPPOneCell(bool newton3, bool doDeleteSomeParticles,
                                                          bool useUnalignedViews) {
  FMCell cellMIPP;

  size_t numParticles = 7;

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::RandomGenerator::fillWithParticles(cellMIPP, defaultParticle, _lowCorner, _highCorner,
                                                               numParticles);

  if (doDeleteSomeParticles) {
    for (auto &particle : cellMIPP) {
      if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
    }
  }

  // copy cells
  FMCell cellNoMIPP(cellMIPP);
  constexpr bool shifting = true;
  constexpr bool mixing = false;
  autopas::LJFunctor<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorNoMIPP(_cutoff);
  ljFunctorNoMIPP.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  autopas::LJFunctorMIPP<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorMIPP(_cutoff);
  ljFunctorMIPP.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

  ASSERT_TRUE(AoSParticlesEqual(cellMIPP, cellNoMIPP)) << "Cells not equal after copy initialization.";

  ljFunctorMIPP.initTraversal();
  ljFunctorNoMIPP.initTraversal();

  ljFunctorNoMIPP.SoALoader(cellNoMIPP, cellNoMIPP._particleSoABuffer, 0);
  ljFunctorMIPP.SoALoader(cellMIPP, cellMIPP._particleSoABuffer, 0);

  ASSERT_TRUE(SoAParticlesEqual(cellMIPP._particleSoABuffer, cellNoMIPP._particleSoABuffer))
      << "Cells not equal after loading.";

  if (useUnalignedViews) {
    ljFunctorNoMIPP.SoAFunctorSingle(cellNoMIPP._particleSoABuffer.constructView(1, cellNoMIPP.numParticles()), newton3);
    ljFunctorMIPP.SoAFunctorSingle(cellMIPP._particleSoABuffer.constructView(1, cellMIPP.numParticles()), newton3);
  } else {
    ljFunctorNoMIPP.SoAFunctorSingle(cellNoMIPP._particleSoABuffer, newton3);
    ljFunctorMIPP.SoAFunctorSingle(cellMIPP._particleSoABuffer, newton3);
  }
  ASSERT_TRUE(SoAParticlesEqual(cellMIPP._particleSoABuffer, cellNoMIPP._particleSoABuffer))
      << "Cells not equal after applying functor.";

  ljFunctorMIPP.SoAExtractor(cellMIPP, cellMIPP._particleSoABuffer, 0);
  ljFunctorMIPP.SoAExtractor(cellNoMIPP, cellNoMIPP._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cellMIPP, cellNoMIPP)) << "Cells 1 not equal after extracting.";

  ljFunctorMIPP.endTraversal(newton3);
  ljFunctorNoMIPP.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorMIPP.getUpot(), ljFunctorNoMIPP.getUpot(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorMIPP.getVirial(), ljFunctorNoMIPP.getVirial(), tolerance) << "global virial";
}

void LJFunctorMIPPTest::testLJFunctorVSLJFunctorMIPPVerlet(bool newton3, bool doDeleteSomeParticles) {
  using namespace autopas::utils::ArrayMath::literals;

  FMCell cellMIPP;

  constexpr size_t numParticles = 7;

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::RandomGenerator::fillWithParticles(cellMIPP, defaultParticle, _lowCorner, _highCorner,
                                                               numParticles);

  if (doDeleteSomeParticles) {
    // mark some particles as deleted to test if the functor handles them correctly
    for (auto &particle : cellMIPP) {
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
      auto dr = cellMIPP[i].getR() - cellMIPP[j].getR();
      double dr2 = autopas::utils::ArrayMath::dot(dr, dr);
      if (dr2 <= _interactionLengthSquare) {
        neighborLists[i].push_back(j);
      }
    }
  }

  // copy cells
  FMCell cellNoMIPP(cellMIPP);
  constexpr bool shifting = true;
  constexpr bool mixing = false;
  constexpr bool calculateGlobals = true;
  autopas::LJFunctor<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, calculateGlobals> ljFunctorNoMIPP(
      _cutoff);
  ljFunctorNoMIPP.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  autopas::LJFunctorMIPP<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, calculateGlobals> ljFunctorMIPP(
      _cutoff);
  ljFunctorMIPP.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

  ASSERT_TRUE(AoSParticlesEqual(cellMIPP, cellNoMIPP)) << "Cells not equal after copy initialization.";

  ljFunctorMIPP.initTraversal();
  ljFunctorNoMIPP.initTraversal();

  ljFunctorNoMIPP.SoALoader(cellNoMIPP, cellNoMIPP._particleSoABuffer, 0);
  ljFunctorMIPP.SoALoader(cellMIPP, cellMIPP._particleSoABuffer, 0);

  ASSERT_TRUE(SoAParticlesEqual(cellMIPP._particleSoABuffer, cellNoMIPP._particleSoABuffer))
      << "Cells not equal after loading.";

  for (size_t i = 0; i < numParticles; ++i) {
    ljFunctorNoMIPP.SoAFunctorVerlet(cellNoMIPP._particleSoABuffer, i, neighborLists[i], newton3);
    ljFunctorMIPP.SoAFunctorVerlet(cellMIPP._particleSoABuffer, i, neighborLists[i], newton3);
  }

  ASSERT_TRUE(SoAParticlesEqual(cellMIPP._particleSoABuffer, cellNoMIPP._particleSoABuffer))
      << "Cells not equal after applying functor.";

  ljFunctorMIPP.SoAExtractor(cellMIPP, cellMIPP._particleSoABuffer, 0);
  ljFunctorMIPP.SoAExtractor(cellNoMIPP, cellNoMIPP._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cellMIPP, cellNoMIPP)) << "Cells not equal after extracting.";

  ljFunctorMIPP.endTraversal(newton3);
  ljFunctorNoMIPP.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorMIPP.getUpot(), ljFunctorNoMIPP.getUpot(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorMIPP.getVirial(), ljFunctorNoMIPP.getVirial(), tolerance) << "global virial";
}

void LJFunctorMIPPTest::testLJFunctorVSLJFunctorMIPPAoS(bool newton3, bool doDeleteSomeParticles) {
  FMCell cellMIPP;

  constexpr size_t numParticles = 7;

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::RandomGenerator::fillWithParticles(cellMIPP, defaultParticle, _lowCorner, _highCorner,
                                                               numParticles);

  if (doDeleteSomeParticles) {
    // mark some particles as deleted to test if the functor handles them correctly
    for (auto &particle : cellMIPP) {
      if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
    }
  }

  // copy cells
  FMCell cellNoMIPP(cellMIPP);
  constexpr bool shifting = true;
  constexpr bool mixing = false;
  autopas::LJFunctor<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorNoMIPP(_cutoff);
  ljFunctorNoMIPP.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  autopas::LJFunctorMIPP<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorMIPP(_cutoff);
  ljFunctorMIPP.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

  ASSERT_TRUE(AoSParticlesEqual(cellMIPP, cellNoMIPP)) << "Cells not equal after copy initialization.";

  ljFunctorMIPP.initTraversal();
  ljFunctorNoMIPP.initTraversal();

  for (size_t i = 0; i < numParticles; ++i) {
    for (size_t j = newton3 ? i + 1 : 0; j < numParticles; ++j) {
      if (i == j) {
        continue;
      }
      ljFunctorNoMIPP.AoSFunctor(cellNoMIPP[i], cellNoMIPP[j], newton3);
      ljFunctorMIPP.AoSFunctor(cellMIPP[i], cellMIPP[j], newton3);
    }
  }

  ASSERT_TRUE(AoSParticlesEqual(cellMIPP, cellNoMIPP)) << "Cells not equal after applying AoSfunctor.";

  ljFunctorMIPP.endTraversal(newton3);
  ljFunctorNoMIPP.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorMIPP.getUpot(), ljFunctorNoMIPP.getUpot(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorMIPP.getVirial(), ljFunctorNoMIPP.getVirial(), tolerance) << "global virial";
}

TEST_P(LJFunctorMIPPTest, testLJFunctorVSLJFunctorMIPPAoS) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testLJFunctorVSLJFunctorMIPPAoS(newton3, doDeleteSomeParticle);
}

TEST_P(LJFunctorMIPPTest, testLJFunctorVSLJFunctorMIPPVerlet) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testLJFunctorVSLJFunctorMIPPVerlet(newton3, doDeleteSomeParticle);
}

TEST_P(LJFunctorMIPPTest, testLJFunctorVSLJFunctorMIPPOneCellAlignedAccess) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testLJFunctorVSLJFunctorMIPPOneCell(newton3, doDeleteSomeParticle, false);
}

TEST_P(LJFunctorMIPPTest, testLJFunctorVSLJFunctorMIPPOneCellUseUnalignedViews) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testLJFunctorVSLJFunctorMIPPOneCell(newton3, doDeleteSomeParticle, true);
}

TEST_P(LJFunctorMIPPTest, testLJFunctorVSLJFunctorMIPPTwoCellsAlignedAccess) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testLJFunctorVSLJFunctorMIPPTwoCells(newton3, doDeleteSomeParticle, false);
}

TEST_P(LJFunctorMIPPTest, testLJFunctorVSLJFunctorMIPPTwoCellsUseUnalignedViews) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testLJFunctorVSLJFunctorMIPPTwoCells(newton3, doDeleteSomeParticle, true);
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

INSTANTIATE_TEST_SUITE_P(Generated, LJFunctorMIPPTest, ::testing::Combine(::testing::Bool(), ::testing::Bool()),
                         toString);
