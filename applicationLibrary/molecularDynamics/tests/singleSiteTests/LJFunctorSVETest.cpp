/**
 * @file LJFunctorSVETest.cpp
 * @author T. Eke
 * @date 15.11.2021
 */

#ifdef __ARM_FEATURE_SVE

#include "LJFunctorSVETest.h"

#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/Particle.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "molecularDynamicsLibrary/LJFunctorSVE.h"

template <class SoAType>
bool LJFunctorSVETest::SoAParticlesEqual(autopas::SoA<SoAType> &soa1, autopas::SoA<SoAType> &soa2) {
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

bool LJFunctorSVETest::particleEqual(Particle &p1, Particle &p2) {
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

bool LJFunctorSVETest::AoSParticlesEqual(FMCell &cell1, FMCell &cell2) {
  EXPECT_GT(cell1.numParticles(), 0);
  EXPECT_EQ(cell1.numParticles(), cell2.numParticles());

  bool ret = true;
  for (size_t i = 0; i < cell1.numParticles(); ++i) {
    ret = ret and particleEqual(cell1._particles[i], cell2._particles[i]);
  }

  return ret;
}

void LJFunctorSVETest::testLJFunctorVSLJFunctorSVETwoCells(bool newton3, bool doDeleteSomeParticles,
                                                           bool useUnalignedViews) {
  FMCell cell1SVE;
  FMCell cell2SVE;

  size_t numParticles = 7;

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::RandomGenerator::fillWithParticles(
      cell1SVE, defaultParticle, _lowCorner, {_highCorner[0] / 2, _highCorner[1], _highCorner[2]}, numParticles);
  autopasTools::generators::RandomGenerator::fillWithParticles(
      cell2SVE, defaultParticle, {_highCorner[0] / 2, _lowCorner[1], _lowCorner[2]}, _highCorner, numParticles);

  if (doDeleteSomeParticles) {
    for (auto &particle : cell1SVE) {
      if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
    }
    for (auto &particle : cell2SVE) {
      if (particle.getID() == 4) autopas::internal::markParticleAsDeleted(particle);
    }
  }

  // copy cells
  FMCell cell1NoSVE(cell1SVE);
  FMCell cell2NoSVE(cell2SVE);

  constexpr bool shifting = true;
  constexpr bool mixing = false;
  mdLib::LJFunctor<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorNoSVE(_cutoff);
  ljFunctorNoSVE.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  mdLib::LJFunctorSVE<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorSVE(_cutoff);
  ljFunctorSVE.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

  ljFunctorSVE.initTraversal();
  ljFunctorNoSVE.initTraversal();

  ASSERT_TRUE(AoSParticlesEqual(cell1SVE, cell1NoSVE)) << "Cells 1 not equal after copy initialization.";
  ASSERT_TRUE(AoSParticlesEqual(cell2SVE, cell2NoSVE)) << "Cells 2 not equal after copy initialization.";

  ljFunctorNoSVE.SoALoader(cell1NoSVE, cell1NoSVE._particleSoABuffer, 0);
  ljFunctorNoSVE.SoALoader(cell2NoSVE, cell2NoSVE._particleSoABuffer, 0);
  ljFunctorSVE.SoALoader(cell1SVE, cell1SVE._particleSoABuffer, 0);
  ljFunctorSVE.SoALoader(cell2SVE, cell2SVE._particleSoABuffer, 0);

  ASSERT_TRUE(SoAParticlesEqual(cell1SVE._particleSoABuffer, cell1NoSVE._particleSoABuffer))
      << "Cells 1 not equal after loading.";
  ASSERT_TRUE(SoAParticlesEqual(cell2SVE._particleSoABuffer, cell2NoSVE._particleSoABuffer))
      << "Cells 2 not equal after loading.";

  if (useUnalignedViews) {
    ljFunctorNoSVE.SoAFunctorPair(cell1NoSVE._particleSoABuffer.constructView(1, cell1NoSVE.numParticles()),
                                  cell2NoSVE._particleSoABuffer.constructView(1, cell2NoSVE.numParticles()), newton3);
    ljFunctorSVE.SoAFunctorPair(cell1SVE._particleSoABuffer.constructView(1, cell1SVE.numParticles()),
                                cell2SVE._particleSoABuffer.constructView(1, cell2SVE.numParticles()), newton3);
  } else {
    ljFunctorNoSVE.SoAFunctorPair(cell1NoSVE._particleSoABuffer, cell2NoSVE._particleSoABuffer, newton3);
    ljFunctorSVE.SoAFunctorPair(cell1SVE._particleSoABuffer, cell2SVE._particleSoABuffer, newton3);
  }
  ASSERT_TRUE(SoAParticlesEqual(cell1SVE._particleSoABuffer, cell1NoSVE._particleSoABuffer))
      << "Cells 1 not equal after applying functor.";
  ASSERT_TRUE(SoAParticlesEqual(cell2SVE._particleSoABuffer, cell2NoSVE._particleSoABuffer))
      << "Cells 2 not equal after applying functor.";

  ljFunctorSVE.SoAExtractor(cell1SVE, cell1SVE._particleSoABuffer, 0);
  ljFunctorSVE.SoAExtractor(cell2SVE, cell2SVE._particleSoABuffer, 0);
  ljFunctorSVE.SoAExtractor(cell1NoSVE, cell1NoSVE._particleSoABuffer, 0);
  ljFunctorSVE.SoAExtractor(cell2NoSVE, cell2NoSVE._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cell1SVE, cell1NoSVE)) << "Cells 1 not equal after extracting.";
  ASSERT_TRUE(AoSParticlesEqual(cell2SVE, cell2NoSVE)) << "Cells 2 not equal after extracting.";

  ljFunctorSVE.endTraversal(newton3);
  ljFunctorNoSVE.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorSVE.getPotentialEnergy(), ljFunctorNoSVE.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorSVE.getVirial(), ljFunctorNoSVE.getVirial(), tolerance) << "global virial";
}

void LJFunctorSVETest::testLJFunctorVSLJFunctorSVEOneCell(bool newton3, bool doDeleteSomeParticles,
                                                          bool useUnalignedViews) {
  FMCell cellSVE;

  size_t numParticles = 7;

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::RandomGenerator::fillWithParticles(cellSVE, defaultParticle, _lowCorner, _highCorner,
                                                               numParticles);

  if (doDeleteSomeParticles) {
    for (auto &particle : cellSVE) {
      if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
    }
  }

  // copy cells
  FMCell cellNoSVE(cellSVE);
  constexpr bool shifting = true;
  constexpr bool mixing = false;
  mdLib::LJFunctor<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorNoSVE(_cutoff);
  ljFunctorNoSVE.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  mdLib::LJFunctorSVE<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorSVE(_cutoff);
  ljFunctorSVE.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

  ASSERT_TRUE(AoSParticlesEqual(cellSVE, cellNoSVE)) << "Cells not equal after copy initialization.";

  ljFunctorSVE.initTraversal();
  ljFunctorNoSVE.initTraversal();

  ljFunctorNoSVE.SoALoader(cellNoSVE, cellNoSVE._particleSoABuffer, 0);
  ljFunctorSVE.SoALoader(cellSVE, cellSVE._particleSoABuffer, 0);

  ASSERT_TRUE(SoAParticlesEqual(cellSVE._particleSoABuffer, cellNoSVE._particleSoABuffer))
      << "Cells not equal after loading.";

  if (useUnalignedViews) {
    ljFunctorNoSVE.SoAFunctorSingle(cellNoSVE._particleSoABuffer.constructView(1, cellNoSVE.numParticles()), newton3);
    ljFunctorSVE.SoAFunctorSingle(cellSVE._particleSoABuffer.constructView(1, cellSVE.numParticles()), newton3);
  } else {
    ljFunctorNoSVE.SoAFunctorSingle(cellNoSVE._particleSoABuffer, newton3);
    ljFunctorSVE.SoAFunctorSingle(cellSVE._particleSoABuffer, newton3);
  }
  ASSERT_TRUE(SoAParticlesEqual(cellSVE._particleSoABuffer, cellNoSVE._particleSoABuffer))
      << "Cells not equal after applying functor.";

  ljFunctorSVE.SoAExtractor(cellSVE, cellSVE._particleSoABuffer, 0);
  ljFunctorSVE.SoAExtractor(cellNoSVE, cellNoSVE._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cellSVE, cellNoSVE)) << "Cells 1 not equal after extracting.";

  ljFunctorSVE.endTraversal(newton3);
  ljFunctorNoSVE.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorSVE.getPotentialEnergy(), ljFunctorNoSVE.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorSVE.getVirial(), ljFunctorNoSVE.getVirial(), tolerance) << "global virial";
}

void LJFunctorSVETest::testLJFunctorVSLJFunctorSVEVerlet(bool newton3, bool doDeleteSomeParticles) {
  using namespace autopas::utils::ArrayMath::literals;

  FMCell cellSVE;

  constexpr size_t numParticles = 7;

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::RandomGenerator::fillWithParticles(cellSVE, defaultParticle, _lowCorner, _highCorner,
                                                               numParticles);

  if (doDeleteSomeParticles) {
    // mark some particles as deleted to test if the functor handles them correctly
    for (auto &particle : cellSVE) {
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
      auto dr = cellSVE[i].getR() - cellSVE[j].getR();
      double dr2 = autopas::utils::ArrayMath::dot(dr, dr);
      if (dr2 <= _interactionLengthSquare) {
        neighborLists[i].push_back(j);
      }
    }
  }

  // copy cells
  FMCell cellNoSVE(cellSVE);
  constexpr bool shifting = true;
  constexpr bool mixing = false;
  constexpr bool calculateGlobals = true;
  mdLib::LJFunctor<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, calculateGlobals> ljFunctorNoSVE(_cutoff);
  ljFunctorNoSVE.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  mdLib::LJFunctorSVE<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, calculateGlobals> ljFunctorSVE(
      _cutoff);
  ljFunctorSVE.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

  ASSERT_TRUE(AoSParticlesEqual(cellSVE, cellNoSVE)) << "Cells not equal after copy initialization.";

  ljFunctorSVE.initTraversal();
  ljFunctorNoSVE.initTraversal();

  ljFunctorNoSVE.SoALoader(cellNoSVE, cellNoSVE._particleSoABuffer, 0);
  ljFunctorSVE.SoALoader(cellSVE, cellSVE._particleSoABuffer, 0);

  ASSERT_TRUE(SoAParticlesEqual(cellSVE._particleSoABuffer, cellNoSVE._particleSoABuffer))
      << "Cells not equal after loading.";

  for (size_t i = 0; i < numParticles; ++i) {
    ljFunctorNoSVE.SoAFunctorVerlet(cellNoSVE._particleSoABuffer, i, neighborLists[i], newton3);
    ljFunctorSVE.SoAFunctorVerlet(cellSVE._particleSoABuffer, i, neighborLists[i], newton3);
  }

  ASSERT_TRUE(SoAParticlesEqual(cellSVE._particleSoABuffer, cellNoSVE._particleSoABuffer))
      << "Cells not equal after applying functor.";

  ljFunctorSVE.SoAExtractor(cellSVE, cellSVE._particleSoABuffer, 0);
  ljFunctorSVE.SoAExtractor(cellNoSVE, cellNoSVE._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cellSVE, cellNoSVE)) << "Cells not equal after extracting.";

  ljFunctorSVE.endTraversal(newton3);
  ljFunctorNoSVE.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorSVE.getPotentialEnergy(), ljFunctorNoSVE.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorSVE.getVirial(), ljFunctorNoSVE.getVirial(), tolerance) << "global virial";
}

void LJFunctorSVETest::testLJFunctorVSLJFunctorSVEAoS(bool newton3, bool doDeleteSomeParticles) {
  FMCell cellSVE;

  constexpr size_t numParticles = 7;

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::RandomGenerator::fillWithParticles(cellSVE, defaultParticle, _lowCorner, _highCorner,
                                                               numParticles);

  if (doDeleteSomeParticles) {
    // mark some particles as deleted to test if the functor handles them correctly
    for (auto &particle : cellSVE) {
      if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
    }
  }

  // copy cells
  FMCell cellNoSVE(cellSVE);
  constexpr bool shifting = true;
  constexpr bool mixing = false;
  mdLib::LJFunctor<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorNoSVE(_cutoff);
  ljFunctorNoSVE.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  mdLib::LJFunctorSVE<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorSVE(_cutoff);
  ljFunctorSVE.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

  ASSERT_TRUE(AoSParticlesEqual(cellSVE, cellNoSVE)) << "Cells not equal after copy initialization.";

  ljFunctorSVE.initTraversal();
  ljFunctorNoSVE.initTraversal();

  for (size_t i = 0; i < numParticles; ++i) {
    for (size_t j = newton3 ? i + 1 : 0; j < numParticles; ++j) {
      if (i == j) {
        continue;
      }
      ljFunctorNoSVE.AoSFunctor(cellNoSVE[i], cellNoSVE[j], newton3);
      ljFunctorSVE.AoSFunctor(cellSVE[i], cellSVE[j], newton3);
    }
  }

  ASSERT_TRUE(AoSParticlesEqual(cellSVE, cellNoSVE)) << "Cells not equal after applying AoSfunctor.";

  ljFunctorSVE.endTraversal(newton3);
  ljFunctorNoSVE.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorSVE.getPotentialEnergy(), ljFunctorNoSVE.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorSVE.getVirial(), ljFunctorNoSVE.getVirial(), tolerance) << "global virial";
}
* /

    TEST_P(LJFunctorSVETest, testLJFunctorVSLJFunctorSVEAoS) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testLJFunctorVSLJFunctorSVEAoS(newton3, doDeleteSomeParticle);
}

TEST_P(LJFunctorSVETest, testLJFunctorVSLJFunctorSVEVerlet) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testLJFunctorVSLJFunctorSVEVerlet(newton3, doDeleteSomeParticle);
}

TEST_P(LJFunctorSVETest, testLJFunctorVSLJFunctorSVEOneCellAlignedAccess) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testLJFunctorVSLJFunctorSVEOneCell(newton3, doDeleteSomeParticle, false);
}

TEST_P(LJFunctorSVETest, testLJFunctorVSLJFunctorSVEOneCellUseUnalignedViews) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testLJFunctorVSLJFunctorSVEOneCell(newton3, doDeleteSomeParticle, true);
}

TEST_P(LJFunctorSVETest, testLJFunctorVSLJFunctorSVETwoCellsAlignedAccess) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testLJFunctorVSLJFunctorSVETwoCells(newton3, doDeleteSomeParticle, false);
}

TEST_P(LJFunctorSVETest, testLJFunctorVSLJFunctorSVETwoCellsUseUnalignedViews) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testLJFunctorVSLJFunctorSVETwoCells(newton3, doDeleteSomeParticle, true);
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

INSTANTIATE_TEST_SUITE_P(Generated, LJFunctorSVETest, ::testing::Combine(::testing::Bool(), ::testing::Bool()),
                         toString);

#endif  // __ARM_FEATURE_SVE
