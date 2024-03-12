#include "MieFunctorSVETest.h"

#ifdef __ARM_FEATURE_SVE

#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/Particle.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "molecularDynamicsLibrary/LJFunctorSVE.h"
#include "molecularDynamicsLibrary/MieFunctorSVE.h"

template <class SoAType>
bool MieFunctorSVETest::SoAParticlesEqual(autopas::SoA<SoAType> &soa1, autopas::SoA<SoAType> &soa2) {
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

bool MieFunctorSVETest::particleEqual(Particle &p1, Particle &p2) {
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

bool MieFunctorSVETest::AoSParticlesEqual(FMCell &cell1, FMCell &cell2) {
  EXPECT_GT(cell1.numParticles(), 0);
  EXPECT_EQ(cell1.numParticles(), cell2.numParticles());

  bool ret = true;
  for (size_t i = 0; i < cell1.numParticles(); ++i) {
    ret = ret and particleEqual(cell1._particles[i], cell2._particles[i]);
  }

  return ret;
}

void MieFunctorSVETest::testMieFunctorSVEVSLJFunctorSVETwoCells(bool newton3, bool doDeleteSomeParticles,
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
  FMCell cell1MieSVE(cell1SVE);
  FMCell cell2MieSVE(cell2SVE);

  constexpr bool shifting = true;
  constexpr bool mixing = false;
  mdLib::MieFunctorSVE<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> mieFunctorSVE(_cutoff, 12, 6);
  mieFunctorSVE.setParticleProperties(_epsilon, _sigma * _sigma);
  mdLib::LJFunctorSVE<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorSVE(_cutoff);
  ljFunctorSVE.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

  ljFunctorSVE.initTraversal();
  mieFunctorSVE.initTraversal();

  ASSERT_TRUE(AoSParticlesEqual(cell1SVE, cell1MieSVE)) << "Cells 1 not equal after copy initialization.";
  ASSERT_TRUE(AoSParticlesEqual(cell2SVE, cell2MieSVE)) << "Cells 2 not equal after copy initialization.";

  mieFunctorSVE.SoALoader(cell1MieSVE, cell1MieSVE._particleSoABuffer, 0);
  mieFunctorSVE.SoALoader(cell2MieSVE, cell2MieSVE._particleSoABuffer, 0);
  ljFunctorSVE.SoALoader(cell1SVE, cell1SVE._particleSoABuffer, 0);
  ljFunctorSVE.SoALoader(cell2SVE, cell2SVE._particleSoABuffer, 0);

  ASSERT_TRUE(SoAParticlesEqual(cell1SVE._particleSoABuffer, cell1MieSVE._particleSoABuffer))
      << "Cells 1 not equal after loading.";
  ASSERT_TRUE(SoAParticlesEqual(cell2SVE._particleSoABuffer, cell2MieSVE._particleSoABuffer))
      << "Cells 2 not equal after loading.";

  if (useUnalignedViews) {
    mieFunctorSVE.SoAFunctorPair(cell1MieSVE._particleSoABuffer.constructView(1, cell1MieSVE.numParticles()),
                                 cell2MieSVE._particleSoABuffer.constructView(1, cell2MieSVE.numParticles()), newton3);
    ljFunctorSVE.SoAFunctorPair(cell1SVE._particleSoABuffer.constructView(1, cell1SVE.numParticles()),
                                cell2SVE._particleSoABuffer.constructView(1, cell2SVE.numParticles()), newton3);
  } else {
    mieFunctorSVE.SoAFunctorPair(cell1MieSVE._particleSoABuffer, cell2MieSVE._particleSoABuffer, newton3);
    ljFunctorSVE.SoAFunctorPair(cell1SVE._particleSoABuffer, cell2SVE._particleSoABuffer, newton3);
  }
  ASSERT_TRUE(SoAParticlesEqual(cell1SVE._particleSoABuffer, cell1MieSVE._particleSoABuffer))
      << "Cells 1 not equal after applying functor.";
  ASSERT_TRUE(SoAParticlesEqual(cell2SVE._particleSoABuffer, cell2MieSVE._particleSoABuffer))
      << "Cells 2 not equal after applying functor.";

  ljFunctorSVE.SoAExtractor(cell1SVE, cell1SVE._particleSoABuffer, 0);
  ljFunctorSVE.SoAExtractor(cell2SVE, cell2SVE._particleSoABuffer, 0);
  ljFunctorSVE.SoAExtractor(cell1MieSVE, cell1MieSVE._particleSoABuffer, 0);
  ljFunctorSVE.SoAExtractor(cell2MieSVE, cell2MieSVE._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cell1SVE, cell1MieSVE)) << "Cells 1 not equal after extracting.";
  ASSERT_TRUE(AoSParticlesEqual(cell2SVE, cell2MieSVE)) << "Cells 2 not equal after extracting.";

  ljFunctorSVE.endTraversal(newton3);
  mieFunctorSVE.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorSVE.getPotentialEnergy(), mieFunctorSVE.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorSVE.getVirial(), mieFunctorSVE.getVirial(), tolerance) << "global virial";
}

void MieFunctorSVETest::testMieFunctorSVEVSLJFunctorSVEOneCell(bool newton3, bool doDeleteSomeParticles,
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
  FMCell cellMieSVE(cellSVE);
  constexpr bool shifting = true;
  constexpr bool mixing = false;
  mdLib::MieFunctorSVE<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> mieFunctorSVE(_cutoff, 12, 6);
  mieFunctorSVE.setParticleProperties(_epsilon, _sigma * _sigma);
  mdLib::LJFunctorSVE<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorSVE(_cutoff);
  ljFunctorSVE.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

  ASSERT_TRUE(AoSParticlesEqual(cellSVE, cellMieSVE)) << "Cells not equal after copy initialization.";

  ljFunctorSVE.initTraversal();
  mieFunctorSVE.initTraversal();

  mieFunctorSVE.SoALoader(cellMieSVE, cellMieSVE._particleSoABuffer, 0);
  ljFunctorSVE.SoALoader(cellSVE, cellSVE._particleSoABuffer, 0);

  ASSERT_TRUE(SoAParticlesEqual(cellSVE._particleSoABuffer, cellMieSVE._particleSoABuffer))
      << "Cells not equal after loading.";

  if (useUnalignedViews) {
    mieFunctorSVE.SoAFunctorSingle(cellMieSVE._particleSoABuffer.constructView(1, cellMieSVE.numParticles()), newton3);
    ljFunctorSVE.SoAFunctorSingle(cellSVE._particleSoABuffer.constructView(1, cellSVE.numParticles()), newton3);
  } else {
    mieFunctorSVE.SoAFunctorSingle(cellMieSVE._particleSoABuffer, newton3);
    ljFunctorSVE.SoAFunctorSingle(cellSVE._particleSoABuffer, newton3);
  }
  ASSERT_TRUE(SoAParticlesEqual(cellSVE._particleSoABuffer, cellMieSVE._particleSoABuffer))
      << "Cells not equal after applying functor.";

  ljFunctorSVE.SoAExtractor(cellSVE, cellSVE._particleSoABuffer, 0);
  ljFunctorSVE.SoAExtractor(cellMieSVE, cellMieSVE._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cellSVE, cellMieSVE)) << "Cells 1 not equal after extracting.";

  ljFunctorSVE.endTraversal(newton3);
  mieFunctorSVE.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorSVE.getPotentialEnergy(), mieFunctorSVE.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorSVE.getVirial(), mieFunctorSVE.getVirial(), tolerance) << "global virial";
}

void MieFunctorSVETest::testMieFunctorSVEVSLJFunctorSVEVerlet(bool newton3, bool doDeleteSomeParticles) {
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
  FMCell cellMieSVE(cellSVE);
  constexpr bool shifting = true;
  constexpr bool mixing = false;
  constexpr bool calculateGlobals = true;
  mdLib::MieFunctorSVE<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, calculateGlobals> mieFunctorSVE(
      _cutoff, 12, 6);
  mieFunctorSVE.setParticleProperties(_epsilon, _sigma * _sigma);
  mdLib::LJFunctorSVE<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, calculateGlobals> ljFunctorSVE(
      _cutoff);
  ljFunctorSVE.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

  ASSERT_TRUE(AoSParticlesEqual(cellSVE, cellMieSVE)) << "Cells not equal after copy initialization.";

  ljFunctorSVE.initTraversal();
  mieFunctorSVE.initTraversal();

  mieFunctorSVE.SoALoader(cellMieSVE, cellMieSVE._particleSoABuffer, 0);
  ljFunctorSVE.SoALoader(cellSVE, cellSVE._particleSoABuffer, 0);

  ASSERT_TRUE(SoAParticlesEqual(cellSVE._particleSoABuffer, cellMieSVE._particleSoABuffer))
      << "Cells not equal after loading.";

  for (size_t i = 0; i < numParticles; ++i) {
    mieFunctorSVE.SoAFunctorVerlet(cellMieSVE._particleSoABuffer, i, neighborLists[i], newton3);
    ljFunctorSVE.SoAFunctorVerlet(cellSVE._particleSoABuffer, i, neighborLists[i], newton3);
  }
  // Note that the calculated forces wont be compared, but the old forces. The knew forces are in t
  ASSERT_TRUE(SoAParticlesEqual(cellSVE._particleSoABuffer, cellMieSVE._particleSoABuffer))
      << "Cells not equal after applying functor.";

  ljFunctorSVE.SoAExtractor(cellSVE, cellSVE._particleSoABuffer, 0);
  mieFunctorSVE.SoAExtractor(cellMieSVE, cellMieSVE._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cellSVE, cellMieSVE)) << "Cells not equal after extracting.";

  ljFunctorSVE.endTraversal(newton3);
  mieFunctorSVE.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorSVE.getPotentialEnergy(), mieFunctorSVE.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorSVE.getVirial(), mieFunctorSVE.getVirial(), tolerance) << "global virial";
}

void MieFunctorSVETest::testMieFunctorSVEAoS(bool newton3, bool doDeleteSomeParticles) {
  FMCell cellSVE;

  constexpr size_t numParticles = 7;

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::RandomGenerator::fillWithParticles(cellSVE, defaultParticle, _lowCorner, _highCorner,
                                                               numParticles);

  if (doDeleteSomeParticles) {
    // mark some particles as deleted to teqt if the functor handles them correctly
    for (auto &particle : cellSVE) {
      if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
    }
  }

  // copy cells
  FMCell cellMieSVE(cellSVE);
  constexpr bool shifting = true;
  constexpr bool mixing = false;
  mdLib::MieFunctorSVE<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> mieFunctorSVE(_cutoff, 12, 6);
  mieFunctorSVE.setParticleProperties(_epsilon, _sigma * _sigma);
  mdLib::LJFunctorSVE<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorSVE(_cutoff);
  ljFunctorSVE.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

  ASSERT_TRUE(AoSParticlesEqual(cellSVE, cellMieSVE)) << "Cells not equal after copy initialization.";

  ljFunctorSVE.initTraversal();
  mieFunctorSVE.initTraversal();

  for (size_t i = 0; i < numParticles; ++i) {
    for (size_t j = newton3 ? i + 1 : 0; j < numParticles; ++j) {
      if (i == j) {
        continue;
      }
      mieFunctorSVE.AoSFunctor(cellMieSVE[i], cellMieSVE[j], newton3);
      ljFunctorSVE.AoSFunctor(cellSVE[i], cellSVE[j], newton3);
    }
  }

  ASSERT_TRUE(AoSParticlesEqual(cellSVE, cellMieSVE)) << "Cells not equal after applying AoSfunctor.";

  ljFunctorSVE.endTraversal(newton3);
  mieFunctorSVE.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorSVE.getPotentialEnergy(), mieFunctorSVE.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorSVE.getVirial(), mieFunctorSVE.getVirial(), tolerance) << "global virial";
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

TEST_P(MieFunctorSVETest, testMieFunctorSVEVSLJFunctorSVETwoCellsUseUnalignedViews) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testMieFunctorSVEAoS(newton3, doDeleteSomeParticle);
}

TEST_P(MieFunctorSVETest, testMieFunctorSVEVSLJFunctorSVEVerlet) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testMieFunctorSVEVSLJFunctorSVEVerlet(newton3, doDeleteSomeParticle);
}

TEST_P(MieFunctorSVETest, testMieFunctorSVEVSLJFunctorSVEOneCellAlignedAccess) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testMieFunctorSVEVSLJFunctorSVEOneCell(newton3, doDeleteSomeParticle, false);
}

TEST_P(MieFunctorSVETest, testMieFunctorSVEVSLJFunctorSVEOneCellUseUnalignedViews) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testMieFunctorSVEVSLJFunctorSVEOneCell(newton3, doDeleteSomeParticle, true);
}

TEST_P(MieFunctorSVETest, testMieFunctorSVEVSLJFunctorSVETwoCellsAlignedAccess) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testMieFunctorSVEVSLJFunctorSVETwoCells(newton3, doDeleteSomeParticle, false);
}

INSTANTIATE_TEST_SUITE_P(Generated, MieFunctorSVETest, ::testing::Combine(::testing::Bool(), ::testing::Bool()),
                         toString);
#endif
