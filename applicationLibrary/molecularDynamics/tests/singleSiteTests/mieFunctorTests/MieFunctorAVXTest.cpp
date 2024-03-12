#include "MieFunctorAVXTest.h"
/**
 * @file MieFunctorAVXTest.cpp
 * @author K. Cole
 * @date 01/10/23
 */

#ifdef __AVX__

#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/Particle.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "molecularDynamicsLibrary/LJFunctorAVX.h"
#include "molecularDynamicsLibrary/MieFunctorAVX.h"

template <class SoAType>
bool MieFunctorAVXTest::SoAParticlesEqual(autopas::SoA<SoAType> &soa1, autopas::SoA<SoAType> &soa2) {
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

bool MieFunctorAVXTest::particleEqual(Particle &p1, Particle &p2) {
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

bool MieFunctorAVXTest::AoSParticlesEqual(FMCell &cell1, FMCell &cell2) {
  EXPECT_GT(cell1.numParticles(), 0);
  EXPECT_EQ(cell1.numParticles(), cell2.numParticles());

  bool ret = true;
  for (size_t i = 0; i < cell1.numParticles(); ++i) {
    ret = ret and particleEqual(cell1._particles[i], cell2._particles[i]);
  }

  return ret;
}

void MieFunctorAVXTest::testMieFunctorAVXVSLJFunctorAVXTwoCells(bool newton3, bool doDeleteSomeParticles,
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
  FMCell cell1MieAVX(cell1AVX);
  FMCell cell2MieAVX(cell2AVX);

  constexpr bool shifting = true;
  constexpr bool mixing = false;
  mdLib::MieFunctorAVX<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> mieFunctorAVX(_cutoff, 12, 6);
  mieFunctorAVX.setParticleProperties(_epsilon, _sigma * _sigma);
  mdLib::LJFunctorAVX<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorAVX(_cutoff);
  ljFunctorAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

  ljFunctorAVX.initTraversal();
  mieFunctorAVX.initTraversal();

  ASSERT_TRUE(AoSParticlesEqual(cell1AVX, cell1MieAVX)) << "Cells 1 not equal after copy initialization.";
  ASSERT_TRUE(AoSParticlesEqual(cell2AVX, cell2MieAVX)) << "Cells 2 not equal after copy initialization.";

  mieFunctorAVX.SoALoader(cell1MieAVX, cell1MieAVX._particleSoABuffer, 0);
  mieFunctorAVX.SoALoader(cell2MieAVX, cell2MieAVX._particleSoABuffer, 0);
  ljFunctorAVX.SoALoader(cell1AVX, cell1AVX._particleSoABuffer, 0);
  ljFunctorAVX.SoALoader(cell2AVX, cell2AVX._particleSoABuffer, 0);

  ASSERT_TRUE(SoAParticlesEqual(cell1AVX._particleSoABuffer, cell1MieAVX._particleSoABuffer))
      << "Cells 1 not equal after loading.";
  ASSERT_TRUE(SoAParticlesEqual(cell2AVX._particleSoABuffer, cell2MieAVX._particleSoABuffer))
      << "Cells 2 not equal after loading.";

  if (useUnalignedViews) {
    mieFunctorAVX.SoAFunctorPair(cell1MieAVX._particleSoABuffer.constructView(1, cell1MieAVX.numParticles()),
                                 cell2MieAVX._particleSoABuffer.constructView(1, cell2MieAVX.numParticles()), newton3);
    ljFunctorAVX.SoAFunctorPair(cell1AVX._particleSoABuffer.constructView(1, cell1AVX.numParticles()),
                                cell2AVX._particleSoABuffer.constructView(1, cell2AVX.numParticles()), newton3);
  } else {
    mieFunctorAVX.SoAFunctorPair(cell1MieAVX._particleSoABuffer, cell2MieAVX._particleSoABuffer, newton3);
    ljFunctorAVX.SoAFunctorPair(cell1AVX._particleSoABuffer, cell2AVX._particleSoABuffer, newton3);
  }
  ASSERT_TRUE(SoAParticlesEqual(cell1AVX._particleSoABuffer, cell1MieAVX._particleSoABuffer))
      << "Cells 1 not equal after applying functor.";
  ASSERT_TRUE(SoAParticlesEqual(cell2AVX._particleSoABuffer, cell2MieAVX._particleSoABuffer))
      << "Cells 2 not equal after applying functor.";

  ljFunctorAVX.SoAExtractor(cell1AVX, cell1AVX._particleSoABuffer, 0);
  ljFunctorAVX.SoAExtractor(cell2AVX, cell2AVX._particleSoABuffer, 0);
  ljFunctorAVX.SoAExtractor(cell1MieAVX, cell1MieAVX._particleSoABuffer, 0);
  ljFunctorAVX.SoAExtractor(cell2MieAVX, cell2MieAVX._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cell1AVX, cell1MieAVX)) << "Cells 1 not equal after extracting.";
  ASSERT_TRUE(AoSParticlesEqual(cell2AVX, cell2MieAVX)) << "Cells 2 not equal after extracting.";

  ljFunctorAVX.endTraversal(newton3);
  mieFunctorAVX.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorAVX.getPotentialEnergy(), mieFunctorAVX.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorAVX.getVirial(), mieFunctorAVX.getVirial(), tolerance) << "global virial";
}

void MieFunctorAVXTest::testMieFunctorAVXVSLJFunctorAVXOneCell(bool newton3, bool doDeleteSomeParticles,
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
  FMCell cellMieAVX(cellAVX);
  constexpr bool shifting = true;
  constexpr bool mixing = false;
  mdLib::MieFunctorAVX<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> mieFunctorAVX(_cutoff, 12, 6);
  mieFunctorAVX.setParticleProperties(_epsilon, _sigma * _sigma);
  mdLib::LJFunctorAVX<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorAVX(_cutoff);
  ljFunctorAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

  ASSERT_TRUE(AoSParticlesEqual(cellAVX, cellMieAVX)) << "Cells not equal after copy initialization.";

  ljFunctorAVX.initTraversal();
  mieFunctorAVX.initTraversal();

  mieFunctorAVX.SoALoader(cellMieAVX, cellMieAVX._particleSoABuffer, 0);
  ljFunctorAVX.SoALoader(cellAVX, cellAVX._particleSoABuffer, 0);

  ASSERT_TRUE(SoAParticlesEqual(cellAVX._particleSoABuffer, cellMieAVX._particleSoABuffer))
      << "Cells not equal after loading.";

  if (useUnalignedViews) {
    mieFunctorAVX.SoAFunctorSingle(cellMieAVX._particleSoABuffer.constructView(1, cellMieAVX.numParticles()), newton3);
    ljFunctorAVX.SoAFunctorSingle(cellAVX._particleSoABuffer.constructView(1, cellAVX.numParticles()), newton3);
  } else {
    mieFunctorAVX.SoAFunctorSingle(cellMieAVX._particleSoABuffer, newton3);
    ljFunctorAVX.SoAFunctorSingle(cellAVX._particleSoABuffer, newton3);
  }
  ASSERT_TRUE(SoAParticlesEqual(cellAVX._particleSoABuffer, cellMieAVX._particleSoABuffer))
      << "Cells not equal after applying functor.";

  ljFunctorAVX.SoAExtractor(cellAVX, cellAVX._particleSoABuffer, 0);
  ljFunctorAVX.SoAExtractor(cellMieAVX, cellMieAVX._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cellAVX, cellMieAVX)) << "Cells 1 not equal after extracting.";

  ljFunctorAVX.endTraversal(newton3);
  mieFunctorAVX.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorAVX.getPotentialEnergy(), mieFunctorAVX.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorAVX.getVirial(), mieFunctorAVX.getVirial(), tolerance) << "global virial";
}

void MieFunctorAVXTest::testMieFunctorAVXVSLJFunctorAVXVerlet(bool newton3, bool doDeleteSomeParticles) {
  using namespace autopas::utils::ArrayMath::literals;

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
      auto dr = cellAVX[i].getR() - cellAVX[j].getR();
      double dr2 = autopas::utils::ArrayMath::dot(dr, dr);
      if (dr2 <= _interactionLengthSquare) {
        neighborLists[i].push_back(j);
      }
    }
  }

  // copy cells
  FMCell cellMieAVX(cellAVX);
  constexpr bool shifting = true;
  constexpr bool mixing = false;
  constexpr bool calculateGlobals = true;
  mdLib::MieFunctorAVX<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, calculateGlobals> mieFunctorAVX(
      _cutoff, 12, 6);
  mieFunctorAVX.setParticleProperties(_epsilon, _sigma * _sigma);
  mdLib::LJFunctorAVX<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, calculateGlobals> ljFunctorAVX(
      _cutoff);
  ljFunctorAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

  ASSERT_TRUE(AoSParticlesEqual(cellAVX, cellMieAVX)) << "Cells not equal after copy initialization.";

  ljFunctorAVX.initTraversal();
  mieFunctorAVX.initTraversal();

  mieFunctorAVX.SoALoader(cellMieAVX, cellMieAVX._particleSoABuffer, 0);
  ljFunctorAVX.SoALoader(cellAVX, cellAVX._particleSoABuffer, 0);

  ASSERT_TRUE(SoAParticlesEqual(cellAVX._particleSoABuffer, cellMieAVX._particleSoABuffer))
      << "Cells not equal after loading.";

  for (size_t i = 0; i < numParticles; ++i) {
    mieFunctorAVX.SoAFunctorVerlet(cellMieAVX._particleSoABuffer, i, neighborLists[i], newton3);
    ljFunctorAVX.SoAFunctorVerlet(cellAVX._particleSoABuffer, i, neighborLists[i], newton3);
  }
  // Note that the calculated forces wont be compared, but the old forces. The knew forces are in t
  ASSERT_TRUE(SoAParticlesEqual(cellAVX._particleSoABuffer, cellMieAVX._particleSoABuffer))
      << "Cells not equal after applying functor.";

  ljFunctorAVX.SoAExtractor(cellAVX, cellAVX._particleSoABuffer, 0);
  mieFunctorAVX.SoAExtractor(cellMieAVX, cellMieAVX._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cellAVX, cellMieAVX)) << "Cells not equal after extracting.";

  ljFunctorAVX.endTraversal(newton3);
  mieFunctorAVX.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorAVX.getPotentialEnergy(), mieFunctorAVX.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorAVX.getVirial(), mieFunctorAVX.getVirial(), tolerance) << "global virial";
}

void MieFunctorAVXTest::testMieFunctorAVXAoS(bool newton3, bool doDeleteSomeParticles) {
  FMCell cellAVX;

  constexpr size_t numParticles = 7;

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::RandomGenerator::fillWithParticles(cellAVX, defaultParticle, _lowCorner, _highCorner,
                                                               numParticles);

  if (doDeleteSomeParticles) {
    // mark some particles as deleted to teqt if the functor handles them correctly
    for (auto &particle : cellAVX) {
      if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
    }
  }

  // copy cells
  FMCell cellMieAVX(cellAVX);
  constexpr bool shifting = true;
  constexpr bool mixing = false;
  mdLib::MieFunctorAVX<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> mieFunctorAVX(_cutoff, 12, 6);
  mieFunctorAVX.setParticleProperties(_epsilon, _sigma * _sigma);
  mdLib::LJFunctorAVX<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorAVX(_cutoff);
  ljFunctorAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

  ASSERT_TRUE(AoSParticlesEqual(cellAVX, cellMieAVX)) << "Cells not equal after copy initialization.";

  ljFunctorAVX.initTraversal();
  mieFunctorAVX.initTraversal();

  for (size_t i = 0; i < numParticles; ++i) {
    for (size_t j = newton3 ? i + 1 : 0; j < numParticles; ++j) {
      if (i == j) {
        continue;
      }
      mieFunctorAVX.AoSFunctor(cellMieAVX[i], cellMieAVX[j], newton3);
      ljFunctorAVX.AoSFunctor(cellAVX[i], cellAVX[j], newton3);
    }
  }

  ASSERT_TRUE(AoSParticlesEqual(cellAVX, cellMieAVX)) << "Cells not equal after applying AoSfunctor.";

  ljFunctorAVX.endTraversal(newton3);
  mieFunctorAVX.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorAVX.getPotentialEnergy(), mieFunctorAVX.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorAVX.getVirial(), mieFunctorAVX.getVirial(), tolerance) << "global virial";
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

TEST_P(MieFunctorAVXTest, testMieFunctorAVXVSLJFunctorAVXTwoCellsUseUnalignedViews) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testMieFunctorAVXAoS(newton3, doDeleteSomeParticle);
}

TEST_P(MieFunctorAVXTest, testMieFunctorAVXVSLJFunctorAVXVerlet) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testMieFunctorAVXVSLJFunctorAVXVerlet(newton3, doDeleteSomeParticle);
}

TEST_P(MieFunctorAVXTest, testMieFunctorAVXVSLJFunctorAVXOneCellAlignedAccess) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testMieFunctorAVXVSLJFunctorAVXOneCell(newton3, doDeleteSomeParticle, false);
}

TEST_P(MieFunctorAVXTest, testMieFunctorAVXVSLJFunctorAVXOneCellUseUnalignedViews) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testMieFunctorAVXVSLJFunctorAVXOneCell(newton3, doDeleteSomeParticle, true);
}

TEST_P(MieFunctorAVXTest, testMieFunctorAVXVSLJFunctorAVXTwoCellsAlignedAccess) {
  auto [newton3, doDeleteSomeParticle] = GetParam();
  testMieFunctorAVXVSLJFunctorAVXTwoCells(newton3, doDeleteSomeParticle, false);
}

INSTANTIATE_TEST_SUITE_P(Generated, MieFunctorAVXTest, ::testing::Combine(::testing::Bool(), ::testing::Bool()),
                         toString);
#endif
