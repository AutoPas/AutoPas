/**
 * @file LJFunctorTestHWY.cpp
 * @author Luis Gall
 * @date 04/23/24
 */

#include "LJFunctorTestHWY.h"

#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/Particle.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "molecularDynamicsLibrary/LJFunctorAVX.h"

template <class SoAType>
bool LJFunctorTestHWY::SoAParticlesEqual(autopas::SoA<SoAType> &soa1, autopas::SoA<SoAType> &soa2) {
    
    EXPECT_GT(soa1.size(), 0);
    EXPECT_EQ(soa1.size(), soa2.size());

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

    for (size_t i = 0; i < soa1.size(); ++i) {
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

bool LJFunctorTestHWY::particleEqual(Particle &p1, Particle &p2) {
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

bool LJFunctorTestHWY::AoSParticlesEqual(FMCell &cell1, FMCell &cell2) {
    EXPECT_GT(cell1.size(), 0);
    EXPECT_EQ(cell1.size(), cell2.size());

    bool ret = true;
    for (size_t i = 0; i < cell1.size(); ++i) {
        ret = ret and particleEqual(cell1._particles[i], cell2._particles[i]);
    }

    return ret;
}

template <VectorizationPattern vecPattern>
void LJFunctorTestHWY::testLJFunctorAVXvsLJFunctorHWYTwoCells(bool newton3, bool doDeleteSomeParticles, bool useUnalignedViews) {
    FMCell cell1HWY;
    FMCell cell2HWY;

    size_t numParticles = 8;

    Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
    autopasTools::generators::RandomGenerator::fillWithParticles(
        cell1HWY, defaultParticle, _lowCorner, {_highCorner[0] / 2, _highCorner[1], _highCorner[2]}, numParticles);
    autopasTools::generators::RandomGenerator::fillWithParticles(
        cell2HWY, defaultParticle, {_highCorner[0] / 2, _lowCorner[1], _lowCorner[2]}, _highCorner, numParticles);

    if (doDeleteSomeParticles) {
        for (auto &particle : cell1HWY) {
        if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
        }
        for (auto &particle : cell2HWY) {
        if (particle.getID() == 4) autopas::internal::markParticleAsDeleted(particle);
        }
    }

    // copy cells
    FMCell cell1NoHWY(cell1HWY);
    FMCell cell2NoHWY(cell2HWY);

    constexpr bool shifting = true;
    constexpr bool mixing = false;
    mdLib::LJFunctorAVX<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorAVX(_cutoff);
    ljFunctorAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    mdLib::LJFunctorHWY<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true, false, vecPattern> ljFunctorHWY(_cutoff);
    ljFunctorHWY.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

    ljFunctorHWY.initTraversal();
    ljFunctorAVX.initTraversal();

    ASSERT_TRUE(AoSParticlesEqual(cell1HWY, cell1NoHWY)) << "Cells 1 not equal after copy initialization.";
    ASSERT_TRUE(AoSParticlesEqual(cell2HWY, cell2NoHWY)) << "Cells 2 not equal after copy initialization.";

    ljFunctorAVX.SoALoader(cell1NoHWY, cell1NoHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
    ljFunctorAVX.SoALoader(cell2NoHWY, cell2NoHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
    ljFunctorHWY.SoALoader(cell1HWY, cell1HWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
    ljFunctorHWY.SoALoader(cell2HWY, cell2HWY._particleSoABuffer, 0, /*skipSoAResize*/ false);

    ASSERT_TRUE(SoAParticlesEqual(cell1HWY._particleSoABuffer, cell1NoHWY._particleSoABuffer))
        << "Cells 1 not equal after loading.";
    ASSERT_TRUE(SoAParticlesEqual(cell2HWY._particleSoABuffer, cell2NoHWY._particleSoABuffer))
        << "Cells 2 not equal after loading.";

    if (useUnalignedViews) {
        ljFunctorAVX.SoAFunctorPair(cell1NoHWY._particleSoABuffer.constructView(1, cell1NoHWY.size()),
                                    cell2NoHWY._particleSoABuffer.constructView(1, cell2NoHWY.size()), newton3);
        ljFunctorHWY.SoAFunctorPair(cell1HWY._particleSoABuffer.constructView(1, cell1HWY.size()),
                                    cell2HWY._particleSoABuffer.constructView(1, cell2HWY.size()), newton3);
    } else {
        ljFunctorAVX.SoAFunctorPair(cell1NoHWY._particleSoABuffer, cell2NoHWY._particleSoABuffer, newton3);
        ljFunctorHWY.SoAFunctorPair(cell1HWY._particleSoABuffer, cell2HWY._particleSoABuffer, newton3);
    }
    ASSERT_TRUE(SoAParticlesEqual(cell1HWY._particleSoABuffer, cell1NoHWY._particleSoABuffer))
        << "Cells 1 not equal after applying functor.";
    ASSERT_TRUE(SoAParticlesEqual(cell2HWY._particleSoABuffer, cell2NoHWY._particleSoABuffer))
        << "Cells 2 not equal after applying functor.";

    ljFunctorHWY.SoAExtractor(cell1HWY, cell1HWY._particleSoABuffer, 0);
    ljFunctorHWY.SoAExtractor(cell2HWY, cell2HWY._particleSoABuffer, 0);
    ljFunctorHWY.SoAExtractor(cell1NoHWY, cell1NoHWY._particleSoABuffer, 0);
    ljFunctorHWY.SoAExtractor(cell2NoHWY, cell2NoHWY._particleSoABuffer, 0);

    ASSERT_TRUE(AoSParticlesEqual(cell1HWY, cell1NoHWY)) << "Cells 1 not equal after extracting.";
    ASSERT_TRUE(AoSParticlesEqual(cell2HWY, cell2NoHWY)) << "Cells 2 not equal after extracting.";

    ljFunctorHWY.endTraversal(newton3);
    ljFunctorAVX.endTraversal(newton3);

    double tolerance = 1e-8;
    EXPECT_NEAR(ljFunctorAVX.getPotentialEnergy(), ljFunctorHWY.getUpot(), tolerance) << "global uPot";
    EXPECT_NEAR(ljFunctorAVX.getVirial(), ljFunctorHWY.getVirial(), tolerance) << "global virial";
}

template <VectorizationPattern vecPattern>
void LJFunctorTestHWY::testLJFunctorAVXvsLJFunctorHWYOneCell(bool newton3, bool doDeleteSomeParticles, bool useUnalignedViews) {
    FMCell cellHWY;

    size_t numParticles = 8;

    Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
    autopasTools::generators::RandomGenerator::fillWithParticles(cellHWY, defaultParticle, _lowCorner, _highCorner,
                                                                numParticles);

    if (doDeleteSomeParticles) {
        for (auto &particle : cellHWY) {
            if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
        }
    }

    // copy cells
    FMCell cellNoHWY(cellHWY);
    constexpr bool shifting = true;
    constexpr bool mixing = false;
    mdLib::LJFunctorAVX<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorAVX(_cutoff);
    ljFunctorAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    mdLib::LJFunctorHWY<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true, false, vecPattern> ljFunctorHWY(_cutoff);
    ljFunctorHWY.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

    ASSERT_TRUE(AoSParticlesEqual(cellHWY, cellNoHWY)) << "Cells not equal after copy initialization.";

    ljFunctorHWY.initTraversal();
    ljFunctorAVX.initTraversal();

    ljFunctorAVX.SoALoader(cellNoHWY, cellNoHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
    ljFunctorHWY.SoALoader(cellHWY, cellHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);

    ASSERT_TRUE(SoAParticlesEqual(cellHWY._particleSoABuffer, cellNoHWY._particleSoABuffer))
        << "Cells not equal after loading.";

    if (useUnalignedViews) {
        ljFunctorAVX.SoAFunctorSingle(cellNoHWY._particleSoABuffer.constructView(1, cellNoHWY.size()), newton3);
        ljFunctorHWY.SoAFunctorSingle(cellHWY._particleSoABuffer.constructView(1, cellHWY.size()), newton3);
    } else {
        ljFunctorAVX.SoAFunctorSingle(cellNoHWY._particleSoABuffer, newton3);
        ljFunctorHWY.SoAFunctorSingle(cellHWY._particleSoABuffer, newton3);
    }
    ASSERT_TRUE(SoAParticlesEqual(cellHWY._particleSoABuffer, cellNoHWY._particleSoABuffer))
        << "Cells not equal after applying functor.";

    ljFunctorHWY.SoAExtractor(cellHWY, cellHWY._particleSoABuffer, 0);
    ljFunctorHWY.SoAExtractor(cellNoHWY, cellNoHWY._particleSoABuffer, 0);

    ASSERT_TRUE(AoSParticlesEqual(cellHWY, cellNoHWY)) << "Cells 1 not equal after extracting.";

    ljFunctorHWY.endTraversal(newton3);
    ljFunctorAVX.endTraversal(newton3);

    double tolerance = 1e-8;
    EXPECT_NEAR(ljFunctorAVX.getPotentialEnergy(), ljFunctorHWY.getUpot(), tolerance) << "global uPot";
    EXPECT_NEAR(ljFunctorAVX.getVirial(), ljFunctorHWY.getVirial(), tolerance) << "global virial";
}

template <VectorizationPattern vecPattern>
void LJFunctorTestHWY::testLJFunctorAVXvsLJFunctorHWYVerlet(bool newton3, bool doDeleteSomeParticles) {
    
    using namespace autopas::utils::ArrayMath::literals;

    FMCell cellAVX;

    constexpr size_t numParticles = 8;

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
    FMCell cellHWY(cellAVX);
    constexpr bool shifting = true;
    constexpr bool mixing = false;
    constexpr bool calculateGlobals = true;
    mdLib::LJFunctorHWY<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true, false, vecPattern> ljFunctorHWY(_cutoff);
    ljFunctorHWY.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    mdLib::LJFunctorAVX<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, calculateGlobals> ljFunctorAVX(
        _cutoff);
    ljFunctorAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

    ASSERT_TRUE(AoSParticlesEqual(cellAVX, cellHWY)) << "Cells not equal after copy initialization.";

    ljFunctorAVX.initTraversal();
    ljFunctorHWY.initTraversal();

    ljFunctorHWY.SoALoader(cellHWY, cellHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
    ljFunctorAVX.SoALoader(cellAVX, cellAVX._particleSoABuffer, 0, /*skipSoAResize*/ false);

    ASSERT_TRUE(SoAParticlesEqual(cellAVX._particleSoABuffer, cellHWY._particleSoABuffer))
        << "Cells not equal after loading.";

    for (size_t i = 0; i < numParticles; ++i) {
        ljFunctorHWY.SoAFunctorVerlet(cellHWY._particleSoABuffer, i, neighborLists[i], newton3);
        ljFunctorAVX.SoAFunctorVerlet(cellAVX._particleSoABuffer, i, neighborLists[i], newton3);
    }

    ASSERT_TRUE(SoAParticlesEqual(cellAVX._particleSoABuffer, cellHWY._particleSoABuffer))
        << "Cells not equal after applying functor.";

    ljFunctorAVX.SoAExtractor(cellAVX, cellAVX._particleSoABuffer, 0);
    ljFunctorAVX.SoAExtractor(cellHWY, cellHWY._particleSoABuffer, 0);

    ASSERT_TRUE(AoSParticlesEqual(cellAVX, cellHWY)) << "Cells not equal after extracting.";

    ljFunctorAVX.endTraversal(newton3);
    ljFunctorHWY.endTraversal(newton3);

    double tolerance = 1e-8;
    EXPECT_NEAR(ljFunctorAVX.getPotentialEnergy(), ljFunctorHWY.getUpot(), tolerance) << "global uPot";
    EXPECT_NEAR(ljFunctorAVX.getVirial(), ljFunctorHWY.getVirial(), tolerance) << "global virial";
}

void LJFunctorTestHWY::testLJFunctorAVXvsLJFunctorHWYAoS(bool newton3, bool doDeleteSomeParticles) {
    FMCell cellHWY;

    constexpr size_t numParticles = 8;

    Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
    autopasTools::generators::RandomGenerator::fillWithParticles(cellHWY, defaultParticle, _lowCorner, _highCorner,
                                                                numParticles);

    if (doDeleteSomeParticles) {
        // mark some particles as deleted to test if the functor handles them correctly
        for (auto &particle : cellHWY) {
        if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
        }
    }

    // copy cells
    FMCell cellNoHWY(cellHWY);
    constexpr bool shifting = true;
    constexpr bool mixing = false;
    mdLib::LJFunctorAVX<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorAVX(_cutoff);
    ljFunctorAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    mdLib::LJFunctorHWY<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true, false> ljFunctorHWY(_cutoff);
    ljFunctorHWY.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

    ASSERT_TRUE(AoSParticlesEqual(cellHWY, cellNoHWY)) << "Cells not equal after copy initialization.";

    ljFunctorHWY.initTraversal();
    ljFunctorAVX.initTraversal();

    for (size_t i = 0; i < numParticles; ++i) {
        for (size_t j = newton3 ? i + 1 : 0; j < numParticles; ++j) {
            if (i == j) {
                continue;
            }
            ljFunctorAVX.AoSFunctor(cellNoHWY[i], cellNoHWY[j], newton3);
            ljFunctorHWY.AoSFunctor(cellHWY[i], cellHWY[j], newton3);
        }
    }

    ASSERT_TRUE(AoSParticlesEqual(cellHWY, cellNoHWY)) << "Cells not equal after applying AoSfunctor.";

    ljFunctorHWY.endTraversal(newton3);
    ljFunctorAVX.endTraversal(newton3);

    double tolerance = 1e-8;
    EXPECT_NEAR(ljFunctorHWY.getUpot(), ljFunctorAVX.getPotentialEnergy(), tolerance) << "global uPot";
    EXPECT_NEAR(ljFunctorHWY.getVirial(), ljFunctorAVX.getVirial(), tolerance) << "global virial";
}

TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYAoS) {
  auto [newton3, doDeleteSomeParticle, _] = GetParam();
  testLJFunctorAVXvsLJFunctorHWYAoS(newton3, doDeleteSomeParticle);
}

TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYVerlet) {
  // different vectorization patterns are currently not supported for Verlet Functor
  auto [newton3, doDeleteSomeParticle, _] = GetParam();
  testLJFunctorAVXvsLJFunctorHWYVerlet<VectorizationPattern::p1xVec>(newton3, doDeleteSomeParticle);
}

TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYOneCellAlignedAccess) {
  auto [newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  switch (vecPattern)
  {
  case VectorizationPattern::p1xVec:
    testLJFunctorAVXvsLJFunctorHWYOneCell<VectorizationPattern::p1xVec>(newton3, doDeleteSomeParticle, false);
    break;
  case VectorizationPattern::p2xVecDiv2:
    testLJFunctorAVXvsLJFunctorHWYOneCell<VectorizationPattern::p2xVecDiv2>(newton3, doDeleteSomeParticle, false);
    break;
  case VectorizationPattern::pVecDiv2x2:
    testLJFunctorAVXvsLJFunctorHWYOneCell<VectorizationPattern::pVecDiv2x2>(newton3, doDeleteSomeParticle, false);
    break;
  case VectorizationPattern::pVecx1:
    testLJFunctorAVXvsLJFunctorHWYOneCell<VectorizationPattern::pVecx1>(newton3, doDeleteSomeParticle, false);
    break;
  case VectorizationPattern::pVecxVec:
    testLJFunctorAVXvsLJFunctorHWYOneCell<VectorizationPattern::pVecxVec>(newton3, doDeleteSomeParticle, false);
    break;
  default:
    throw std::runtime_error("No vectorization pattern matched");
  }
}

TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYOneCellUseUnalignedViews) {
  auto [newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  switch (vecPattern)
  {
  case VectorizationPattern::p1xVec:
    testLJFunctorAVXvsLJFunctorHWYOneCell<VectorizationPattern::p1xVec>(newton3, doDeleteSomeParticle, true);
    break;
  case VectorizationPattern::p2xVecDiv2:
    testLJFunctorAVXvsLJFunctorHWYOneCell<VectorizationPattern::p2xVecDiv2>(newton3, doDeleteSomeParticle, true);
    break;
  case VectorizationPattern::pVecDiv2x2:
    testLJFunctorAVXvsLJFunctorHWYOneCell<VectorizationPattern::pVecDiv2x2>(newton3, doDeleteSomeParticle, true);
    break;
  case VectorizationPattern::pVecx1:
    testLJFunctorAVXvsLJFunctorHWYOneCell<VectorizationPattern::pVecx1>(newton3, doDeleteSomeParticle, true);
    break;
  case VectorizationPattern::pVecxVec:
    testLJFunctorAVXvsLJFunctorHWYOneCell<VectorizationPattern::pVecxVec>(newton3, doDeleteSomeParticle, true);
    break;
  default:
    throw std::runtime_error("No vectorization pattern matched");
  }
  
}

TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsAlignedAccess) {
  auto [newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  switch (vecPattern)
  {
  case VectorizationPattern::p1xVec:
    testLJFunctorAVXvsLJFunctorHWYTwoCells<VectorizationPattern::p1xVec>(newton3, doDeleteSomeParticle, false);
    break;
  case VectorizationPattern::p2xVecDiv2:
    testLJFunctorAVXvsLJFunctorHWYTwoCells<VectorizationPattern::p2xVecDiv2>(newton3, doDeleteSomeParticle, false);
    break;
  case VectorizationPattern::pVecDiv2x2:
    testLJFunctorAVXvsLJFunctorHWYTwoCells<VectorizationPattern::pVecDiv2x2>(newton3, doDeleteSomeParticle, false);
    break;
  case VectorizationPattern::pVecx1:
    testLJFunctorAVXvsLJFunctorHWYTwoCells<VectorizationPattern::pVecx1>(newton3, doDeleteSomeParticle, false);
    break;
  case VectorizationPattern::pVecxVec:
    testLJFunctorAVXvsLJFunctorHWYTwoCells<VectorizationPattern::pVecxVec>(newton3, doDeleteSomeParticle, false);
    break;
  default:
    throw std::runtime_error("No vectorization pattern matched");
  }
}

TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsUseUnalignedViews) {
  auto [newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  switch (vecPattern)
  {
  case VectorizationPattern::p1xVec:
    testLJFunctorAVXvsLJFunctorHWYTwoCells<VectorizationPattern::p1xVec>(newton3, doDeleteSomeParticle, true);
    break;
  case VectorizationPattern::p2xVecDiv2:
    testLJFunctorAVXvsLJFunctorHWYTwoCells<VectorizationPattern::p2xVecDiv2>(newton3, doDeleteSomeParticle, true);
    break;
  case VectorizationPattern::pVecDiv2x2:
    testLJFunctorAVXvsLJFunctorHWYTwoCells<VectorizationPattern::pVecDiv2x2>(newton3, doDeleteSomeParticle, true);
    break;
  case VectorizationPattern::pVecx1:
    testLJFunctorAVXvsLJFunctorHWYTwoCells<VectorizationPattern::pVecx1>(newton3, doDeleteSomeParticle, true);
    break;
  case VectorizationPattern::pVecxVec:
    testLJFunctorAVXvsLJFunctorHWYTwoCells<VectorizationPattern::pVecxVec>(newton3, doDeleteSomeParticle, true);
    break;
  default:
    throw std::runtime_error("No vectorization pattern matched");
  }
}

std::vector<VectorizationPattern> patterns {
    VectorizationPattern::p1xVec,
    VectorizationPattern::p2xVecDiv2,
    // VectorizationPattern::pVecDiv2x2,
    // VectorizationPattern::pVecx1,
};

std::map<VectorizationPattern, std::string> patternsToString {
    { VectorizationPattern::p1xVec, "1xVec"},
    { VectorizationPattern::p2xVecDiv2, "2xVec_2"},
    { VectorizationPattern::pVecDiv2x2, "Vec_2x2" },
    { VectorizationPattern::pVecx1, "Vecx1" },
};

static auto toString = [](const auto &info) {
  auto [newton3, doDeleteSomeParticle, vecPattern] = info.param;
  std::stringstream resStream;
  resStream << patternsToString[vecPattern] << (newton3 ? "N3" : "noN3") << "_" << (doDeleteSomeParticle ? "withDeletions" : "noDeletions");
  std::string res = resStream.str();
  std::replace(res.begin(), res.end(), '-', '_');
  std::replace(res.begin(), res.end(), '.', '_');
  return res;
};

INSTANTIATE_TEST_SUITE_P(Generated, LJFunctorTestHWY, ::testing::Combine(::testing::Bool(), ::testing::Bool(), ::testing::ValuesIn(patterns)), toString);