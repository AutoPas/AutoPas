/**
 * @file LJFunctorAVXTest.cpp
 * @author F. Gratl
 * @date 12/17/18
 */

#ifdef __AVX__

#include "LJFunctorAVXTest.h"

#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/ParticleDefinitions.h"
#include "autopasTools/generators/UniformGenerator.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "molecularDynamicsLibrary/LJFunctorAVX.h"

template <class SoAType>
bool LJFunctorAVXTest::SoAParticlesEqual(autopas::SoA<SoAType> &soa1, autopas::SoA<SoAType> &soa2) {
  EXPECT_GT(soa1.size(), 0);
  EXPECT_EQ(soa1.size(), soa2.size());

  unsigned long *const __restrict idptr1 = soa1.template begin<Molecule::AttributeNames::id>();
  unsigned long *const __restrict idptr2 = soa2.template begin<Molecule::AttributeNames::id>();

  double *const __restrict xptr1 = soa1.template begin<Molecule::AttributeNames::posX>();
  double *const __restrict yptr1 = soa1.template begin<Molecule::AttributeNames::posY>();
  double *const __restrict zptr1 = soa1.template begin<Molecule::AttributeNames::posZ>();
  double *const __restrict xptr2 = soa2.template begin<Molecule::AttributeNames::posX>();
  double *const __restrict yptr2 = soa2.template begin<Molecule::AttributeNames::posY>();
  double *const __restrict zptr2 = soa2.template begin<Molecule::AttributeNames::posZ>();

  double *const __restrict fxptr1 = soa1.template begin<Molecule::AttributeNames::forceX>();
  double *const __restrict fyptr1 = soa1.template begin<Molecule::AttributeNames::forceY>();
  double *const __restrict fzptr1 = soa1.template begin<Molecule::AttributeNames::forceZ>();
  double *const __restrict fxptr2 = soa2.template begin<Molecule::AttributeNames::forceX>();
  double *const __restrict fyptr2 = soa2.template begin<Molecule::AttributeNames::forceY>();
  double *const __restrict fzptr2 = soa2.template begin<Molecule::AttributeNames::forceZ>();

  for (size_t i = 0; i < soa1.size(); ++i) {
    EXPECT_EQ(idptr1[i], idptr2[i]);

    double tolerance = 2e-8;
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

bool LJFunctorAVXTest::particleEqual(Molecule &p1, Molecule &p2) {
  EXPECT_EQ(p1.getID(), p2.getID());

  double tolerance = 2e-8;

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
  EXPECT_GT(cell1.size(), 0);
  EXPECT_EQ(cell1.size(), cell2.size());

  bool ret = true;
  for (size_t i = 0; i < cell1.size(); ++i) {
    ret = ret and particleEqual(cell1._particles[i], cell2._particles[i]);
  }

  return ret;
}

template <bool mixing>
void LJFunctorAVXTest::testLJFunctorVSLJFunctorAVXTwoCells(bool newton3, bool doDeleteSomeParticles,
                                                           bool useUnalignedViews) {
  FMCell cell1AVX;
  FMCell cell2AVX;

  size_t numParticles = 7;

  ParticlePropertiesLibrary<double, size_t> PPL{_cutoff};
  if constexpr (mixing) {
    PPL.addSiteType(0, 1.);
    PPL.addLJParametersToSite(0, 1., 1.);
    PPL.addSiteType(1, 1.5);
    PPL.addLJParametersToSite(1, 2., 1.);
    PPL.addSiteType(2, 2.);
    PPL.addLJParametersToSite(2, 1., 1.);
    PPL.addSiteType(3, 2.5);
    PPL.addLJParametersToSite(3, 2., 1.);
    PPL.addSiteType(4, 3.);
    PPL.addLJParametersToSite(4, 1., 1.);
    PPL.calculateMixingCoefficients();
  }

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(
      cell1AVX, defaultParticle, _lowCorner, {_highCorner[0] / 2, _highCorner[1], _highCorner[2]}, numParticles);
  autopasTools::generators::UniformGenerator::fillWithParticles(
      cell2AVX, defaultParticle, {_highCorner[0] / 2, _lowCorner[1], _lowCorner[2]}, _highCorner, numParticles);

  for (auto &particle : cell1AVX) {
    if (doDeleteSomeParticles) {
      if (particle.getID() == 3) {
        autopas::internal::markParticleAsDeleted(particle);
      }
    }
    if constexpr (mixing) {
      particle.setTypeId(particle.getID() % 5);
    }
  }
  for (auto &particle : cell2AVX) {
    if (doDeleteSomeParticles) {
      if (particle.getID() == 4) {
        autopas::internal::markParticleAsDeleted(particle);
      }
    }
    if constexpr (mixing) {
      particle.setTypeId(particle.getID() % 5);
    }
  }

  // copy cells
  FMCell cell1NoAVX(cell1AVX);
  FMCell cell2NoAVX(cell2AVX);

  constexpr bool shifting = true;

  auto ljFunctorNoAVX = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctor<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctor<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto ljFunctorAVX = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorAVX<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctorAVX<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  if constexpr (not mixing) {
    ljFunctorNoAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    ljFunctorAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  }

  ljFunctorAVX.initTraversal();
  ljFunctorNoAVX.initTraversal();

  ASSERT_TRUE(AoSParticlesEqual(cell1AVX, cell1NoAVX)) << "Cells 1 not equal after copy initialization.";
  ASSERT_TRUE(AoSParticlesEqual(cell2AVX, cell2NoAVX)) << "Cells 2 not equal after copy initialization.";

  ljFunctorNoAVX.SoALoader(cell1NoAVX, cell1NoAVX._particleSoABuffer, 0, /*skipSoAResize*/ false);
  ljFunctorNoAVX.SoALoader(cell2NoAVX, cell2NoAVX._particleSoABuffer, 0, /*skipSoAResize*/ false);
  ljFunctorAVX.SoALoader(cell1AVX, cell1AVX._particleSoABuffer, 0, /*skipSoAResize*/ false);
  ljFunctorAVX.SoALoader(cell2AVX, cell2AVX._particleSoABuffer, 0, /*skipSoAResize*/ false);

  ASSERT_TRUE(SoAParticlesEqual(cell1AVX._particleSoABuffer, cell1NoAVX._particleSoABuffer))
      << "Cells 1 not equal after loading.";
  ASSERT_TRUE(SoAParticlesEqual(cell2AVX._particleSoABuffer, cell2NoAVX._particleSoABuffer))
      << "Cells 2 not equal after loading.";

  if (useUnalignedViews) {
    ljFunctorNoAVX.SoAFunctorPair(cell1NoAVX._particleSoABuffer.constructView(1, cell1NoAVX.size()),
                                  cell2NoAVX._particleSoABuffer.constructView(1, cell2NoAVX.size()), newton3);
    ljFunctorAVX.SoAFunctorPair(cell1AVX._particleSoABuffer.constructView(1, cell1AVX.size()),
                                cell2AVX._particleSoABuffer.constructView(1, cell2AVX.size()), newton3);
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
  EXPECT_NEAR(ljFunctorAVX.getPotentialEnergy(), ljFunctorNoAVX.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorAVX.getVirial(), ljFunctorNoAVX.getVirial(), tolerance) << "global virial";
}

template <bool mixing>
void LJFunctorAVXTest::testLJFunctorVSLJFunctorAVXOneCell(bool newton3, bool doDeleteSomeParticles,
                                                          bool useUnalignedViews) {
  FMCell cellAVX;

  size_t numParticles = 7;

  ParticlePropertiesLibrary<double, size_t> PPL{_cutoff};
  if constexpr (mixing) {
    PPL.addSiteType(0, 1.);
    PPL.addLJParametersToSite(0, 1., 1.);
    PPL.addSiteType(1, 1.5);
    PPL.addLJParametersToSite(1, 2., 1.);
    PPL.addSiteType(2, 2.);
    PPL.addLJParametersToSite(2, 1., 1.);
    PPL.addSiteType(3, 2.5);
    PPL.addLJParametersToSite(3, 2., 1.);
    PPL.addSiteType(4, 3.);
    PPL.addLJParametersToSite(4, 1., 1.);
    PPL.calculateMixingCoefficients();
  }

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cellAVX, defaultParticle, _lowCorner, _highCorner,
                                                                numParticles);

  for (auto &particle : cellAVX) {
    if (doDeleteSomeParticles) {
      if (particle.getID() == 3) {
        autopas::internal::markParticleAsDeleted(particle);
      }
    }
    if constexpr (mixing) {
      particle.setTypeId(particle.getID() % 5);
    }
  }

  // copy cells
  FMCell cellNoAVX(cellAVX);
  constexpr bool shifting = true;

  auto ljFunctorNoAVX = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctor<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctor<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto ljFunctorAVX = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorAVX<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctorAVX<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  if constexpr (not mixing) {
    ljFunctorNoAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    ljFunctorAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  }

  ASSERT_TRUE(AoSParticlesEqual(cellAVX, cellNoAVX)) << "Cells not equal after copy initialization.";

  ljFunctorAVX.initTraversal();
  ljFunctorNoAVX.initTraversal();

  ljFunctorNoAVX.SoALoader(cellNoAVX, cellNoAVX._particleSoABuffer, 0, /*skipSoAResize*/ false);
  ljFunctorAVX.SoALoader(cellAVX, cellAVX._particleSoABuffer, 0, /*skipSoAResize*/ false);

  ASSERT_TRUE(SoAParticlesEqual(cellAVX._particleSoABuffer, cellNoAVX._particleSoABuffer))
      << "Cells not equal after loading.";

  if (useUnalignedViews) {
    ljFunctorNoAVX.SoAFunctorSingle(cellNoAVX._particleSoABuffer.constructView(1, cellNoAVX.size()), newton3);
    ljFunctorAVX.SoAFunctorSingle(cellAVX._particleSoABuffer.constructView(1, cellAVX.size()), newton3);
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
  EXPECT_NEAR(ljFunctorAVX.getPotentialEnergy(), ljFunctorNoAVX.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorAVX.getVirial(), ljFunctorNoAVX.getVirial(), tolerance) << "global virial";
}

template <bool mixing>
void LJFunctorAVXTest::testLJFunctorVSLJFunctorAVXVerlet(bool newton3, bool doDeleteSomeParticles) {
  using namespace autopas::utils::ArrayMath::literals;

  FMCell cellAVX;

  constexpr size_t numParticles = 7;

  ParticlePropertiesLibrary<double, size_t> PPL{_cutoff};
  if constexpr (mixing) {
    PPL.addSiteType(0, 1.);
    PPL.addLJParametersToSite(0, 1., 1.);
    PPL.addSiteType(1, 1.5);
    PPL.addLJParametersToSite(1, 2., 1.);
    PPL.addSiteType(2, 2.);
    PPL.addLJParametersToSite(2, 1., 1.);
    PPL.addSiteType(3, 2.5);
    PPL.addLJParametersToSite(3, 2., 1.);
    PPL.addSiteType(4, 3.);
    PPL.addLJParametersToSite(4, 1., 1.);
    PPL.calculateMixingCoefficients();
  }

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cellAVX, defaultParticle, _lowCorner, _highCorner,
                                                                numParticles);

  for (auto &particle : cellAVX) {
    if (doDeleteSomeParticles) {
      if (particle.getID() == 3) {
        autopas::internal::markParticleAsDeleted(particle);
      }
    }
    if constexpr (mixing) {
      particle.setTypeId(particle.getID() % 5);
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
  FMCell cellNoAVX(cellAVX);

  constexpr bool shifting = true;
  constexpr bool calculateGlobals = true;

  auto ljFunctorNoAVX = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctor<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctor<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto ljFunctorAVX = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorAVX<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctorAVX<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  if constexpr (not mixing) {
    ljFunctorNoAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    ljFunctorAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  }

  ASSERT_TRUE(AoSParticlesEqual(cellAVX, cellNoAVX)) << "Cells not equal after copy initialization.";

  ljFunctorAVX.initTraversal();
  ljFunctorNoAVX.initTraversal();

  ljFunctorNoAVX.SoALoader(cellNoAVX, cellNoAVX._particleSoABuffer, 0, /*skipSoAResize*/ false);
  ljFunctorAVX.SoALoader(cellAVX, cellAVX._particleSoABuffer, 0, /*skipSoAResize*/ false);

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
  EXPECT_NEAR(ljFunctorAVX.getPotentialEnergy(), ljFunctorNoAVX.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorAVX.getVirial(), ljFunctorNoAVX.getVirial(), tolerance) << "global virial";
}

template <bool mixing>
void LJFunctorAVXTest::testLJFunctorVSLJFunctorAVXAoS(bool newton3, bool doDeleteSomeParticles) {
  FMCell cellAVX;

  constexpr size_t numParticles = 7;

  ParticlePropertiesLibrary<double, size_t> PPL{_cutoff};
  if constexpr (mixing) {
    PPL.addSiteType(0, 1.);
    PPL.addLJParametersToSite(0, 1., 1.);
    PPL.addSiteType(1, 1.5);
    PPL.addLJParametersToSite(1, 2., 1.);
    PPL.addSiteType(2, 2.);
    PPL.addLJParametersToSite(2, 1., 1.);
    PPL.addSiteType(3, 2.5);
    PPL.addLJParametersToSite(3, 2., 1.);
    PPL.addSiteType(4, 3.);
    PPL.addLJParametersToSite(4, 1., 1.);
    PPL.calculateMixingCoefficients();
  }

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cellAVX, defaultParticle, _lowCorner, _highCorner,
                                                                numParticles);

  for (auto &particle : cellAVX) {
    if (doDeleteSomeParticles) {
      if (particle.getID() == 3) {
        autopas::internal::markParticleAsDeleted(particle);
      }
    }
    if constexpr (mixing) {
      particle.setTypeId(particle.getID() % 5);
    }
  }

  // copy cells
  FMCell cellNoAVX(cellAVX);
  constexpr bool shifting = true;

  auto ljFunctorNoAVX = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctor<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctor<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto ljFunctorAVX = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorAVX<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctorAVX<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  if constexpr (not mixing) {
    ljFunctorNoAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    ljFunctorAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  }

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
  EXPECT_NEAR(ljFunctorAVX.getPotentialEnergy(), ljFunctorNoAVX.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorAVX.getVirial(), ljFunctorNoAVX.getVirial(), tolerance) << "global virial";
}

TEST_P(LJFunctorAVXTest, testLJFunctorVSLJFunctorAVXAoS) {
  const auto [mixing, newton3, doDeleteSomeParticle] = GetParam();
  if (mixing) {
    testLJFunctorVSLJFunctorAVXAoS<true>(newton3, doDeleteSomeParticle);
  } else {
    testLJFunctorVSLJFunctorAVXAoS<false>(newton3, doDeleteSomeParticle);
  }
}

TEST_P(LJFunctorAVXTest, testLJFunctorVSLJFunctorAVXVerlet) {
  const auto [mixing, newton3, doDeleteSomeParticle] = GetParam();
  if (mixing) {
    testLJFunctorVSLJFunctorAVXVerlet<true>(newton3, doDeleteSomeParticle);
  } else {
    testLJFunctorVSLJFunctorAVXVerlet<false>(newton3, doDeleteSomeParticle);
  }
}

TEST_P(LJFunctorAVXTest, testLJFunctorVSLJFunctorAVXOneCellAlignedAccess) {
  const auto [mixing, newton3, doDeleteSomeParticle] = GetParam();
  if (mixing) {
    testLJFunctorVSLJFunctorAVXOneCell<true>(newton3, doDeleteSomeParticle, false);
  } else {
    testLJFunctorVSLJFunctorAVXOneCell<false>(newton3, doDeleteSomeParticle, false);
  }
}

TEST_P(LJFunctorAVXTest, testLJFunctorVSLJFunctorAVXOneCellUseUnalignedViews) {
  const auto [mixing, newton3, doDeleteSomeParticle] = GetParam();
  if (mixing) {
    testLJFunctorVSLJFunctorAVXOneCell<true>(newton3, doDeleteSomeParticle, true);
  } else {
    testLJFunctorVSLJFunctorAVXOneCell<false>(newton3, doDeleteSomeParticle, true);
  }
}

TEST_P(LJFunctorAVXTest, testLJFunctorVSLJFunctorAVXTwoCellsAlignedAccess) {
  const auto [mixing, newton3, doDeleteSomeParticle] = GetParam();
  if (mixing) {
    testLJFunctorVSLJFunctorAVXTwoCells<true>(newton3, doDeleteSomeParticle, false);
  } else {
    testLJFunctorVSLJFunctorAVXTwoCells<false>(newton3, doDeleteSomeParticle, false);
  }
}

TEST_P(LJFunctorAVXTest, testLJFunctorVSLJFunctorAVXTwoCellsUseUnalignedViews) {
  const auto [mixing, newton3, doDeleteSomeParticle] = GetParam();
  if (mixing) {
    testLJFunctorVSLJFunctorAVXTwoCells<true>(newton3, doDeleteSomeParticle, true);
  } else {
    testLJFunctorVSLJFunctorAVXTwoCells<false>(newton3, doDeleteSomeParticle, true);
  }
}

/**
 * Lambda to generate a readable string out of the parameters of this test.
 */
static auto toString = [](const auto &info) {
  const auto [mixing, newton3, doDeleteSomeParticle] = info.param;
  std::stringstream resStream;
  resStream << (mixing ? "mixingEnabled" : "mixingDisabled") << "_" << (newton3 ? "N3" : "noN3") << "_"
            << (doDeleteSomeParticle ? "withDeletions" : "noDeletions");
  std::string res = resStream.str();
  std::replace(res.begin(), res.end(), '-', '_');
  std::replace(res.begin(), res.end(), '.', '_');
  return res;
};

INSTANTIATE_TEST_SUITE_P(Generated, LJFunctorAVXTest,
                         ::testing::Combine(::testing::Bool(), ::testing::Bool(), ::testing::Bool()), toString);

#endif  // __AVX__