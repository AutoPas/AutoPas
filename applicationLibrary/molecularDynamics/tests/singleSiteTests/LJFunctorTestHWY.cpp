/**
 * @file LJFunctorTestHWY.cpp
 * @author Luis Gall
 * @date 04/23/24
 */

#include "LJFunctorTestHWY.h"

#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/ParticleDefinitions.h"
#include "autopasTools/generators/UniformGenerator.h"
#include "molecularDynamicsLibrary/LJFunctorAVX.h"

template <class SoAType>
bool LJFunctorTestHWY::SoAParticlesEqual(autopas::SoA<SoAType> &soa1, autopas::SoA<SoAType> &soa2) {
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

bool LJFunctorTestHWY::particleEqual(Molecule &p1, Molecule &p2) {
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

bool LJFunctorTestHWY::AoSParticlesEqual(FMCell &cell1, FMCell &cell2) {
  EXPECT_GT(cell1.size(), 0);
  EXPECT_EQ(cell1.size(), cell2.size());

  bool ret = true;

  int counter{0};

  for (size_t i = 0; i < cell1.size(); ++i) {
    if (!particleEqual(cell1._particles[i], cell2._particles[i])) ++counter;
  }

  return counter == 0;

  return ret;
}

template <bool mixing>
void LJFunctorTestHWY::testLJFunctorAVXvsLJFunctorHWYTwoCells(bool newton3, bool doDeleteSomeParticles,
                                                              bool useUnalignedViews, VectorizationPattern pattern) {
  FMCell cell1HWY;
  FMCell cell2HWY;

  size_t numParticles = 23;

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
      cell1HWY, defaultParticle, _lowCorner, {_highCorner[0] / 2, _highCorner[1], _highCorner[2]}, numParticles);
  autopasTools::generators::UniformGenerator::fillWithParticles(
      cell2HWY, defaultParticle, {_highCorner[0] / 2, _lowCorner[1], _lowCorner[2]}, _highCorner, numParticles);

  for (auto &particle : cell1HWY) {
    if (doDeleteSomeParticles) {
      if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
      if (particle.getID() == 11) autopas::internal::markParticleAsDeleted(particle);
      if (particle.getID() == 12) autopas::internal::markParticleAsDeleted(particle);
    }
    if constexpr (mixing) {
      particle.setTypeId(particle.getID() % 5);
    }
  }

  for (auto &particle : cell2HWY) {
    if (doDeleteSomeParticles) {
      if (particle.getID() == 4) autopas::internal::markParticleAsDeleted(particle);
      if (particle.getID() == 20) autopas::internal::markParticleAsDeleted(particle);
      if (particle.getID() == 17) autopas::internal::markParticleAsDeleted(particle);
    }
    if constexpr (mixing) {
      particle.setTypeId(particle.getID() % 5);
    }
  }

  // copy cells
  FMCell cell1NoHWY(cell1HWY);
  FMCell cell2NoHWY(cell2HWY);

  constexpr bool shifting = true;

  auto ljFunctorAVX = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorAVX<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctorAVX<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto ljFunctorHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorHWY<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctorHWY<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();
  ljFunctorHWY.setVecPattern(pattern);

  if constexpr (not mixing) {
    ljFunctorAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    ljFunctorHWY.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  }

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
  // TODO : uncomment before release
  // EXPECT_NEAR(ljFunctorAVX.getPotentialEnergy(), ljFunctorHWY.getUpot(), tolerance) << "global uPot";
  // EXPECT_NEAR(ljFunctorAVX.getVirial(), ljFunctorHWY.getVirial(), tolerance) << "global virial";
}

template <bool mixing>
void LJFunctorTestHWY::testLJFunctorAVXvsLJFunctorHWYOneCell(bool newton3, bool doDeleteSomeParticles,
                                                             bool useUnalignedViews, VectorizationPattern pattern) {
  FMCell cellHWY;

  size_t numParticles = 23;

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
  autopasTools::generators::UniformGenerator::fillWithParticles(cellHWY, defaultParticle, _lowCorner, _highCorner,
                                                                numParticles);

  for (auto &particle : cellHWY) {
    if (doDeleteSomeParticles) {
      if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
      if (particle.getID() == 11) autopas::internal::markParticleAsDeleted(particle);
      if (particle.getID() == 12) autopas::internal::markParticleAsDeleted(particle);
    }
    if constexpr (mixing) {
      particle.setTypeId(particle.getID() % 5);
    }
  }

  // copy cells
  FMCell cellNoHWY(cellHWY);
  constexpr bool shifting = true;

  auto ljFunctorAVX = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorAVX<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctorAVX<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto ljFunctorHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorHWY<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctorHWY<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();
  ljFunctorHWY.setVecPattern(pattern);

  if constexpr (not mixing) {
    ljFunctorAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    ljFunctorHWY.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  }

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
  // TODO : uncomment before release
  // EXPECT_NEAR(ljFunctorAVX.getPotentialEnergy(), ljFunctorHWY.getUpot(), tolerance) << "global uPot";
  // EXPECT_NEAR(ljFunctorAVX.getVirial(), ljFunctorHWY.getVirial(), tolerance) << "global virial";
}

template <bool mixing>
void LJFunctorTestHWY::testLJFunctorAVXvsLJFunctorHWYVerlet(bool newton3, bool doDeleteSomeParticles) {
  using namespace autopas::utils::ArrayMath::literals;

  FMCell cellAVX;

  constexpr size_t numParticles = 23;

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
      if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
      if (particle.getID() == 11) autopas::internal::markParticleAsDeleted(particle);
      if (particle.getID() == 12) autopas::internal::markParticleAsDeleted(particle);
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
  FMCell cellHWY(cellAVX);
  constexpr bool shifting = true;
  constexpr bool calculateGlobals = true;
  auto ljFunctorAVX = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorAVX<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctorAVX<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto ljFunctorHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorHWY<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctorHWY<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  if constexpr (not mixing) {
    ljFunctorAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    ljFunctorHWY.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  }

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
  // TODO : uncomment before release
  // EXPECT_NEAR(ljFunctorAVX.getPotentialEnergy(), ljFunctorHWY.getUpot(), tolerance) << "global uPot";
  // EXPECT_NEAR(ljFunctorAVX.getVirial(), ljFunctorHWY.getVirial(), tolerance) << "global virial";
}

template <bool mixing>
void LJFunctorTestHWY::testLJFunctorAVXvsLJFunctorHWYAoS(bool newton3, bool doDeleteSomeParticles) {
  FMCell cellHWY;

  constexpr size_t numParticles = 23;

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
  autopasTools::generators::UniformGenerator::fillWithParticles(cellHWY, defaultParticle, _lowCorner, _highCorner,
                                                                numParticles);

  for (auto &particle : cellHWY) {
    if (doDeleteSomeParticles) {
      if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
      if (particle.getID() == 11) autopas::internal::markParticleAsDeleted(particle);
      if (particle.getID() == 12) autopas::internal::markParticleAsDeleted(particle);
    }
    if constexpr (mixing) {
      particle.setTypeId(particle.getID() % 5);
    }
  }

  // copy cells
  FMCell cellNoHWY(cellHWY);
  constexpr bool shifting = true;

  auto ljFunctorAVX = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorAVX<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctorAVX<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto ljFunctorHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorHWY<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctorHWY<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  if constexpr (not mixing) {
    ljFunctorAVX.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    ljFunctorHWY.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  }

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
  // TODO : uncomment before release
  // EXPECT_NEAR(ljFunctorHWY.getUpot(), ljFunctorAVX.getPotentialEnergy(), tolerance) << "global uPot";
  // EXPECT_NEAR(ljFunctorHWY.getVirial(), ljFunctorAVX.getVirial(), tolerance) << "global virial";
}
template <bool mixing>
void LJFunctorTestHWY::testPatternSelection(bool newton3) {
  constexpr size_t benchmarkSize = autopas::AutoTuner::_benchmarkSize;
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

  auto ljFunctorHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorHWY<Molecule, false, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctorHWY<Molecule, false, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();
  if constexpr (not mixing) {
    ljFunctorHWY.setParticleProperties(_epsilon, _sigma);
  }
  ljFunctorHWY.setVecPattern(mdLib::VectorizationPattern::p2xVecDiv2);
  std::array<autopas::VectorizationPatternOption::Value, benchmarkSize * benchmarkSize> benchmarkResults;
  for (size_t fcs = 1; fcs <= benchmarkSize; ++fcs) {
    for (size_t scs = 1; scs <= benchmarkSize; ++scs) {
      if (fcs == benchmarkSize && scs == benchmarkSize) {
        benchmarkResults[(fcs - 1) + benchmarkSize * (scs - 1)] = autopas::VectorizationPatternOption::Value::pVecx1;
      } else if (fcs == benchmarkSize) {
        benchmarkResults[(fcs - 1) + benchmarkSize * (scs - 1)] =
            autopas::VectorizationPatternOption::Value::p2xVecDiv2;
      } else if (scs == benchmarkSize) {
        benchmarkResults[(fcs - 1) + benchmarkSize * (scs - 1)] =
            autopas::VectorizationPatternOption::Value::pVecDiv2x2;
      } else {
        benchmarkResults[(fcs - 1) + benchmarkSize * (scs - 1)] = autopas::VectorizationPatternOption::Value::p1xVec;
      }
    }
  }
  ljFunctorHWY.setPatternSelection(&benchmarkResults, &benchmarkResults);

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);

  // not edge case
  FMCell cell1;
  FMCell cell2;
  autopasTools::generators::UniformGenerator::fillWithParticles(cell1, defaultParticle, _lowCorner, _highCorner, 1);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell2, defaultParticle, _lowCorner, _highCorner, 1);
  ljFunctorHWY.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
  ljFunctorHWY.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
  ljFunctorHWY.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, newton3);
  // ASSERT_EQ(ljFunctorHWY.getVecPattern(), mdLib::VectorizationPattern::p1xVec);
  //  right edge
  FMCell cell3;
  FMCell cell4;
  autopasTools::generators::UniformGenerator::fillWithParticles(cell3, defaultParticle, _lowCorner, _highCorner,
                                                                2 * benchmarkSize);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell4, defaultParticle, _lowCorner, _highCorner, 1);
  ljFunctorHWY.SoALoader(cell3, cell3._particleSoABuffer, 0, false);
  ljFunctorHWY.SoALoader(cell4, cell4._particleSoABuffer, 0, false);
  ljFunctorHWY.SoAFunctorPair(cell3._particleSoABuffer, cell4._particleSoABuffer, newton3);
  // ASSERT_EQ(ljFunctorHWY.getVecPattern(), mdLib::VectorizationPattern::p2xVecDiv2);

  // upper edge
  FMCell cell5;
  FMCell cell6;
  autopasTools::generators::UniformGenerator::fillWithParticles(cell5, defaultParticle, _lowCorner, _highCorner, 1);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell6, defaultParticle, _lowCorner, _highCorner,
                                                                2 * benchmarkSize);
  ljFunctorHWY.SoALoader(cell5, cell5._particleSoABuffer, 0, false);
  ljFunctorHWY.SoALoader(cell6, cell6._particleSoABuffer, 0, false);
  ljFunctorHWY.SoAFunctorPair(cell5._particleSoABuffer, cell6._particleSoABuffer, newton3);
  // ASSERT_EQ(ljFunctorHWY.getVecPattern(), mdLib::VectorizationPattern::pVecDiv2x2);
  //  top right corner
  FMCell cell7;
  FMCell cell8;
  autopasTools::generators::UniformGenerator::fillWithParticles(cell7, defaultParticle, _lowCorner, _highCorner,
                                                                2 * benchmarkSize);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell8, defaultParticle, _lowCorner, _highCorner,
                                                                2 * benchmarkSize);
  ljFunctorHWY.SoALoader(cell7, cell7._particleSoABuffer, 0, false);
  ljFunctorHWY.SoALoader(cell8, cell8._particleSoABuffer, 0, false);
  ljFunctorHWY.SoAFunctorPair(cell7._particleSoABuffer, cell8._particleSoABuffer, newton3);
  // ASSERT_EQ(ljFunctorHWY.getVecPattern(), mdLib::VectorizationPattern::pVecx1);
}

TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYAoS) {
  auto [mixing, newton3, doDeleteSomeParticle, _] = GetParam();
  if (mixing) {
    testLJFunctorAVXvsLJFunctorHWYAoS<true>(newton3, doDeleteSomeParticle);
  } else {
    testLJFunctorAVXvsLJFunctorHWYAoS<false>(newton3, doDeleteSomeParticle);
  }
}

TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYVerlet) {
  // different vectorization patterns are currently not supported for Verlet Functor
  auto [mixing, newton3, doDeleteSomeParticle, _] = GetParam();
  if (mixing) {
    testLJFunctorAVXvsLJFunctorHWYVerlet<true>(newton3, doDeleteSomeParticle);
  } else {
    testLJFunctorAVXvsLJFunctorHWYVerlet<false>(newton3, doDeleteSomeParticle);
  }
}

TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYOneCellAlignedAccess) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorAVXvsLJFunctorHWYOneCell<true>(newton3, doDeleteSomeParticle, false, vecPattern);
  } else {
    testLJFunctorAVXvsLJFunctorHWYOneCell<false>(newton3, doDeleteSomeParticle, false, vecPattern);
  }
}

TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYOneCellUseUnalignedViews) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorAVXvsLJFunctorHWYOneCell<true>(newton3, doDeleteSomeParticle, true, vecPattern);
  } else {
    testLJFunctorAVXvsLJFunctorHWYOneCell<false>(newton3, doDeleteSomeParticle, true, vecPattern);
  }
}

TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsAlignedAccess) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorAVXvsLJFunctorHWYTwoCells<true>(newton3, doDeleteSomeParticle, false, vecPattern);
  } else {
    testLJFunctorAVXvsLJFunctorHWYTwoCells<false>(newton3, doDeleteSomeParticle, false, vecPattern);
  }
}

TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsUseUnalignedViews) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorAVXvsLJFunctorHWYTwoCells<true>(newton3, doDeleteSomeParticle, true, vecPattern);
  } else {
    testLJFunctorAVXvsLJFunctorHWYTwoCells<false>(newton3, doDeleteSomeParticle, true, vecPattern);
  }
}

TEST_P(LJFunctorTestHWY, selectionTest) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testPatternSelection<true>(newton3);
  } else {
    testPatternSelection<false>(newton3);
  }
}

std::vector<VectorizationPattern> patterns{VectorizationPattern::p1xVec, VectorizationPattern::p2xVecDiv2,
                                           VectorizationPattern::pVecDiv2x2, VectorizationPattern::pVecx1};

std::map<VectorizationPattern, std::string> patternsToString{{VectorizationPattern::p1xVec, "1xVec"},
                                                             {VectorizationPattern::p2xVecDiv2, "2xVec_2"},
                                                             {VectorizationPattern::pVecDiv2x2, "Vec_2x2"},
                                                             {VectorizationPattern::pVecx1, "Vecx1"}};

static auto toString = [](const auto &info) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = info.param;
  std::stringstream resStream;
  resStream << patternsToString[vecPattern] << (mixing ? "mixing" : "NoMixing") << (newton3 ? "N3" : "NoN3") << "_"
            << (doDeleteSomeParticle ? "withDeletions" : "NoDeletions");
  std::string res = resStream.str();
  std::replace(res.begin(), res.end(), '-', '_');
  std::replace(res.begin(), res.end(), '.', '_');
  return res;
};

INSTANTIATE_TEST_SUITE_P(Generated, LJFunctorTestHWY,
                         ::testing::Combine(::testing::Bool(), ::testing::Bool(), ::testing::Bool(),
                                            ::testing::ValuesIn(patterns)),
                         toString);