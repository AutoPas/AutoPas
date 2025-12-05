/**
 * @file LJFunctorTestHWY.cpp
 * @author Luis Gall
 * @date 23/04/24
 */

#include "LJFunctorTestHWY.h"

#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/ParticleDefinitions.h"
#include "autopasTools/generators/UniformGenerator.h"
#include "molecularDynamicsLibrary/LJFunctor.h"

template <class SoAType>
bool LJFunctorTestHWY::checkSoAParticlesAreEqual(const autopas::SoA<SoAType> &soa1, const autopas::SoA<SoAType> &soa2) {
  EXPECT_GT(soa1.size(), 0);
  EXPECT_EQ(soa1.size(), soa2.size());

  const unsigned long *const __restrict idptr1 = soa1.template begin<Molecule::AttributeNames::id>();
  const unsigned long *const __restrict idptr2 = soa2.template begin<Molecule::AttributeNames::id>();

  const double *const __restrict xptr1 = soa1.template begin<Molecule::AttributeNames::posX>();
  const double *const __restrict yptr1 = soa1.template begin<Molecule::AttributeNames::posY>();
  const double *const __restrict zptr1 = soa1.template begin<Molecule::AttributeNames::posZ>();
  const double *const __restrict xptr2 = soa2.template begin<Molecule::AttributeNames::posX>();
  const double *const __restrict yptr2 = soa2.template begin<Molecule::AttributeNames::posY>();
  const double *const __restrict zptr2 = soa2.template begin<Molecule::AttributeNames::posZ>();

  const double *const __restrict fxptr1 = soa1.template begin<Molecule::AttributeNames::forceX>();
  const double *const __restrict fyptr1 = soa1.template begin<Molecule::AttributeNames::forceY>();
  const double *const __restrict fzptr1 = soa1.template begin<Molecule::AttributeNames::forceZ>();
  const double *const __restrict fxptr2 = soa2.template begin<Molecule::AttributeNames::forceX>();
  const double *const __restrict fyptr2 = soa2.template begin<Molecule::AttributeNames::forceY>();
  const double *const __restrict fzptr2 = soa2.template begin<Molecule::AttributeNames::forceZ>();

  for (size_t i = 0; i < soa1.size(); ++i) {
    EXPECT_EQ(idptr1[i], idptr2[i]);

    constexpr double tolerance = 2e-8;
    EXPECT_NEAR(xptr1[i], xptr2[i], tolerance) << "for particle pair " << idptr1[i] << " and i=" << i;
    EXPECT_NEAR(yptr1[i], yptr2[i], tolerance) << "for particle pair " << idptr1[i] << " and i=" << i;
    EXPECT_NEAR(zptr1[i], zptr2[i], tolerance) << "for particle pair " << idptr1[i] << " and i=" << i;
    EXPECT_NEAR(fxptr1[i], fxptr2[i], tolerance) << "for particle pair " << idptr1[i] << " and i=" << i;
    EXPECT_NEAR(fyptr1[i], fyptr2[i], tolerance) << "for particle pair " << idptr1[i] << " and i=" << i;
    EXPECT_NEAR(fzptr1[i], fzptr2[i], tolerance) << "for particle pair " << idptr1[i] << " and i=" << i;
  }

  return not HasFailure();
}

bool LJFunctorTestHWY::checkParticlesAreEqual(const Molecule &p1, const Molecule &p2) {
  EXPECT_EQ(p1.getID(), p2.getID());

  double tolerance = 2e-8;
  EXPECT_NEAR(p1.getR()[0], p2.getR()[0], tolerance) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getR()[1], p2.getR()[1], tolerance) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getR()[2], p2.getR()[2], tolerance) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getF()[0], p2.getF()[0], tolerance) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getF()[1], p2.getF()[1], tolerance) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getF()[2], p2.getF()[2], tolerance) << "for particle pair " << p1.getID();

  return not HasFailure();
}

bool LJFunctorTestHWY::checkAoSParticlesAreEqual(const FMCell &cell1, const FMCell &cell2) {
  EXPECT_GT(cell1.size(), 0);
  EXPECT_EQ(cell1.size(), cell2.size());

  bool anyNotEqual = false;

  for (size_t i = 0; i < cell1.size(); ++i) {
    if (!checkParticlesAreEqual(cell1._particles[i], cell2._particles[i])) {
      anyNotEqual = true;
    }
  }

  return !anyNotEqual;
}

template <bool mixing>
void LJFunctorTestHWY::testLJFunctorvsLJFunctorHWYTwoCells(bool newton3, bool doDeleteSomeParticles,
                                                           bool useUnalignedViews, VectorizationPattern pattern) {
  FMCell cell1HWY;
  FMCell cell2HWY;

  const size_t numParticles = 23;

  ParticlePropertiesLibrary PPL{_cutoff};
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

  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(
      cell1HWY, defaultParticle, _lowCorner, {_highCorner[0] / 2, _highCorner[1], _highCorner[2]}, numParticles);
  autopasTools::generators::UniformGenerator::fillWithParticles(
      cell2HWY, defaultParticle, {_highCorner[0] / 2, _lowCorner[1], _lowCorner[2]}, _highCorner, numParticles);

  for (auto &particle : cell1HWY) {
    if (doDeleteSomeParticles) {
      // pick some random particles to be marked as deleted
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
      // pick some random particles to be marked as deleted
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

  auto ljFunctor = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctor<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctor<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto ljFunctorHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorHWY<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, &PPL);
    } else {
      return mdLib::LJFunctorHWY<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();
  ljFunctorHWY.setVecPattern(pattern);

  if constexpr (not mixing) {
    ljFunctor.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    ljFunctorHWY.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  }

  ljFunctorHWY.initTraversal();
  ljFunctor.initTraversal();

  ASSERT_TRUE(checkAoSParticlesAreEqual(cell1HWY, cell1NoHWY)) << "Cells 1 not equal after copy initialization.";
  ASSERT_TRUE(checkAoSParticlesAreEqual(cell2HWY, cell2NoHWY)) << "Cells 2 not equal after copy initialization.";

  ljFunctor.SoALoader(cell1NoHWY, cell1NoHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  ljFunctor.SoALoader(cell2NoHWY, cell2NoHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  ljFunctorHWY.SoALoader(cell1HWY, cell1HWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  ljFunctorHWY.SoALoader(cell2HWY, cell2HWY._particleSoABuffer, 0, /*skipSoAResize*/ false);

  ASSERT_TRUE(checkSoAParticlesAreEqual(cell1HWY._particleSoABuffer, cell1NoHWY._particleSoABuffer))
      << "Cells 1 not equal after loading.";
  ASSERT_TRUE(checkSoAParticlesAreEqual(cell2HWY._particleSoABuffer, cell2NoHWY._particleSoABuffer))
      << "Cells 2 not equal after loading.";

  if (useUnalignedViews) {
    ljFunctor.SoAFunctorPair(cell1NoHWY._particleSoABuffer.constructView(1, cell1NoHWY.size()),
                             cell2NoHWY._particleSoABuffer.constructView(1, cell2NoHWY.size()), newton3);
    ljFunctorHWY.SoAFunctorPair(cell1HWY._particleSoABuffer.constructView(1, cell1HWY.size()),
                                cell2HWY._particleSoABuffer.constructView(1, cell2HWY.size()), newton3);
  } else {
    ljFunctor.SoAFunctorPair(cell1NoHWY._particleSoABuffer, cell2NoHWY._particleSoABuffer, newton3);
    ljFunctorHWY.SoAFunctorPair(cell1HWY._particleSoABuffer, cell2HWY._particleSoABuffer, newton3);
  }
  ASSERT_TRUE(checkSoAParticlesAreEqual(cell1HWY._particleSoABuffer, cell1NoHWY._particleSoABuffer))
      << "Cells 1 not equal after applying functor.";
  ASSERT_TRUE(checkSoAParticlesAreEqual(cell2HWY._particleSoABuffer, cell2NoHWY._particleSoABuffer))
      << "Cells 2 not equal after applying functor.";

  ljFunctorHWY.SoAExtractor(cell1HWY, cell1HWY._particleSoABuffer, 0);
  ljFunctorHWY.SoAExtractor(cell2HWY, cell2HWY._particleSoABuffer, 0);
  ljFunctorHWY.SoAExtractor(cell1NoHWY, cell1NoHWY._particleSoABuffer, 0);
  ljFunctorHWY.SoAExtractor(cell2NoHWY, cell2NoHWY._particleSoABuffer, 0);

  ASSERT_TRUE(checkAoSParticlesAreEqual(cell1HWY, cell1NoHWY)) << "Cells 1 not equal after extracting.";
  ASSERT_TRUE(checkAoSParticlesAreEqual(cell2HWY, cell2NoHWY)) << "Cells 2 not equal after extracting.";

  ljFunctorHWY.endTraversal(newton3);
  ljFunctor.endTraversal(newton3);

  const double tolerance = 1e-7;

  EXPECT_NEAR(ljFunctor.getPotentialEnergy(), ljFunctorHWY.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctor.getVirial(), ljFunctorHWY.getVirial(), tolerance) << "global virial";
}

template <bool mixing>
void LJFunctorTestHWY::testLJFunctorvsLJFunctorHWYOneCell(bool newton3, bool doDeleteSomeParticles,
                                                          bool useUnalignedViews, VectorizationPattern pattern) {
  FMCell cellHWY;

  const size_t numParticles = 23;

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

  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cellHWY, defaultParticle, _lowCorner, _highCorner,
                                                                numParticles);

  for (auto &particle : cellHWY) {
    if (doDeleteSomeParticles) {
      // pick some random particles to be marked as deleted
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

  auto ljFunctor = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctor<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctor<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto ljFunctorHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorHWY<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, &PPL);
    } else {
      return mdLib::LJFunctorHWY<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();
  ljFunctorHWY.setVecPattern(pattern);

  if constexpr (not mixing) {
    ljFunctor.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    ljFunctorHWY.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  }

  ASSERT_TRUE(checkAoSParticlesAreEqual(cellHWY, cellNoHWY)) << "Cells not equal after copy initialization.";

  ljFunctorHWY.initTraversal();
  ljFunctor.initTraversal();

  ljFunctor.SoALoader(cellNoHWY, cellNoHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  ljFunctorHWY.SoALoader(cellHWY, cellHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);

  ASSERT_TRUE(checkSoAParticlesAreEqual(cellHWY._particleSoABuffer, cellNoHWY._particleSoABuffer))
      << "Cells not equal after loading.";

  if (useUnalignedViews) {
    ljFunctor.SoAFunctorSingle(cellNoHWY._particleSoABuffer.constructView(1, cellNoHWY.size()), newton3);
    ljFunctorHWY.SoAFunctorSingle(cellHWY._particleSoABuffer.constructView(1, cellHWY.size()), newton3);
  } else {
    ljFunctor.SoAFunctorSingle(cellNoHWY._particleSoABuffer, newton3);
    ljFunctorHWY.SoAFunctorSingle(cellHWY._particleSoABuffer, newton3);
  }
  ASSERT_TRUE(checkSoAParticlesAreEqual(cellHWY._particleSoABuffer, cellNoHWY._particleSoABuffer))
      << "Cells not equal after applying functor.";

  ljFunctorHWY.SoAExtractor(cellHWY, cellHWY._particleSoABuffer, 0);
  ljFunctorHWY.SoAExtractor(cellNoHWY, cellNoHWY._particleSoABuffer, 0);

  ASSERT_TRUE(checkAoSParticlesAreEqual(cellHWY, cellNoHWY)) << "Cells not equal after extracting.";

  ljFunctorHWY.endTraversal(newton3);
  ljFunctor.endTraversal(newton3);

  const double tolerance = 1e-8;

  EXPECT_NEAR(ljFunctor.getPotentialEnergy(), ljFunctorHWY.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctor.getVirial(), ljFunctorHWY.getVirial(), tolerance) << "global virial";
}

template <bool mixing>
void LJFunctorTestHWY::testLJFunctorvsLJFunctorHWYVerlet(bool newton3, bool doDeleteSomeParticles) {
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

  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cellAVX, defaultParticle, _lowCorner, _highCorner,
                                                                numParticles);

  for (auto &particle : cellAVX) {
    if (doDeleteSomeParticles) {
      // pick some random particles to be marked as deleted
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
  auto ljFunctor = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctor<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctor<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto ljFunctorHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorHWY<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, &PPL);
    } else {
      return mdLib::LJFunctorHWY<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  if constexpr (not mixing) {
    ljFunctor.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    ljFunctorHWY.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  }

  ASSERT_TRUE(checkAoSParticlesAreEqual(cellAVX, cellHWY)) << "Cells not equal after copy initialization.";

  ljFunctor.initTraversal();
  ljFunctorHWY.initTraversal();

  ljFunctorHWY.SoALoader(cellHWY, cellHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  ljFunctor.SoALoader(cellAVX, cellAVX._particleSoABuffer, 0, /*skipSoAResize*/ false);

  ASSERT_TRUE(checkSoAParticlesAreEqual(cellAVX._particleSoABuffer, cellHWY._particleSoABuffer))
      << "Cells not equal after loading.";

  for (size_t i = 0; i < numParticles; ++i) {
    ljFunctorHWY.SoAFunctorVerlet(cellHWY._particleSoABuffer, i, neighborLists[i], newton3);
    ljFunctor.SoAFunctorVerlet(cellAVX._particleSoABuffer, i, neighborLists[i], newton3);
  }

  ASSERT_TRUE(checkSoAParticlesAreEqual(cellAVX._particleSoABuffer, cellHWY._particleSoABuffer))
      << "Cells not equal after applying functor.";

  ljFunctor.SoAExtractor(cellAVX, cellAVX._particleSoABuffer, 0);
  ljFunctor.SoAExtractor(cellHWY, cellHWY._particleSoABuffer, 0);

  ASSERT_TRUE(checkAoSParticlesAreEqual(cellAVX, cellHWY)) << "Cells not equal after extracting.";

  ljFunctor.endTraversal(newton3);
  ljFunctorHWY.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctor.getPotentialEnergy(), ljFunctorHWY.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctor.getVirial(), ljFunctorHWY.getVirial(), tolerance) << "global virial";
}

template <bool mixing>
void LJFunctorTestHWY::testLJFunctorvsLJFunctorHWYAoS(bool newton3, bool doDeleteSomeParticles) {
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

  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cellHWY, defaultParticle, _lowCorner, _highCorner,
                                                                numParticles);

  for (auto &particle : cellHWY) {
    if (doDeleteSomeParticles) {
      // pick some random particles to be marked as deleted
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

  auto ljFunctor = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctor<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctor<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto ljFunctorHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorHWY<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, &PPL);
    } else {
      return mdLib::LJFunctorHWY<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  if constexpr (not mixing) {
    ljFunctor.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    ljFunctorHWY.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  }

  ASSERT_TRUE(checkAoSParticlesAreEqual(cellHWY, cellNoHWY)) << "Cells not equal after copy initialization.";

  ljFunctorHWY.initTraversal();
  ljFunctor.initTraversal();

  for (size_t i = 0; i < numParticles; ++i) {
    for (size_t j = newton3 ? i + 1 : 0; j < numParticles; ++j) {
      if (i == j) {
        continue;
      }
      ljFunctor.AoSFunctor(cellNoHWY[i], cellNoHWY[j], newton3);
      ljFunctorHWY.AoSFunctor(cellHWY[i], cellHWY[j], newton3);
    }
  }

  ASSERT_TRUE(checkAoSParticlesAreEqual(cellHWY, cellNoHWY)) << "Cells not equal after applying AoSfunctor.";

  ljFunctorHWY.endTraversal(newton3);
  ljFunctor.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorHWY.getPotentialEnergy(), ljFunctor.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorHWY.getVirial(), ljFunctor.getVirial(), tolerance) << "global virial";
}

TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYAoS) {
  auto [mixing, newton3, doDeleteSomeParticle, _] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYAoS<true>(newton3, doDeleteSomeParticle);
  } else {
    testLJFunctorvsLJFunctorHWYAoS<false>(newton3, doDeleteSomeParticle);
  }
}

TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYVerlet) {
  // different vectorization patterns are currently not supported for Verlet Functor
  auto [mixing, newton3, doDeleteSomeParticle, _] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYVerlet<true>(newton3, doDeleteSomeParticle);
  } else {
    testLJFunctorvsLJFunctorHWYVerlet<false>(newton3, doDeleteSomeParticle);
  }
}

TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYOneCellAlignedAccess) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYOneCell<true>(newton3, doDeleteSomeParticle, false, vecPattern);
  } else {
    testLJFunctorvsLJFunctorHWYOneCell<false>(newton3, doDeleteSomeParticle, false, vecPattern);
  }
}

TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYOneCellUseUnalignedViews) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYOneCell<true>(newton3, doDeleteSomeParticle, true, vecPattern);
  } else {
    testLJFunctorvsLJFunctorHWYOneCell<false>(newton3, doDeleteSomeParticle, true, vecPattern);
  }
}

TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsAlignedAccess) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYTwoCells<true>(newton3, doDeleteSomeParticle, false, vecPattern);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false>(newton3, doDeleteSomeParticle, false, vecPattern);
  }
}

TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsUseUnalignedViews) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYTwoCells<true>(newton3, doDeleteSomeParticle, true, vecPattern);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false>(newton3, doDeleteSomeParticle, true, vecPattern);
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