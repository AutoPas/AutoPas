/**
 * @file LJFunctorTestHWY.cpp
 * @author Luis Gall
 * @date 23/04/24
 */

#include "LJFunctorTestHWY.h"

#include "autopas/baseFunctors/CellFunctor.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/ParticleDefinitions.h"
#include "autopasTools/generators/UniformGenerator.h"
#include "molecularDynamicsLibrary/LJFunctor.h"

void LJFunctorTestHWY::SetUp() {
  bool mixing = std::get<0>(GetParam());
  if (mixing) {
    setupPPL();
  }
}

void LJFunctorTestHWY::setupPPL() {
  _PPL.addSiteType(0, 1.0);
  _PPL.addLJParametersToSite(0, 1.0, 1.4);
  _PPL.addSiteType(1, 1.5);
  _PPL.addLJParametersToSite(1, 1.1, 1.3);
  _PPL.addSiteType(2, 2.);
  _PPL.addLJParametersToSite(2, 1.2, 1.2);
  _PPL.addSiteType(3, 2.5);
  _PPL.addLJParametersToSite(3, 1.3, 1.1);
  _PPL.addSiteType(4, 3.);
  _PPL.addLJParametersToSite(4, 1.5, 1.0);
  _PPL.calculateMixingCoefficients();
}

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

    EXPECT_NEAR(xptr1[i], xptr2[i], _maxError) << "for particle pair " << idptr1[i] << " and i=" << i;
    EXPECT_NEAR(yptr1[i], yptr2[i], _maxError) << "for particle pair " << idptr1[i] << " and i=" << i;
    EXPECT_NEAR(zptr1[i], zptr2[i], _maxError) << "for particle pair " << idptr1[i] << " and i=" << i;
    EXPECT_NEAR(fxptr1[i], fxptr2[i], _maxError) << "for particle pair " << idptr1[i] << " and i=" << i;
    EXPECT_NEAR(fyptr1[i], fyptr2[i], _maxError) << "for particle pair " << idptr1[i] << " and i=" << i;
    EXPECT_NEAR(fzptr1[i], fzptr2[i], _maxError) << "for particle pair " << idptr1[i] << " and i=" << i;
  }

  return not HasFailure();
}

bool LJFunctorTestHWY::checkParticlesAreEqual(const Molecule &p1, const Molecule &p2) {
  EXPECT_EQ(p1.getID(), p2.getID());

  EXPECT_NEAR(p1.getR()[0], p2.getR()[0], _maxError) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getR()[1], p2.getR()[1], _maxError) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getR()[2], p2.getR()[2], _maxError) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getF()[0], p2.getF()[0], _maxError) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getF()[1], p2.getF()[1], _maxError) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getF()[2], p2.getF()[2], _maxError) << "for particle pair " << p1.getID();

  return not HasFailure();
}

bool LJFunctorTestHWY::checkAoSParticlesAreEqual(const FMCell &cell1, const FMCell &cell2) {
  EXPECT_GT(cell1.size(), 0);
  EXPECT_EQ(cell1.size(), cell2.size());

  bool allEqual = true;

  for (size_t i = 0; i < std::min(cell1.size(), cell2.size()); ++i) {
    if (not checkParticlesAreEqual(cell1._particles[i], cell2._particles[i])) {
      allEqual = false;
    }
  }

  return allEqual;
}

template <bool mixing, bool sorted>
void LJFunctorTestHWY::testLJFunctorvsLJFunctorHWYTwoCells(bool newton3, bool doDeleteSomeParticles,
                                                           VectorizationPattern pattern, CellLayout layout) {
  std::array<double, 3> cell1Low{}, cell1High{}, cell2Low{}, cell2High{}, sortingDirection{};

  switch (layout) {
    case CellLayout::face:
      cell1Low = _lowCorner;
      cell1High = _highCorner;
      cell2Low = {_highCorner[0], _lowCorner[1], _lowCorner[2]};
      cell2High = {2 * _highCorner[0], _highCorner[1], _highCorner[2]};
      sortingDirection = {1.0, 0.0, 0.0};
      break;
    case CellLayout::edge:
      cell1Low = _lowCorner;
      cell1High = _highCorner;
      cell2Low = {_highCorner[0], _highCorner[1], _lowCorner[2]};
      cell2High = {2 * _highCorner[0], 2 * _highCorner[1], _highCorner[2]};
      sortingDirection = {1.0 / std::sqrt(2.0), 1.0 / std::sqrt(2.0), 0.0};
      break;
    case CellLayout::corner:
      cell1Low = _lowCorner;
      cell1High = _highCorner;
      cell2Low = _highCorner;
      cell2High = {2 * _highCorner[0], 2 * _highCorner[1], 2 * _highCorner[2]};
      sortingDirection = {1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
      break;
    case CellLayout::faceReversed:
      cell1Low = {_highCorner[0], _lowCorner[1], _lowCorner[2]};
      cell1High = {2 * _highCorner[0], _highCorner[1], _highCorner[2]};
      cell2Low = _lowCorner;
      cell2High = _highCorner;
      sortingDirection = {1.0, 0.0, 0.0};
      break;
    case CellLayout::edgeReversed:
      cell1Low = {_highCorner[0], _highCorner[1], _lowCorner[2]};
      cell1High = {2 * _highCorner[0], 2 * _highCorner[1], _highCorner[2]};
      cell2Low = _lowCorner;
      cell2High = _highCorner;
      sortingDirection = {1.0 / std::sqrt(2.0), 1.0 / std::sqrt(2.0), 0.0};
      break;
    case CellLayout::cornerReversed:
      cell1Low = _highCorner;
      cell1High = {2 * _highCorner[0], 2 * _highCorner[1], 2 * _highCorner[2]};
      cell2Low = _lowCorner;
      cell2High = _highCorner;
      sortingDirection = {1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
      break;
  }

  const size_t numParticles = 23;
  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);

  FMCell cell1Ref, cell2Ref;
  autopasTools::generators::UniformGenerator::fillWithParticles(cell1Ref, defaultParticle, cell1Low, cell1High,
                                                                numParticles);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell2Ref, defaultParticle, cell2Low, cell2High,
                                                                numParticles);

  for (auto &particle : cell1Ref) {
    if (doDeleteSomeParticles) {
      if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
      if (particle.getID() == 11) autopas::internal::markParticleAsDeleted(particle);
      if (particle.getID() == 12) autopas::internal::markParticleAsDeleted(particle);
    }
    if constexpr (mixing) particle.setTypeId(particle.getID() % 5);
  }
  for (auto &particle : cell2Ref) {
    if (doDeleteSomeParticles) {
      if (particle.getID() == 4) autopas::internal::markParticleAsDeleted(particle);
      if (particle.getID() == 20) autopas::internal::markParticleAsDeleted(particle);
      if (particle.getID() == 17) autopas::internal::markParticleAsDeleted(particle);
    }
    if constexpr (mixing) particle.setTypeId(particle.getID() % 5);
  }

  FMCell cell1HWY(cell1Ref), cell2HWY(cell2Ref);

  constexpr bool shifting = true;

  auto ljFunctor = [&]() {
    if constexpr (mixing)
      return mdLib::LJFunctor<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, _PPL);
    else
      return mdLib::LJFunctor<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
  }();

  auto ljFunctorHWY = [&]() {
    if constexpr (mixing)
      return mdLib::LJFunctorHWY<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff,
                                                                                                std::ref(_PPL));
    else
      return mdLib::LJFunctorHWY<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
  }();
  ljFunctorHWY.setVecPattern(pattern);

  if constexpr (not mixing) {
    ljFunctor.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    ljFunctorHWY.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  }

  ljFunctor.initTraversal();
  ljFunctorHWY.initTraversal();

  ljFunctor.SoALoader(cell1Ref, cell1Ref._particleSoABuffer, 0, /*skipSoAResize*/ false);
  ljFunctor.SoALoader(cell2Ref, cell2Ref._particleSoABuffer, 0, /*skipSoAResize*/ false);
  ljFunctorHWY.SoALoader(cell1HWY, cell1HWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  ljFunctorHWY.SoALoader(cell2HWY, cell2HWY._particleSoABuffer, 0, /*skipSoAResize*/ false);

  ljFunctor.SoAFunctorPair(cell1Ref._particleSoABuffer, cell2Ref._particleSoABuffer, newton3);

  if constexpr (sorted) {
    autopas::internal::CellFunctor<FMCell, decltype(ljFunctorHWY), /*bidirectional=*/false> cf(
        ljFunctorHWY, _cutoff, autopas::DataLayoutOption::soa, newton3);
    cf.setSoASortingThreshold(0);
    cf.processCellPair(cell1HWY, cell2HWY, sortingDirection);
  } else {
    ljFunctorHWY.SoAFunctorPair(cell1HWY._particleSoABuffer, cell2HWY._particleSoABuffer, newton3);
  }

  EXPECT_TRUE(checkSoAParticlesAreEqual(cell1HWY._particleSoABuffer, cell1Ref._particleSoABuffer))
      << "cell1 mismatch after applying functor.";
  EXPECT_TRUE(checkSoAParticlesAreEqual(cell2HWY._particleSoABuffer, cell2Ref._particleSoABuffer))
      << "cell2 mismatch after applying functor.";

  ljFunctor.endTraversal(newton3);
  ljFunctorHWY.endTraversal(newton3);

  constexpr double energyTol = sorted ? 2 * _maxError : _maxError;
  EXPECT_NEAR(ljFunctor.getPotentialEnergy(), ljFunctorHWY.getPotentialEnergy(), energyTol) << "global uPot";
  EXPECT_NEAR(ljFunctor.getVirial(), ljFunctorHWY.getVirial(), energyTol) << "global virial";
}

template <bool mixing>
void LJFunctorTestHWY::testLJFunctorvsLJFunctorHWYOneCell(bool newton3, bool doDeleteSomeParticles,
                                                          bool useUnalignedViews, VectorizationPattern pattern) {
  FMCell cellHWY;

  const size_t numParticles = 23;

  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cellHWY, defaultParticle, _lowCorner, _highCorner,
                                                                numParticles);

  for (auto &particle : cellHWY) {
    if (doDeleteSomeParticles) {
      // pick some arbitrary particles to be marked as deleted
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
      return mdLib::LJFunctor<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, _PPL);
    } else {
      return mdLib::LJFunctor<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto ljFunctorHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorHWY<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff,
                                                                                                std::ref(_PPL));
    } else {
      return mdLib::LJFunctorHWY<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();
  ljFunctorHWY.setVecPattern(pattern);

  if constexpr (not mixing) {
    ljFunctor.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    ljFunctorHWY.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  }

  EXPECT_TRUE(checkAoSParticlesAreEqual(cellHWY, cellNoHWY)) << "Cells not equal after copy initialization.";

  ljFunctorHWY.initTraversal();
  ljFunctor.initTraversal();

  ljFunctor.SoALoader(cellNoHWY, cellNoHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  ljFunctorHWY.SoALoader(cellHWY, cellHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);

  EXPECT_TRUE(checkSoAParticlesAreEqual(cellHWY._particleSoABuffer, cellNoHWY._particleSoABuffer))
      << "Cells not equal after loading.";

  if (useUnalignedViews) {
    ljFunctor.SoAFunctorSingle(cellNoHWY._particleSoABuffer.constructView(1, cellNoHWY.size()), newton3);
    ljFunctorHWY.SoAFunctorSingle(cellHWY._particleSoABuffer.constructView(1, cellHWY.size()), newton3);
  } else {
    ljFunctor.SoAFunctorSingle(cellNoHWY._particleSoABuffer, newton3);
    ljFunctorHWY.SoAFunctorSingle(cellHWY._particleSoABuffer, newton3);
  }
  EXPECT_TRUE(checkSoAParticlesAreEqual(cellHWY._particleSoABuffer, cellNoHWY._particleSoABuffer))
      << "Cells not equal after applying functor.";

  ljFunctorHWY.SoAExtractor(cellHWY, cellHWY._particleSoABuffer, 0);
  ljFunctorHWY.SoAExtractor(cellNoHWY, cellNoHWY._particleSoABuffer, 0);

  EXPECT_TRUE(checkAoSParticlesAreEqual(cellHWY, cellNoHWY)) << "Cells not equal after extracting.";

  ljFunctorHWY.endTraversal(newton3);
  ljFunctor.endTraversal(newton3);

  EXPECT_NEAR(ljFunctor.getPotentialEnergy(), ljFunctorHWY.getPotentialEnergy(), _maxError) << "global uPot";
  EXPECT_NEAR(ljFunctor.getVirial(), ljFunctorHWY.getVirial(), _maxError) << "global virial";
}

template <bool mixing>
void LJFunctorTestHWY::testLJFunctorvsLJFunctorHWYVerlet(bool newton3, bool doDeleteSomeParticles) {
  using namespace autopas::utils::ArrayMath::literals;

  FMCell cellAVX;

  constexpr size_t numParticles = 23;

  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cellAVX, defaultParticle, _lowCorner, _highCorner,
                                                                numParticles);

  for (auto &particle : cellAVX) {
    if (doDeleteSomeParticles) {
      // pick some arbitrary particles to be marked as deleted
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
      return mdLib::LJFunctor<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, _PPL);
    } else {
      return mdLib::LJFunctor<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto ljFunctorHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorHWY<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff,
                                                                                                std::ref(_PPL));
    } else {
      return mdLib::LJFunctorHWY<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  if constexpr (not mixing) {
    ljFunctor.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    ljFunctorHWY.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  }

  EXPECT_TRUE(checkAoSParticlesAreEqual(cellAVX, cellHWY)) << "Cells not equal after copy initialization.";

  ljFunctor.initTraversal();
  ljFunctorHWY.initTraversal();

  ljFunctorHWY.SoALoader(cellHWY, cellHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  ljFunctor.SoALoader(cellAVX, cellAVX._particleSoABuffer, 0, /*skipSoAResize*/ false);

  EXPECT_TRUE(checkSoAParticlesAreEqual(cellAVX._particleSoABuffer, cellHWY._particleSoABuffer))
      << "Cells not equal after loading.";

  for (size_t i = 0; i < numParticles; ++i) {
    ljFunctorHWY.SoAFunctorVerlet(cellHWY._particleSoABuffer, i, neighborLists[i], newton3);
    ljFunctor.SoAFunctorVerlet(cellAVX._particleSoABuffer, i, neighborLists[i], newton3);
  }

  EXPECT_TRUE(checkSoAParticlesAreEqual(cellAVX._particleSoABuffer, cellHWY._particleSoABuffer))
      << "Cells not equal after applying functor.";

  ljFunctor.SoAExtractor(cellAVX, cellAVX._particleSoABuffer, 0);
  ljFunctor.SoAExtractor(cellHWY, cellHWY._particleSoABuffer, 0);

  EXPECT_TRUE(checkAoSParticlesAreEqual(cellAVX, cellHWY)) << "Cells not equal after extracting.";

  ljFunctor.endTraversal(newton3);
  ljFunctorHWY.endTraversal(newton3);

  EXPECT_NEAR(ljFunctor.getPotentialEnergy(), ljFunctorHWY.getPotentialEnergy(), _maxError) << "global uPot";
  EXPECT_NEAR(ljFunctor.getVirial(), ljFunctorHWY.getVirial(), _maxError) << "global virial";
}

template <bool mixing>
void LJFunctorTestHWY::testLJFunctorvsLJFunctorHWYAoS(bool newton3, bool doDeleteSomeParticles) {
  FMCell cellHWY;

  constexpr size_t numParticles = 23;

  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cellHWY, defaultParticle, _lowCorner, _highCorner,
                                                                numParticles);

  for (auto &particle : cellHWY) {
    if (doDeleteSomeParticles) {
      // pick some arbitrary particles to be marked as deleted
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
      return mdLib::LJFunctor<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, _PPL);
    } else {
      return mdLib::LJFunctor<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto ljFunctorHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorHWY<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff,
                                                                                                std::ref(_PPL));
    } else {
      return mdLib::LJFunctorHWY<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  if constexpr (not mixing) {
    ljFunctor.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    ljFunctorHWY.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  }

  EXPECT_TRUE(checkAoSParticlesAreEqual(cellHWY, cellNoHWY)) << "Cells not equal after copy initialization.";

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

  EXPECT_TRUE(checkAoSParticlesAreEqual(cellHWY, cellNoHWY)) << "Cells not equal after applying AoSfunctor.";

  ljFunctorHWY.endTraversal(newton3);
  ljFunctor.endTraversal(newton3);

  EXPECT_NEAR(ljFunctorHWY.getPotentialEnergy(), ljFunctor.getPotentialEnergy(), _maxError) << "global uPot";
  EXPECT_NEAR(ljFunctorHWY.getVirial(), ljFunctor.getVirial(), _maxError) << "global virial";
}

/**
 * Checks that the HWY Functor computes forces that match the AutoVec functor in the AoS case.
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYAoS) {
  auto [mixing, newton3, doDeleteSomeParticle, _] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYAoS<true>(newton3, doDeleteSomeParticle);
  } else {
    testLJFunctorvsLJFunctorHWYAoS<false>(newton3, doDeleteSomeParticle);
  }
}

/**
 * Checks that the HWY Functor computes forces that match the AutoVec functor in the SoA Verlet case.
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYVerlet) {
  // different vectorization patterns are currently not supported for Verlet Functor
  auto [mixing, newton3, doDeleteSomeParticle, _] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYVerlet<true>(newton3, doDeleteSomeParticle);
  } else {
    testLJFunctorvsLJFunctorHWYVerlet<false>(newton3, doDeleteSomeParticle);
  }
}

/**
 * Checks that the HWY Functor computes forces that match the AutoVec functor in the SoA Single (with aligned views)
 * case.
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYOneCellAlignedAccess) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYOneCell<true>(newton3, doDeleteSomeParticle, false, vecPattern);
  } else {
    testLJFunctorvsLJFunctorHWYOneCell<false>(newton3, doDeleteSomeParticle, false, vecPattern);
  }
}

/**
 * Checks that the HWY Functor computes forces that match the AutoVec functor in the SoA Single (with unaligned views)
 * case.
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYOneCellUseUnalignedViews) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYOneCell<true>(newton3, doDeleteSomeParticle, true, vecPattern);
  } else {
    testLJFunctorvsLJFunctorHWYOneCell<false>(newton3, doDeleteSomeParticle, true, vecPattern);
  }
}

/**
 * Checks that the HWY SoAFunctorPair matches the autovec reference for face-adjacent cells
 * (sorting axis {1,0,0}).
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsFace) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYTwoCells<true, false>(newton3, doDeleteSomeParticle, vecPattern, CellLayout::face);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, false>(newton3, doDeleteSomeParticle, vecPattern, CellLayout::face);
  }
}

/**
 * Checks that the HWY SoAFunctorPairSorted matches the autovec reference for face-adjacent cells
 * (sorting axis {1,0,0}).
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsSortedFace) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYTwoCells<true, true>(newton3, doDeleteSomeParticle, vecPattern, CellLayout::face);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, true>(newton3, doDeleteSomeParticle, vecPattern, CellLayout::face);
  }
}

/**
 * Checks that the HWY SoAFunctorPair matches the autovec reference for edge-adjacent cells
 * (sorting axis {1/√2, 1/√2, 0}).
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsEdge) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYTwoCells<true, false>(newton3, doDeleteSomeParticle, vecPattern, CellLayout::edge);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, false>(newton3, doDeleteSomeParticle, vecPattern, CellLayout::edge);
  }
}

/**
 * Checks that the HWY SoAFunctorPairSorted matches the autovec reference for edge-adjacent cells
 * (sorting axis {1/√2, 1/√2, 0}).
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsSortedEdge) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYTwoCells<true, true>(newton3, doDeleteSomeParticle, vecPattern, CellLayout::edge);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, true>(newton3, doDeleteSomeParticle, vecPattern, CellLayout::edge);
  }
}

/**
 * Checks that the HWY SoAFunctorPair matches the autovec reference for corner-adjacent cells
 * (sorting axis {1/√3, 1/√3, 1/√3}).
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsCorner) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYTwoCells<true, false>(newton3, doDeleteSomeParticle, vecPattern, CellLayout::corner);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, false>(newton3, doDeleteSomeParticle, vecPattern, CellLayout::corner);
  }
}

/**
 * Checks that the HWY SoAFunctorPairSorted matches the autovec reference for corner-adjacent cells
 * (sorting axis {1/√3, 1/√3, 1/√3}).
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsSortedCorner) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYTwoCells<true, true>(newton3, doDeleteSomeParticle, vecPattern, CellLayout::corner);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, true>(newton3, doDeleteSomeParticle, vecPattern, CellLayout::corner);
  }
}

/**
 * Checks that the HWY SoAFunctorPair matches the autovec reference when cell2 is to the LEFT of cell1
 * (face adjacency, sorting axis {1,0,0}).
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsFaceReversed) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYTwoCells<true, false>(newton3, doDeleteSomeParticle, vecPattern,
                                                     CellLayout::faceReversed);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, false>(newton3, doDeleteSomeParticle, vecPattern,
                                                      CellLayout::faceReversed);
  }
}

/**
 * Checks that the HWY SoAFunctorPairSorted matches the autovec reference when cell2 is to the LEFT of cell1
 * (face adjacency, sorting axis {1,0,0}). Exercises minIndex > 0 (left-side j pruning).
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsSortedFaceReversed) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYTwoCells<true, true>(newton3, doDeleteSomeParticle, vecPattern,
                                                    CellLayout::faceReversed);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, true>(newton3, doDeleteSomeParticle, vecPattern,
                                                     CellLayout::faceReversed);
  }
}

/**
 * Checks that the HWY SoAFunctorPair matches the autovec reference when cell2 is to the LEFT of cell1
 * (edge adjacency, sorting axis {1/√2, 1/√2, 0}).
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsEdgeReversed) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYTwoCells<true, false>(newton3, doDeleteSomeParticle, vecPattern,
                                                     CellLayout::edgeReversed);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, false>(newton3, doDeleteSomeParticle, vecPattern,
                                                      CellLayout::edgeReversed);
  }
}

/**
 * Checks that the HWY SoAFunctorPairSorted matches the autovec reference when cell2 is to the LEFT of cell1
 * (edge adjacency, sorting axis {1/√2, 1/√2, 0}). Exercises minIndex > 0 (left-side j pruning).
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsSortedEdgeReversed) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYTwoCells<true, true>(newton3, doDeleteSomeParticle, vecPattern,
                                                    CellLayout::edgeReversed);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, true>(newton3, doDeleteSomeParticle, vecPattern,
                                                     CellLayout::edgeReversed);
  }
}

/**
 * Checks that the HWY SoAFunctorPair matches the autovec reference when cell2 is to the LEFT of cell1
 * (corner adjacency, sorting axis {1/√3, 1/√3, 1/√3}).
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsCornerReversed) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYTwoCells<true, false>(newton3, doDeleteSomeParticle, vecPattern,
                                                     CellLayout::cornerReversed);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, false>(newton3, doDeleteSomeParticle, vecPattern,
                                                      CellLayout::cornerReversed);
  }
}

/**
 * Checks that the HWY SoAFunctorPairSorted matches the autovec reference when cell2 is to the LEFT of cell1
 * (corner adjacency, sorting axis {1/√3, 1/√3, 1/√3}). Exercises minIndex > 0 (left-side j pruning).
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsSortedCornerReversed) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYTwoCells<true, true>(newton3, doDeleteSomeParticle, vecPattern,
                                                    CellLayout::cornerReversed);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, true>(newton3, doDeleteSomeParticle, vecPattern,
                                                     CellLayout::cornerReversed);
  }
}

std::vector<VectorizationPattern> patterns{VectorizationPattern::p1xVec, VectorizationPattern::p2xVecDiv2,
                                           VectorizationPattern::pVecDiv2x2, VectorizationPattern::pVecx1};

static auto toString = [](const auto &info) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = info.param;

  const auto &names = autopas::VectorizationPatternOption::getOptionNames();
  const std::string &vecStr = names.at(vecPattern);

  std::stringstream resStream;
  resStream << vecStr << (mixing ? "mixing" : "NoMixing") << (newton3 ? "N3" : "NoN3") << "_"
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
