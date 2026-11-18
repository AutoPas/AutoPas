/**
 * @file LJFunctorTestHWY.cpp
 * @author Luis Gall
 * @date 23/04/24
 */

#include "LJFunctorTestHWY.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include "autopas/baseFunctors/CellFunctor.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/ParticleDefinitions.h"
#include "autopas/utils/generators/UniformGenerator.h"
#include "molecularDynamicsLibrary/HWYHelpers.h"
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
                                                           VectorizationPattern pattern,
                                                           autopas::SortingDirectionOption direction, bool reversed) {
  std::array<double, 3> cell1Low{}, cell1High{}, cell2Low{}, cell2High{}, sortingDirection{};

  switch (direction) {
    case autopas::SortingDirectionOption::face:
      cell1Low = _lowCorner;
      cell1High = _highCorner;
      cell2Low = {_highCorner[0], _lowCorner[1], _lowCorner[2]};
      cell2High = {2 * _highCorner[0], _highCorner[1], _highCorner[2]};
      sortingDirection = {1.0, 0.0, 0.0};
      break;
    case autopas::SortingDirectionOption::edge:
      cell1Low = _lowCorner;
      cell1High = _highCorner;
      cell2Low = {_highCorner[0], _highCorner[1], _lowCorner[2]};
      cell2High = {2 * _highCorner[0], 2 * _highCorner[1], _highCorner[2]};
      sortingDirection = {1.0 / std::sqrt(2.0), 1.0 / std::sqrt(2.0), 0.0};
      break;
    case autopas::SortingDirectionOption::corner:
      cell1Low = _lowCorner;
      cell1High = _highCorner;
      cell2Low = _highCorner;
      cell2High = {2 * _highCorner[0], 2 * _highCorner[1], 2 * _highCorner[2]};
      sortingDirection = {1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
      break;
  }

  // Swapping cell1 and cell2 leaves the sorting axis unchanged but makes cell1 the one with higher projections.
  if (reversed) {
    std::swap(cell1Low, cell2Low);
    std::swap(cell1High, cell2High);
  }

  const size_t numParticles = 23;
  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);

  FMCell cell1Ref, cell2Ref;
  autopas::generators::UniformGenerator::fillWithParticles(cell1Ref, defaultParticle, cell1Low, cell1High,
                                                           numParticles);
  autopas::generators::UniformGenerator::fillWithParticles(cell2Ref, defaultParticle, cell2Low, cell2High,
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
    cf.setSoASortingThresholds(autopas::SortingThresholdInfoSingle{0});
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

  EXPECT_NEAR(ljFunctor.getPotentialEnergy(), ljFunctorHWY.getPotentialEnergy(), _maxError) << "global uPot";
  EXPECT_NEAR(ljFunctor.getVirial(), ljFunctorHWY.getVirial(), _maxError) << "global virial";
}

template <bool mixing>
void LJFunctorTestHWY::testLJFunctorvsLJFunctorHWYTwoCellsUseUnalignedViews(bool newton3, bool doDeleteSomeParticles,
                                                                            VectorizationPattern pattern) {
  const std::array<double, 3> cell1Low = _lowCorner;
  const std::array<double, 3> cell1High = _highCorner;
  const std::array<double, 3> cell2Low = {_highCorner[0], _lowCorner[1], _lowCorner[2]};
  const std::array<double, 3> cell2High = {2 * _highCorner[0], _highCorner[1], _highCorner[2]};

  const size_t numParticles = 23;
  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);

  FMCell cell1Ref, cell2Ref;
  autopas::generators::UniformGenerator::fillWithParticles(cell1Ref, defaultParticle, cell1Low, cell1High,
                                                           numParticles);
  autopas::generators::UniformGenerator::fillWithParticles(cell2Ref, defaultParticle, cell2Low, cell2High,
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

  // Both the reference and the HWY functor operate on sub-views starting at index 1 (dropping particle 0), so
  // neither the pointers passed to the kernels nor the trip count are aligned to the SIMD vector width.
  auto refView1 = cell1Ref._particleSoABuffer.constructView(1, cell1Ref.size());
  auto refView2 = cell2Ref._particleSoABuffer.constructView(1, cell2Ref.size());
  ljFunctor.SoAFunctorPair(refView1, refView2, newton3);

  auto hwyView1 = cell1HWY._particleSoABuffer.constructView(1, cell1HWY.size());
  auto hwyView2 = cell2HWY._particleSoABuffer.constructView(1, cell2HWY.size());
  ljFunctorHWY.SoAFunctorPair(hwyView1, hwyView2, newton3);

  EXPECT_TRUE(checkSoAParticlesAreEqual(cell1HWY._particleSoABuffer, cell1Ref._particleSoABuffer))
      << "cell1 mismatch after applying functor.";
  EXPECT_TRUE(checkSoAParticlesAreEqual(cell2HWY._particleSoABuffer, cell2Ref._particleSoABuffer))
      << "cell2 mismatch after applying functor.";

  ljFunctor.endTraversal(newton3);
  ljFunctorHWY.endTraversal(newton3);

  EXPECT_NEAR(ljFunctor.getPotentialEnergy(), ljFunctorHWY.getPotentialEnergy(), _maxError) << "global uPot";
  EXPECT_NEAR(ljFunctor.getVirial(), ljFunctorHWY.getVirial(), _maxError) << "global virial";
}

template <bool mixing>
void LJFunctorTestHWY::testLJFunctorvsLJFunctorHWYOneCell(bool newton3, bool doDeleteSomeParticles,
                                                          bool useUnalignedViews, VectorizationPattern pattern) {
  FMCell cellHWY;

  const size_t numParticles = 23;

  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopas::generators::UniformGenerator::fillWithParticles(cellHWY, defaultParticle, _lowCorner, _highCorner,
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
  autopas::generators::UniformGenerator::fillWithParticles(cellAVX, defaultParticle, _lowCorner, _highCorner,
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
  autopas::generators::UniformGenerator::fillWithParticles(cellHWY, defaultParticle, _lowCorner, _highCorner,
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
    testLJFunctorvsLJFunctorHWYTwoCells<true, false>(newton3, doDeleteSomeParticle, vecPattern,
                                                     autopas::SortingDirectionOption::face, false);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, false>(newton3, doDeleteSomeParticle, vecPattern,
                                                      autopas::SortingDirectionOption::face, false);
  }
}

/**
 * Checks that the HWY SoAFunctorPairSorted matches the autovec reference for face-adjacent cells
 * (sorting axis {1,0,0}).
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsSortedFace) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYTwoCells<true, true>(newton3, doDeleteSomeParticle, vecPattern,
                                                    autopas::SortingDirectionOption::face, false);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, true>(newton3, doDeleteSomeParticle, vecPattern,
                                                     autopas::SortingDirectionOption::face, false);
  }
}

/**
 * Checks that the HWY SoAFunctorPair matches the autovec reference for edge-adjacent cells
 * (sorting axis {1/√2, 1/√2, 0}).
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsEdge) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYTwoCells<true, false>(newton3, doDeleteSomeParticle, vecPattern,
                                                     autopas::SortingDirectionOption::edge, false);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, false>(newton3, doDeleteSomeParticle, vecPattern,
                                                      autopas::SortingDirectionOption::edge, false);
  }
}

/**
 * Checks that the HWY SoAFunctorPairSorted matches the autovec reference for edge-adjacent cells
 * (sorting axis {1/√2, 1/√2, 0}).
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsSortedEdge) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYTwoCells<true, true>(newton3, doDeleteSomeParticle, vecPattern,
                                                    autopas::SortingDirectionOption::edge, false);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, true>(newton3, doDeleteSomeParticle, vecPattern,
                                                     autopas::SortingDirectionOption::edge, false);
  }
}

/**
 * Checks that the HWY SoAFunctorPair matches the autovec reference for corner-adjacent cells
 * (sorting axis {1/√3, 1/√3, 1/√3}).
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsCorner) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYTwoCells<true, false>(newton3, doDeleteSomeParticle, vecPattern,
                                                     autopas::SortingDirectionOption::corner, false);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, false>(newton3, doDeleteSomeParticle, vecPattern,
                                                      autopas::SortingDirectionOption::corner, false);
  }
}

/**
 * Checks that the HWY SoAFunctorPairSorted matches the autovec reference for corner-adjacent cells
 * (sorting axis {1/√3, 1/√3, 1/√3}).
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsSortedCorner) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYTwoCells<true, true>(newton3, doDeleteSomeParticle, vecPattern,
                                                    autopas::SortingDirectionOption::corner, false);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, true>(newton3, doDeleteSomeParticle, vecPattern,
                                                     autopas::SortingDirectionOption::corner, false);
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
                                                     autopas::SortingDirectionOption::face, true);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, false>(newton3, doDeleteSomeParticle, vecPattern,
                                                      autopas::SortingDirectionOption::face, true);
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
                                                    autopas::SortingDirectionOption::face, true);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, true>(newton3, doDeleteSomeParticle, vecPattern,
                                                     autopas::SortingDirectionOption::face, true);
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
                                                     autopas::SortingDirectionOption::edge, true);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, false>(newton3, doDeleteSomeParticle, vecPattern,
                                                      autopas::SortingDirectionOption::edge, true);
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
                                                    autopas::SortingDirectionOption::edge, true);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, true>(newton3, doDeleteSomeParticle, vecPattern,
                                                     autopas::SortingDirectionOption::edge, true);
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
                                                     autopas::SortingDirectionOption::corner, true);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, false>(newton3, doDeleteSomeParticle, vecPattern,
                                                      autopas::SortingDirectionOption::corner, true);
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
                                                    autopas::SortingDirectionOption::corner, true);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCells<false, true>(newton3, doDeleteSomeParticle, vecPattern,
                                                     autopas::SortingDirectionOption::corner, true);
  }
}

/**
 * Checks that the HWY SoAFunctorPair matches the autovec reference in the pair-cell (with unaligned views) case,
 * i.e. both cells are passed as sub-views starting at a non-zero, non-vector-aligned offset.
 */
TEST_P(LJFunctorTestHWY, testLJFunctorVSLJFunctorHWYTwoCellsUseUnalignedViews) {
  auto [mixing, newton3, doDeleteSomeParticle, vecPattern] = GetParam();
  if (mixing) {
    testLJFunctorvsLJFunctorHWYTwoCellsUseUnalignedViews<true>(newton3, doDeleteSomeParticle, vecPattern);
  } else {
    testLJFunctorvsLJFunctorHWYTwoCellsUseUnalignedViews<false>(newton3, doDeleteSomeParticle, vecPattern);
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

namespace {

template <bool forceGenericFallback, class D>
void verifyReductionHelpersForTag(D d, double errorTol) {
  const size_t lanes = mdLib::highway::Lanes(d);
  const size_t allocLanes = std::max<size_t>(mdLib::highway::MaxLanes(D{}), lanes);

  auto makeInputs = [&](double base, double step) {
    std::vector<double> values(allocLanes, 0.0);
    double sum = 0.0;
    for (size_t k = 0; k < lanes; ++k) {
      const double sign = (k % 2 == 0) ? 1.0 : -1.0;
      const double val = sign * (base + step * static_cast<double>(k + 1) + 0.0625 * static_cast<double>(k * k));
      values[k] = val;
      sum += val;
    }
    return std::pair{values, sum};
  };

  auto [v0Data, sum0] = makeInputs(1.25, 0.3);
  auto [v1Data, sum1] = makeInputs(2.50, 0.7);
  auto [v2Data, sum2] = makeInputs(0.75, 0.4);
  auto [v3Data, sum3] = makeInputs(3.10, 0.2);
  auto [v4Data, sum4] = makeInputs(1.80, 0.5);
  auto [v5Data, sum5] = makeInputs(2.20, 0.6);

  const auto v0 = mdLib::highway::LoadU(d, v0Data.data());
  const auto v1 = mdLib::highway::LoadU(d, v1Data.data());
  const auto v2 = mdLib::highway::LoadU(d, v2Data.data());
  const auto v3 = mdLib::highway::LoadU(d, v3Data.data());
  const auto v4 = mdLib::highway::LoadU(d, v4Data.data());
  const auto v5 = mdLib::highway::LoadU(d, v5Data.data());

  EXPECT_NEAR(mdLib::hwy_helpers::reduceSumHalving<forceGenericFallback>(d, v0), sum0, errorTol);
  EXPECT_NEAR(mdLib::hwy_helpers::reduceSumHalving<forceGenericFallback>(d, v5), sum5, errorTol);

  double a = 10.5;
  double b = -4.25;
  mdLib::hwy_helpers::reduceAdd2_InterleaveFirst<forceGenericFallback>(d, v0, v1, a, b);
  EXPECT_NEAR(a, 10.5 + sum0, errorTol);
  EXPECT_NEAR(b, -4.25 + sum1, errorTol);
  mdLib::hwy_helpers::reduceSub2_InterleaveFirst<forceGenericFallback>(d, v2, v3, a, b);
  EXPECT_NEAR(a, 10.5 + sum0 - sum2, errorTol);
  EXPECT_NEAR(b, -4.25 + sum1 - sum3, errorTol);

  double x = -1.5;
  double y = 3.25;
  double z = 7.75;
  mdLib::hwy_helpers::reduceAdd3_InterleaveFirst<forceGenericFallback>(d, v0, v1, v2, x, y, z);
  EXPECT_NEAR(x, -1.5 + sum0, errorTol);
  EXPECT_NEAR(y, 3.25 + sum1, errorTol);
  EXPECT_NEAR(z, 7.75 + sum2, errorTol);
  mdLib::hwy_helpers::reduceSub3_InterleaveFirst<forceGenericFallback>(d, v3, v4, v5, x, y, z);
  EXPECT_NEAR(x, -1.5 + sum0 - sum3, errorTol);
  EXPECT_NEAR(y, 3.25 + sum1 - sum4, errorTol);
  EXPECT_NEAR(z, 7.75 + sum2 - sum5, errorTol);

  double g0 = 0.5;
  double g1 = -1.25;
  double g2 = 2.75;
  double g3 = 4.5;
  mdLib::hwy_helpers::reduceAdd4_Interleaved<forceGenericFallback>(d, v0, v1, v2, v3, g0, g1, g2, g3);
  EXPECT_NEAR(g0, 0.5 + sum0, errorTol);
  EXPECT_NEAR(g1, -1.25 + sum1, errorTol);
  EXPECT_NEAR(g2, 2.75 + sum2, errorTol);
  EXPECT_NEAR(g3, 4.5 + sum3, errorTol);

  // Verify Vec2Store only modifies ptr[0..1] and leaves neighboring sentinel elements untouched.
  constexpr size_t kBufSize = 16;
  constexpr size_t kTargetIdx = 4;
  auto initCanaryBuf = [&](double offset) {
    std::array<double, kBufSize> buf{};
    for (size_t idx = 0; idx < kBufSize; ++idx) {
      buf[idx] = offset + static_cast<double>(idx);
    }
    return buf;
  };
  const auto fxInit = initCanaryBuf(100.0);
  const auto fyInit = initCanaryBuf(200.0);
  const auto fzInit = initCanaryBuf(300.0);
  auto fxBuf = fxInit;
  auto fyBuf = fyInit;
  auto fzBuf = fzInit;

  mdLib::hwy_helpers::reduceAdd6Half_Vec2Store<forceGenericFallback>(d, v0, v1, v2, v3, v4, v5, &fxBuf[kTargetIdx],
                                                                     &fyBuf[kTargetIdx], &fzBuf[kTargetIdx]);
  EXPECT_NEAR(fxBuf[kTargetIdx + 0], fxInit[kTargetIdx + 0] + sum0, errorTol);
  EXPECT_NEAR(fxBuf[kTargetIdx + 1], fxInit[kTargetIdx + 1] + sum3, errorTol);
  EXPECT_NEAR(fyBuf[kTargetIdx + 0], fyInit[kTargetIdx + 0] + sum1, errorTol);
  EXPECT_NEAR(fyBuf[kTargetIdx + 1], fyInit[kTargetIdx + 1] + sum4, errorTol);
  EXPECT_NEAR(fzBuf[kTargetIdx + 0], fzInit[kTargetIdx + 0] + sum2, errorTol);
  EXPECT_NEAR(fzBuf[kTargetIdx + 1], fzInit[kTargetIdx + 1] + sum5, errorTol);

  for (size_t idx = 0; idx < kBufSize; ++idx) {
    if (idx != kTargetIdx && idx != kTargetIdx + 1) {
      EXPECT_DOUBLE_EQ(fxBuf[idx], fxInit[idx]) << "Add6Half overwrote fxBuf[" << idx << "]";
      EXPECT_DOUBLE_EQ(fyBuf[idx], fyInit[idx]) << "Add6Half overwrote fyBuf[" << idx << "]";
      EXPECT_DOUBLE_EQ(fzBuf[idx], fzInit[idx]) << "Add6Half overwrote fzBuf[" << idx << "]";
    }
  }

  mdLib::hwy_helpers::reduceSub6Half_Vec2Store<forceGenericFallback>(d, v0, v1, v2, v3, v4, v5, &fxBuf[kTargetIdx],
                                                                     &fyBuf[kTargetIdx], &fzBuf[kTargetIdx]);
  EXPECT_NEAR(fxBuf[kTargetIdx + 0], fxInit[kTargetIdx + 0], errorTol);
  EXPECT_NEAR(fxBuf[kTargetIdx + 1], fxInit[kTargetIdx + 1], errorTol);
  EXPECT_NEAR(fyBuf[kTargetIdx + 0], fyInit[kTargetIdx + 0], errorTol);
  EXPECT_NEAR(fyBuf[kTargetIdx + 1], fyInit[kTargetIdx + 1], errorTol);
  EXPECT_NEAR(fzBuf[kTargetIdx + 0], fzInit[kTargetIdx + 0], errorTol);
  EXPECT_NEAR(fzBuf[kTargetIdx + 1], fzInit[kTargetIdx + 1], errorTol);

  for (size_t idx = 0; idx < kBufSize; ++idx) {
    if (idx != kTargetIdx && idx != kTargetIdx + 1) {
      EXPECT_DOUBLE_EQ(fxBuf[idx], fxInit[idx]) << "Sub6Half overwrote fxBuf[" << idx << "]";
      EXPECT_DOUBLE_EQ(fyBuf[idx], fyInit[idx]) << "Sub6Half overwrote fyBuf[" << idx << "]";
      EXPECT_DOUBLE_EQ(fzBuf[idx], fzInit[idx]) << "Sub6Half overwrote fzBuf[" << idx << "]";
    }
  }
}

template <class D>
void verifyReductionBothModes(D d) {
  constexpr double kTol = 1e-11;
  verifyReductionHelpersForTag<false>(d, kTol);
  verifyReductionHelpersForTag<true>(d, kTol);
}

template <size_t kLanes, size_t kMaxLanes = mdLib::_maxVecLengthDouble>
void verifyFixedTagIfSupported() {
  if constexpr (kMaxLanes >= kLanes) {
    verifyReductionBothModes(mdLib::highway::FixedTag<double, kLanes>{});
  }
}

template <class SoAType>
void expectSoAForcesNear(const autopas::SoA<SoAType> &soa1, const autopas::SoA<SoAType> &soa2, double tol,
                         const std::string &context) {
  ASSERT_EQ(soa1.size(), soa2.size()) << context;
  const double *const __restrict fx1 = soa1.template begin<Molecule::AttributeNames::forceX>();
  const double *const __restrict fy1 = soa1.template begin<Molecule::AttributeNames::forceY>();
  const double *const __restrict fz1 = soa1.template begin<Molecule::AttributeNames::forceZ>();
  const double *const __restrict fx2 = soa2.template begin<Molecule::AttributeNames::forceX>();
  const double *const __restrict fy2 = soa2.template begin<Molecule::AttributeNames::forceY>();
  const double *const __restrict fz2 = soa2.template begin<Molecule::AttributeNames::forceZ>();
  for (size_t i = 0; i < soa1.size(); ++i) {
    EXPECT_NEAR(fx1[i], fx2[i], tol) << context << " i=" << i << " fx";
    EXPECT_NEAR(fy1[i], fy2[i], tol) << context << " i=" << i << " fy";
    EXPECT_NEAR(fz1[i], fz2[i], tol) << context << " i=" << i << " fz";
  }
}

}  // namespace

/**
 * Verifies all custom Highway reduction helpers across native full/half tags, 1-lane tags, and multiple fixed/capped
 * vector widths under both the native fast path and the generic ReduceSum fallback, including strict 2-element
 * store bounds for Vec2Store.
 */
TEST(LJFunctorTestHWYReduction, FixedAndFallbackMatchReference) {
  verifyReductionBothModes(mdLib::tag_double);
  verifyReductionBothModes(mdLib::tag_double_half);
  verifyReductionBothModes(mdLib::highway::FixedTag<double, 1>{});
  verifyReductionBothModes(mdLib::highway::Half<mdLib::highway::CappedTag<double, 2>>{});
  verifyReductionBothModes(mdLib::highway::CappedTag<double, 1>{});
  verifyReductionBothModes(mdLib::highway::CappedTag<double, 2>{});
  verifyReductionBothModes(mdLib::highway::CappedTag<double, 4>{});
  verifyReductionBothModes(mdLib::highway::CappedTag<double, 8>{});

  verifyFixedTagIfSupported<2>();
  verifyFixedTagIfSupported<4>();
  verifyFixedTagIfSupported<8>();
}

/**
 * Verifies pattern allowlist checks and compares LJFunctorHWY against LJFunctor for boundary particle counts
 * spanning 128-bit single-lane halves up to >2048-bit (35 > 32 lanes) across Single, Pair, and Verlet kernels.
 */
TEST(LJFunctorTestHWYBoundaryCounts, MatchesReferenceForEdgeParticleCounts) {
  using namespace autopas::utils::ArrayMath::literals;

  constexpr double cutoff = 6.0;
  constexpr double skin = 2.0;
  constexpr double interactionLengthSq = (cutoff + skin) * (cutoff + skin);
  constexpr double epsilon = 1.0;
  constexpr double sigma = 1.0;
  constexpr double tol = LJFunctorTestHWY::_maxError;

  EXPECT_TRUE(mdLib::LJFunctorHWY<Molecule>(cutoff).isVecPatternAllowed(VectorizationPattern::p1xVec));
  EXPECT_TRUE(mdLib::LJFunctorHWY<Molecule>(cutoff).isVecPatternAllowed(VectorizationPattern::pVecx1));
  EXPECT_EQ(mdLib::LJFunctorHWY<Molecule>(cutoff).isVecPatternAllowed(VectorizationPattern::p2xVecDiv2),
            mdLib::_maxVecLengthDouble >= 2);
  EXPECT_EQ(mdLib::LJFunctorHWY<Molecule>(cutoff).isVecPatternAllowed(VectorizationPattern::pVecDiv2x2),
            mdLib::_maxVecLengthDouble >= 2);

  auto fillDeterministicCell = [](FMCell &cell, size_t count, const std::array<double, 3> &origin, size_t idOffset) {
    cell.clear();
    for (size_t idx = 0; idx < count; ++idx) {
      const double x = origin[0] + 0.5 + 1.1 * static_cast<double>(idx % 4);
      const double y = origin[1] + 0.5 + 1.1 * static_cast<double>((idx / 4) % 4);
      const double z = origin[2] + 0.5 + 1.1 * static_cast<double>(idx / 16);
      const Molecule p({x, y, z}, {0.0, 0.0, 0.0}, idOffset + idx, 0);
      cell.addParticle(p);
    }
  };

  const std::array<size_t, 7> particleCounts{1, 2, 3, 5, 9, 17, 35};
  const std::array<VectorizationPattern, 4> allPatterns{VectorizationPattern::p1xVec, VectorizationPattern::p2xVecDiv2,
                                                        VectorizationPattern::pVecDiv2x2, VectorizationPattern::pVecx1};

  for (const size_t n : particleCounts) {
    {
      FMCell cellRef;
      fillDeterministicCell(cellRef, n, {0.0, 0.0, 0.0}, 0);
      FMCell cellHWY(cellRef);

      mdLib::LJFunctor<Molecule, true, false, autopas::FunctorN3Modes::Both, true> ljRef(cutoff);
      mdLib::LJFunctorHWY<Molecule, true, false, autopas::FunctorN3Modes::Both, true> ljHWY(cutoff);
      ljRef.setParticleProperties(epsilon * 24.0, sigma * sigma);
      ljHWY.setParticleProperties(epsilon * 24.0, sigma * sigma);

      ljRef.initTraversal();
      ljHWY.initTraversal();
      ljRef.SoALoader(cellRef, cellRef._particleSoABuffer, 0, false);
      ljHWY.SoALoader(cellHWY, cellHWY._particleSoABuffer, 0, false);
      ljRef.SoAFunctorSingle(cellRef._particleSoABuffer, true);
      ljHWY.SoAFunctorSingle(cellHWY._particleSoABuffer, true);
      ljRef.endTraversal(true);
      ljHWY.endTraversal(true);

      EXPECT_NEAR(ljHWY.getPotentialEnergy(), ljRef.getPotentialEnergy(), tol) << "Single N=" << n;
      EXPECT_NEAR(ljHWY.getVirial(), ljRef.getVirial(), tol) << "Single N=" << n;
      expectSoAForcesNear(cellHWY._particleSoABuffer, cellRef._particleSoABuffer, tol, "Single N=" + std::to_string(n));
    }

    for (const auto pattern : allPatterns) {
      mdLib::LJFunctorHWY<Molecule, true, false, autopas::FunctorN3Modes::Both, true> patternProbe(cutoff);
      if (not patternProbe.isVecPatternAllowed(pattern)) {
        continue;
      }
      for (const bool newton3 : {false, true}) {
        FMCell cell1Ref, cell2Ref;
        fillDeterministicCell(cell1Ref, n, {0.0, 0.0, 0.0}, 0);
        fillDeterministicCell(cell2Ref, n, {5.0, 0.0, 0.0}, 1000);
        FMCell cell1HWY(cell1Ref), cell2HWY(cell2Ref);

        mdLib::LJFunctor<Molecule, true, false, autopas::FunctorN3Modes::Both, true> ljRef(cutoff);
        mdLib::LJFunctorHWY<Molecule, true, false, autopas::FunctorN3Modes::Both, true> ljHWY(cutoff);
        ljHWY.setVecPattern(pattern);
        ljRef.setParticleProperties(epsilon * 24.0, sigma * sigma);
        ljHWY.setParticleProperties(epsilon * 24.0, sigma * sigma);

        ljRef.initTraversal();
        ljHWY.initTraversal();
        ljRef.SoALoader(cell1Ref, cell1Ref._particleSoABuffer, 0, false);
        ljRef.SoALoader(cell2Ref, cell2Ref._particleSoABuffer, 0, false);
        ljHWY.SoALoader(cell1HWY, cell1HWY._particleSoABuffer, 0, false);
        ljHWY.SoALoader(cell2HWY, cell2HWY._particleSoABuffer, 0, false);
        ljRef.SoAFunctorPair(cell1Ref._particleSoABuffer, cell2Ref._particleSoABuffer, newton3);
        ljHWY.SoAFunctorPair(cell1HWY._particleSoABuffer, cell2HWY._particleSoABuffer, newton3);
        ljRef.endTraversal(newton3);
        ljHWY.endTraversal(newton3);

        const std::string ctx = "Pair N=" + std::to_string(n) +
                                " pattern=" + autopas::VectorizationPatternOption::getOptionNames().at(pattern) +
                                " n3=" + std::to_string(newton3);
        EXPECT_NEAR(ljHWY.getPotentialEnergy(), ljRef.getPotentialEnergy(), tol) << ctx;
        EXPECT_NEAR(ljHWY.getVirial(), ljRef.getVirial(), tol) << ctx;
        expectSoAForcesNear(cell1HWY._particleSoABuffer, cell1Ref._particleSoABuffer, tol, ctx + " cell1");
        expectSoAForcesNear(cell2HWY._particleSoABuffer, cell2Ref._particleSoABuffer, tol, ctx + " cell2");
      }
    }

    for (const bool newton3 : {false, true}) {
      FMCell cellRef;
      fillDeterministicCell(cellRef, n, {0.0, 0.0, 0.0}, 0);
      FMCell cellHWY(cellRef);

      std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> neighborLists(n);
      for (size_t i = 0; i < n; ++i) {
        for (size_t j = newton3 ? i + 1 : 0; j < n; ++j) {
          if (i == j) {
            continue;
          }
          const auto dr = cellRef[i].getR() - cellRef[j].getR();
          if (autopas::utils::ArrayMath::dot(dr, dr) <= interactionLengthSq) {
            neighborLists[i].push_back(j);
          }
        }
      }

      mdLib::LJFunctor<Molecule, true, false, autopas::FunctorN3Modes::Both, true> ljRef(cutoff);
      mdLib::LJFunctorHWY<Molecule, true, false, autopas::FunctorN3Modes::Both, true> ljHWY(cutoff);
      ljRef.setParticleProperties(epsilon * 24.0, sigma * sigma);
      ljHWY.setParticleProperties(epsilon * 24.0, sigma * sigma);

      ljRef.initTraversal();
      ljHWY.initTraversal();
      ljRef.SoALoader(cellRef, cellRef._particleSoABuffer, 0, false);
      ljHWY.SoALoader(cellHWY, cellHWY._particleSoABuffer, 0, false);
      for (size_t i = 0; i < n; ++i) {
        ljRef.SoAFunctorVerlet(cellRef._particleSoABuffer, i, neighborLists[i], newton3);
        ljHWY.SoAFunctorVerlet(cellHWY._particleSoABuffer, i, neighborLists[i], newton3);
      }
      ljRef.endTraversal(newton3);
      ljHWY.endTraversal(newton3);

      const std::string ctx = "Verlet N=" + std::to_string(n) + " n3=" + std::to_string(newton3);
      EXPECT_NEAR(ljHWY.getPotentialEnergy(), ljRef.getPotentialEnergy(), tol) << ctx;
      EXPECT_NEAR(ljHWY.getVirial(), ljRef.getVirial(), tol) << ctx;
      expectSoAForcesNear(cellHWY._particleSoABuffer, cellRef._particleSoABuffer, tol, ctx);
    }
  }
}
