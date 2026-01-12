/**
 * @file ATMFunctorHWYTest.cpp
 * @author D. Martin
 * @date 11/17/25
 */

#include "ATMFunctorHWYTest.h"

#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/ParticleDefinitions.h"
#include "autopasTools/generators/UniformGenerator.h"
#include "molecularDynamicsLibrary/AxilrodTellerMutoFunctor.h"
#include "molecularDynamicsLibrary/AxilrodTellerMutoFunctorHWY.h"

/**
 * Is called before each test case
 */
void ATMFunctorHWYTest::SetUp() {
  bool mixing = std::get<0>(GetParam());
  if (mixing) {
    setupPPL();
  }
}

/**
 * Called from SetUp to initialize the particle properties library with common values.
 */
void ATMFunctorHWYTest::setupPPL() {
  _PPL.addSiteType(0, 1.);
  _PPL.addATMParametersToSite(0, 0.6);
  _PPL.addSiteType(1, 1.5);
  _PPL.addATMParametersToSite(1, 0.8);
  _PPL.addSiteType(2, 2.);
  _PPL.addATMParametersToSite(2, 1.0);
  _PPL.addSiteType(3, 2.5);
  _PPL.addATMParametersToSite(3, 1.2);
  _PPL.addSiteType(4, 3.);
  _PPL.addATMParametersToSite(4, 1.4);
  _PPL.calculateMixingCoefficients();
}

template <class SoAType>
bool ATMFunctorHWYTest::SoAParticlesEqual(const autopas::SoA<SoAType> &soa1, const autopas::SoA<SoAType> &soa2) {
  bool ok = true;

  if (!(soa1.size() > 0)) {
    EXPECT_GT(soa1.size(), 0);
    ok = false;
  }

  if (soa1.size() != soa2.size()) {
    EXPECT_EQ(soa1.size(), soa2.size());
    ok = false;
  }

  const auto *idptr1 = soa1.template begin<Molecule::AttributeNames::id>();
  const auto *idptr2 = soa2.template begin<Molecule::AttributeNames::id>();

  const auto *xptr1 = soa1.template begin<Molecule::AttributeNames::posX>();
  const auto *yptr1 = soa1.template begin<Molecule::AttributeNames::posY>();
  const auto *zptr1 = soa1.template begin<Molecule::AttributeNames::posZ>();
  const auto *xptr2 = soa2.template begin<Molecule::AttributeNames::posX>();
  const auto *yptr2 = soa2.template begin<Molecule::AttributeNames::posY>();
  const auto *zptr2 = soa2.template begin<Molecule::AttributeNames::posZ>();

  const auto *fxptr1 = soa1.template begin<Molecule::AttributeNames::forceX>();
  const auto *fyptr1 = soa1.template begin<Molecule::AttributeNames::forceY>();
  const auto *fzptr1 = soa1.template begin<Molecule::AttributeNames::forceZ>();
  const auto *fxptr2 = soa2.template begin<Molecule::AttributeNames::forceX>();
  const auto *fyptr2 = soa2.template begin<Molecule::AttributeNames::forceY>();
  const auto *fzptr2 = soa2.template begin<Molecule::AttributeNames::forceZ>();

  const size_t n = std::min(soa1.size(), soa2.size());
  for (size_t i = 0; i < n; ++i) {
    if (idptr1[i] != idptr2[i]) {
      EXPECT_EQ(idptr1[i], idptr2[i]);
      ok = false;
    }

    auto near = [&](double a, double b, const char *msg) {
      if (std::abs(a - b) > _maxError) {
        EXPECT_NEAR(a, b, _maxError) << msg << " (i=" << i << ", id=" << idptr1[i] << ")";
        ok = false;
      }
    };

    near(xptr1[i], xptr2[i], "posX");
    near(yptr1[i], yptr2[i], "posY");
    near(zptr1[i], zptr2[i], "posZ");
    near(fxptr1[i], fxptr2[i], "forceX");
    near(fyptr1[i], fyptr2[i], "forceY");
    near(fzptr1[i], fzptr2[i], "forceZ");
  }

  return ok;
}

bool ATMFunctorHWYTest::particleEqual(const Molecule &p1, const Molecule &p2) {
  EXPECT_EQ(p1.getID(), p2.getID());

  EXPECT_NEAR(p1.getR()[0], p2.getR()[0], _maxError) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getR()[1], p2.getR()[1], _maxError) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getR()[2], p2.getR()[2], _maxError) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getF()[0], p2.getF()[0], _maxError) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getF()[1], p2.getF()[1], _maxError) << "for particle pair " << p1.getID();
  EXPECT_NEAR(p1.getF()[2], p2.getF()[2], _maxError) << "for particle pair " << p1.getID();

  // clang-format off
  return not ::testing::Test::HasFailure();
  // clang-format on
}

bool ATMFunctorHWYTest::AoSParticlesEqual(const FMCell &cell1, const FMCell &cell2) {
  EXPECT_GT(cell1.size(), 0);
  EXPECT_EQ(cell1.size(), cell2.size());

  bool ret = true;
  for (size_t i = 0; i < std::min(cell1.size(), cell2.size()); ++i) {
    ret = ret and particleEqual(cell1._particles[i], cell2._particles[i]);
  }

  return ret;
}

template <bool mixing>
void ATMFunctorHWYTest::testATMFunctorVSATMFunctorHWYThreeCells(bool newton3, bool doDeleteSomeParticles) {
  FMCell cell1HWY;
  FMCell cell2HWY;
  FMCell cell3HWY;

  size_t numParticles = 7;

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(
      cell1HWY, defaultParticle, _lowCorner, {_highCorner[0] / 3.0, _highCorner[1], _highCorner[2]}, numParticles);
  autopasTools::generators::UniformGenerator::fillWithParticles(
      cell2HWY, defaultParticle, {_highCorner[0] / 3.0, _lowCorner[1], _lowCorner[2]},
      {(_highCorner[0] / 3.0) * 2.0, _highCorner[1], _highCorner[2]}, numParticles);
  autopasTools::generators::UniformGenerator::fillWithParticles(
      cell3HWY, defaultParticle, {(_highCorner[0] / 3.0) * 2.0, _lowCorner[1], _lowCorner[2]}, _highCorner,
      numParticles);

  for (auto &particle : cell1HWY) {
    if (doDeleteSomeParticles) {
      if (particle.getID() == 3) {
        autopas::internal::markParticleAsDeleted(particle);
      }
    }
    if constexpr (mixing) {
      particle.setTypeId(particle.getID() % 5);
    }
  }
  for (auto &particle : cell2HWY) {
    if (doDeleteSomeParticles) {
      if (particle.getID() == 4) {
        autopas::internal::markParticleAsDeleted(particle);
      }
    }
    if constexpr (mixing) {
      particle.setTypeId(particle.getID() % 5);
    }
    particle.setTypeId(particle.getID() % 5);
  }
  for (auto &particle : cell3HWY) {
    if (doDeleteSomeParticles) {
      if (particle.getID() == 6) {
        autopas::internal::markParticleAsDeleted(particle);
      }
    }
    if constexpr (mixing) {
      particle.setTypeId(particle.getID() % 5);
    }
    particle.setTypeId(particle.getID() % 5);
  }

  // copy cells
  FMCell cell1NoHWY(cell1HWY);
  FMCell cell2NoHWY(cell2HWY);
  FMCell cell3NoHWY(cell3HWY);

  auto atmFunctorNoHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::AxilrodTellerMutoFunctor<Molecule, true, autopas::FunctorN3Modes::Both, true>(_cutoff, _PPL);
    } else {
      return mdLib::AxilrodTellerMutoFunctor<Molecule, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto atmFunctorHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::AxilrodTellerMutoFunctorHWY<Molecule, true, autopas::FunctorN3Modes::Both, true>(_cutoff, _PPL);
    } else {
      return mdLib::AxilrodTellerMutoFunctorHWY<Molecule, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  if constexpr (not mixing) {
    atmFunctorNoHWY.setParticleProperties(_nu);
    atmFunctorHWY.setParticleProperties(_nu);
  }

  atmFunctorHWY.initTraversal();
  atmFunctorNoHWY.initTraversal();

  assertTripletEqual(cell1HWY, cell1NoHWY, cell2HWY, cell2NoHWY, cell3HWY, cell3NoHWY, AoSParticlesEqual,
                     "after copy initialization.");

  atmFunctorNoHWY.SoALoader(cell1NoHWY, cell1NoHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctorNoHWY.SoALoader(cell2NoHWY, cell2NoHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctorNoHWY.SoALoader(cell3NoHWY, cell3NoHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctorHWY.SoALoader(cell1HWY, cell1HWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctorHWY.SoALoader(cell2HWY, cell2HWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctorHWY.SoALoader(cell3HWY, cell3HWY._particleSoABuffer, 0, /*skipSoAResize*/ false);

  assertTripletEqual(
      cell1HWY._particleSoABuffer, cell1NoHWY._particleSoABuffer, cell2HWY._particleSoABuffer,
      cell2NoHWY._particleSoABuffer, cell3HWY._particleSoABuffer, cell3NoHWY._particleSoABuffer,
      [](const auto &a, const auto &b) { return ATMFunctorHWYTest::SoAParticlesEqual(a, b); }, "after loading.");

  atmFunctorNoHWY.SoAFunctorTriple(cell1NoHWY._particleSoABuffer, cell2NoHWY._particleSoABuffer,
                                   cell3NoHWY._particleSoABuffer, newton3);
  atmFunctorHWY.SoAFunctorTriple(cell1HWY._particleSoABuffer, cell2HWY._particleSoABuffer, cell3HWY._particleSoABuffer,
                                 newton3);

  assertTripletEqual(
      cell1HWY._particleSoABuffer, cell1NoHWY._particleSoABuffer, cell2HWY._particleSoABuffer,
      cell2NoHWY._particleSoABuffer, cell3HWY._particleSoABuffer, cell3NoHWY._particleSoABuffer,
      [](const auto &a, const auto &b) { return ATMFunctorHWYTest::SoAParticlesEqual(a, b); },
      "after applying functor.");

  atmFunctorHWY.SoAExtractor(cell1HWY, cell1HWY._particleSoABuffer, 0);
  atmFunctorHWY.SoAExtractor(cell2HWY, cell2HWY._particleSoABuffer, 0);
  atmFunctorHWY.SoAExtractor(cell3HWY, cell3HWY._particleSoABuffer, 0);
  atmFunctorHWY.SoAExtractor(cell1NoHWY, cell1NoHWY._particleSoABuffer, 0);
  atmFunctorHWY.SoAExtractor(cell2NoHWY, cell2NoHWY._particleSoABuffer, 0);
  atmFunctorHWY.SoAExtractor(cell3NoHWY, cell3NoHWY._particleSoABuffer, 0);

  assertTripletEqual(cell1HWY, cell1NoHWY, cell2HWY, cell2NoHWY, cell3HWY, cell3NoHWY, AoSParticlesEqual,
                     "after extracting.");

  atmFunctorHWY.endTraversal(newton3);
  atmFunctorNoHWY.endTraversal(newton3);

  EXPECT_NEAR(atmFunctorHWY.getPotentialEnergy(), atmFunctorNoHWY.getPotentialEnergy(), _maxError) << "global uPot";
  EXPECT_NEAR(atmFunctorHWY.getVirial(), atmFunctorNoHWY.getVirial(), _maxError) << "global virial";
}

template <bool mixing>
void ATMFunctorHWYTest::testATMFunctorVSATMFunctorHWYTwoCells(bool newton3, bool doDeleteSomeParticles) {
  FMCell cell1HWY;
  FMCell cell2HWY;

  size_t numParticles = 7;

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(
      cell1HWY, defaultParticle, _lowCorner, {_highCorner[0] / 2, _highCorner[1], _highCorner[2]}, numParticles);
  autopasTools::generators::UniformGenerator::fillWithParticles(
      cell2HWY, defaultParticle, {_highCorner[0] / 2, _lowCorner[1], _lowCorner[2]}, _highCorner, numParticles);

  for (auto &particle : cell1HWY) {
    if (doDeleteSomeParticles) {
      if (particle.getID() == 3) {
        autopas::internal::markParticleAsDeleted(particle);
      }
    }
    if constexpr (mixing) {
      particle.setTypeId(particle.getID() % 5);
    }
  }
  for (auto &particle : cell2HWY) {
    if (doDeleteSomeParticles) {
      if (particle.getID() == 4) {
        autopas::internal::markParticleAsDeleted(particle);
      }
    }
    if constexpr (mixing) {
      particle.setTypeId(particle.getID() % 5);
    }
    particle.setTypeId(particle.getID() % 5);
  }

  // copy cells
  FMCell cell1NoHWY(cell1HWY);
  FMCell cell2NoHWY(cell2HWY);

  auto atmFunctorNoHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::AxilrodTellerMutoFunctor<Molecule, true, autopas::FunctorN3Modes::Both, true>(_cutoff, _PPL);
    } else {
      return mdLib::AxilrodTellerMutoFunctor<Molecule, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto atmFunctorHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::AxilrodTellerMutoFunctorHWY<Molecule, true, autopas::FunctorN3Modes::Both, true>(_cutoff, _PPL);
    } else {
      return mdLib::AxilrodTellerMutoFunctorHWY<Molecule, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  if constexpr (not mixing) {
    atmFunctorNoHWY.setParticleProperties(_nu);
    atmFunctorHWY.setParticleProperties(_nu);
  }

  atmFunctorHWY.initTraversal();
  atmFunctorNoHWY.initTraversal();

  assertPairEqual(cell1HWY, cell1NoHWY, cell2HWY, cell2NoHWY, AoSParticlesEqual, "after copy initialization.");

  atmFunctorNoHWY.SoALoader(cell1NoHWY, cell1NoHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctorNoHWY.SoALoader(cell2NoHWY, cell2NoHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctorHWY.SoALoader(cell1HWY, cell1HWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctorHWY.SoALoader(cell2HWY, cell2HWY._particleSoABuffer, 0, /*skipSoAResize*/ false);

  assertPairEqual(
      cell1HWY._particleSoABuffer, cell1NoHWY._particleSoABuffer, cell2HWY._particleSoABuffer,
      cell2NoHWY._particleSoABuffer,
      [](const auto &a, const auto &b) { return ATMFunctorHWYTest::SoAParticlesEqual(a, b); }, "after loading.");

  atmFunctorNoHWY.SoAFunctorPair(cell1NoHWY._particleSoABuffer, cell2NoHWY._particleSoABuffer, newton3);
  atmFunctorHWY.SoAFunctorPair(cell1HWY._particleSoABuffer, cell2HWY._particleSoABuffer, newton3);

  assertPairEqual(
      cell1HWY._particleSoABuffer, cell1NoHWY._particleSoABuffer, cell2HWY._particleSoABuffer,
      cell2NoHWY._particleSoABuffer,
      [](const auto &a, const auto &b) { return ATMFunctorHWYTest::SoAParticlesEqual(a, b); },
      "after applying functor.");

  atmFunctorHWY.SoAExtractor(cell1HWY, cell1HWY._particleSoABuffer, 0);
  atmFunctorHWY.SoAExtractor(cell2HWY, cell2HWY._particleSoABuffer, 0);
  atmFunctorHWY.SoAExtractor(cell1NoHWY, cell1NoHWY._particleSoABuffer, 0);
  atmFunctorHWY.SoAExtractor(cell2NoHWY, cell2NoHWY._particleSoABuffer, 0);

  assertPairEqual(cell1HWY, cell1NoHWY, cell2HWY, cell2NoHWY, AoSParticlesEqual, "after extracting.");

  atmFunctorHWY.endTraversal(newton3);
  atmFunctorNoHWY.endTraversal(newton3);

  EXPECT_NEAR(atmFunctorHWY.getPotentialEnergy(), atmFunctorNoHWY.getPotentialEnergy(), _maxError) << "global uPot";
  EXPECT_NEAR(atmFunctorHWY.getVirial(), atmFunctorNoHWY.getVirial(), _maxError) << "global virial";
}

template <bool mixing>
void ATMFunctorHWYTest::testATMFunctorVSATMFunctorHWYOneCell(bool newton3, bool doDeleteSomeParticles) {
  FMCell cellHWY;

  size_t numParticles = 7;

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cellHWY, defaultParticle, _lowCorner, _highCorner,
                                                                numParticles);

  for (auto &particle : cellHWY) {
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
  FMCell cellNoHWY(cellHWY);

  auto atmFunctorNoHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::AxilrodTellerMutoFunctor<Molecule, true, autopas::FunctorN3Modes::Both, true>(_cutoff, _PPL);
    } else {
      return mdLib::AxilrodTellerMutoFunctor<Molecule, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto atmFunctorHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::AxilrodTellerMutoFunctorHWY<Molecule, true, autopas::FunctorN3Modes::Both, true>(_cutoff, _PPL);
    } else {
      return mdLib::AxilrodTellerMutoFunctorHWY<Molecule, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  if constexpr (not mixing) {
    atmFunctorNoHWY.setParticleProperties(_nu);
    atmFunctorHWY.setParticleProperties(_nu);
  }

  ASSERT_TRUE(AoSParticlesEqual(cellHWY, cellNoHWY)) << "Cells not equal after copy initialization.";

  atmFunctorHWY.initTraversal();
  atmFunctorNoHWY.initTraversal();

  atmFunctorNoHWY.SoALoader(cellNoHWY, cellNoHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctorHWY.SoALoader(cellHWY, cellHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);

  ASSERT_TRUE(SoAParticlesEqual(cellHWY._particleSoABuffer, cellNoHWY._particleSoABuffer))
      << "Cells not equal after loading.";

  atmFunctorNoHWY.SoAFunctorSingle(cellNoHWY._particleSoABuffer, newton3);
  atmFunctorHWY.SoAFunctorSingle(cellHWY._particleSoABuffer, newton3);

  ASSERT_TRUE(SoAParticlesEqual(cellHWY._particleSoABuffer, cellNoHWY._particleSoABuffer))
      << "Cells not equal after applying functor.";

  atmFunctorHWY.SoAExtractor(cellHWY, cellHWY._particleSoABuffer, 0);
  atmFunctorHWY.SoAExtractor(cellNoHWY, cellNoHWY._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cellHWY, cellNoHWY)) << "Cells 1 not equal after extracting.";

  atmFunctorHWY.endTraversal(newton3);
  atmFunctorNoHWY.endTraversal(newton3);

  EXPECT_NEAR(atmFunctorHWY.getPotentialEnergy(), atmFunctorNoHWY.getPotentialEnergy(), _maxError) << "global uPot";
  EXPECT_NEAR(atmFunctorHWY.getVirial(), atmFunctorNoHWY.getVirial(), _maxError) << "global virial";
}

template <bool mixing>
void ATMFunctorHWYTest::testATMFunctorVSATMFunctorHWYAoS(bool newton3, bool doDeleteSomeParticles) {
  FMCell cellHWY;

  constexpr size_t numParticles = 7;

  Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cellHWY, defaultParticle, _lowCorner, _highCorner,
                                                                numParticles);

  for (auto &particle : cellHWY) {
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
  FMCell cellNoHWY(cellHWY);

  auto atmFunctorNoHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::AxilrodTellerMutoFunctor<Molecule, true, autopas::FunctorN3Modes::Both, true>(_cutoff, _PPL);
    } else {
      return mdLib::AxilrodTellerMutoFunctor<Molecule, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto atmFunctorHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::AxilrodTellerMutoFunctorHWY<Molecule, true, autopas::FunctorN3Modes::Both, true>(_cutoff, _PPL);
    } else {
      return mdLib::AxilrodTellerMutoFunctorHWY<Molecule, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  if constexpr (not mixing) {
    atmFunctorNoHWY.setParticleProperties(_nu);
    atmFunctorHWY.setParticleProperties(_nu);
  }

  ASSERT_TRUE(AoSParticlesEqual(cellHWY, cellNoHWY)) << "Cells not equal after copy initialization.";

  atmFunctorHWY.initTraversal();
  atmFunctorNoHWY.initTraversal();

  for (size_t i = 0; i < numParticles; ++i) {
    for (size_t j = newton3 ? i + 1 : 0; j < numParticles; ++j) {
      if (j == i) continue;
      for (size_t k = newton3 ? j + 1 : 0; k < numParticles; ++k) {
        if (k == i || k == j) continue;

        atmFunctorNoHWY.AoSFunctor(cellNoHWY[i], cellNoHWY[j], cellNoHWY[k], newton3);
        atmFunctorHWY.AoSFunctor(cellHWY[i], cellHWY[j], cellHWY[k], newton3);
      }
    }
  }

  ASSERT_TRUE(AoSParticlesEqual(cellHWY, cellNoHWY)) << "Cells not equal after applying AoSfunctor.";

  atmFunctorHWY.endTraversal(newton3);
  atmFunctorNoHWY.endTraversal(newton3);

  EXPECT_NEAR(atmFunctorHWY.getPotentialEnergy(), atmFunctorNoHWY.getPotentialEnergy(), _maxError) << "global uPot";
  EXPECT_NEAR(atmFunctorHWY.getVirial(), atmFunctorNoHWY.getVirial(), _maxError) << "global virial";
}

TEST_P(ATMFunctorHWYTest, testATMFunctorVSATMFunctorHWYAoS) {
  const auto [mixing, newton3, doDeleteSomeParticle] = GetParam();
  if (mixing) {
    testATMFunctorVSATMFunctorHWYAoS<true>(newton3, doDeleteSomeParticle);
  } else {
    testATMFunctorVSATMFunctorHWYAoS<false>(newton3, doDeleteSomeParticle);
  }
}

TEST_P(ATMFunctorHWYTest, testATMFunctorVSATMFunctorHWYOneCellAlignedAccess) {
  const auto [mixing, newton3, doDeleteSomeParticle] = GetParam();
  if (mixing) {
    testATMFunctorVSATMFunctorHWYOneCell<true>(newton3, doDeleteSomeParticle);
  } else {
    testATMFunctorVSATMFunctorHWYOneCell<false>(newton3, doDeleteSomeParticle);
  }
}

TEST_P(ATMFunctorHWYTest, testATMFunctorVSATMFunctorHWYTwoCellsAlignedAccess) {
  const auto [mixing, newton3, doDeleteSomeParticle] = GetParam();
  if (mixing) {
    testATMFunctorVSATMFunctorHWYTwoCells<true>(newton3, doDeleteSomeParticle);
  } else {
    testATMFunctorVSATMFunctorHWYTwoCells<false>(newton3, doDeleteSomeParticle);
  }
}

TEST_P(ATMFunctorHWYTest, testATMFunctorVSATMFunctorHWYThreeCellsAlignedAccess) {
  const auto [mixing, newton3, doDeleteSomeParticle] = GetParam();
  if (mixing) {
    testATMFunctorVSATMFunctorHWYThreeCells<true>(newton3, doDeleteSomeParticle);
  } else {
    testATMFunctorVSATMFunctorHWYThreeCells<false>(newton3, doDeleteSomeParticle);
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

INSTANTIATE_TEST_SUITE_P(Generated, ATMFunctorHWYTest,
                         ::testing::Combine(::testing::Bool(), ::testing::Bool(), ::testing::Bool()), toString);