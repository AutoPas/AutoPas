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

template <class SoAType>
bool ATMFunctorHWYTest::SoAParticlesEqual(autopas::SoA<SoAType> &soa1, autopas::SoA<SoAType> &soa2) {
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
    EXPECT_NEAR(xptr1[i], xptr2[i], tolerance) << "for particle pair " << idptr1[i] << " and i=" << i;
    EXPECT_NEAR(yptr1[i], yptr2[i], tolerance) << "for particle pair " << idptr1[i] << " and i=" << i;
    EXPECT_NEAR(zptr1[i], zptr2[i], tolerance) << "for particle pair " << idptr1[i] << " and i=" << i;
    EXPECT_NEAR(fxptr1[i], fxptr2[i], tolerance) << "for particle pair " << idptr1[i] << " and i=" << i;
    EXPECT_NEAR(fyptr1[i], fyptr2[i], tolerance) << "for particle pair " << idptr1[i] << " and i=" << i;
    EXPECT_NEAR(fzptr1[i], fzptr2[i], tolerance) << "for particle pair " << idptr1[i] << " and i=" << i;
  }
  // clang-format off
  return not ::testing::Test::HasFailure();
  // clang-format on
}

bool ATMFunctorHWYTest::particleEqual(Molecule &p1, Molecule &p2) {
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

bool ATMFunctorHWYTest::AoSParticlesEqual(FMCell &cell1, FMCell &cell2) {
  EXPECT_GT(cell1.size(), 0);
  EXPECT_EQ(cell1.size(), cell2.size());

  bool ret = true;
  for (size_t i = 0; i < cell1.size(); ++i) {
    ret = ret and particleEqual(cell1._particles[i], cell2._particles[i]);
  }

  return ret;
}

template <bool mixing>
void ATMFunctorHWYTest::testATMFunctorVSATMFunctorHWYThreeCells(bool newton3, bool doDeleteSomeParticles,
                                                                bool useUnalignedViews) {
  FMCell cell1HWY;
  FMCell cell2HWY;
  FMCell cell3HWY;

  size_t numParticles = 7;

  ParticlePropertiesLibrary<double, size_t> PPL{_cutoff};
  if constexpr (mixing) {
    PPL.addSiteType(0, 1.);
    PPL.addATMParametersToSite(0, 1.);
    PPL.addSiteType(1, 1.5);
    PPL.addATMParametersToSite(1, 1.);
    PPL.addSiteType(2, 2.);
    PPL.addATMParametersToSite(2, 1.);
    PPL.addSiteType(3, 2.5);
    PPL.addATMParametersToSite(3, 1.);
    PPL.addSiteType(4, 3.);
    PPL.addATMParametersToSite(4, 1.);
    PPL.calculateMixingCoefficients();
  }

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
      return mdLib::AxilrodTellerMutoFunctor<Molecule, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::AxilrodTellerMutoFunctor<Molecule, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto atmFunctorHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::AxilrodTellerMutoFunctorHWY<Molecule, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
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

  ASSERT_TRUE(AoSParticlesEqual(cell1HWY, cell1NoHWY)) << "Cells 1 not equal after copy initialization.";
  ASSERT_TRUE(AoSParticlesEqual(cell2HWY, cell2NoHWY)) << "Cells 2 not equal after copy initialization.";
  ASSERT_TRUE(AoSParticlesEqual(cell3HWY, cell3NoHWY)) << "Cells 3 not equal after copy initialization.";

  atmFunctorNoHWY.SoALoader(cell1NoHWY, cell1NoHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctorNoHWY.SoALoader(cell2NoHWY, cell2NoHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctorNoHWY.SoALoader(cell3NoHWY, cell3NoHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctorHWY.SoALoader(cell1HWY, cell1HWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctorHWY.SoALoader(cell2HWY, cell2HWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctorHWY.SoALoader(cell3HWY, cell3HWY._particleSoABuffer, 0, /*skipSoAResize*/ false);

  ASSERT_TRUE(SoAParticlesEqual(cell1HWY._particleSoABuffer, cell1NoHWY._particleSoABuffer))
      << "Cells 1 not equal after loading.";
  ASSERT_TRUE(SoAParticlesEqual(cell2HWY._particleSoABuffer, cell2NoHWY._particleSoABuffer))
      << "Cells 2 not equal after loading.";
  ASSERT_TRUE(SoAParticlesEqual(cell3HWY._particleSoABuffer, cell3NoHWY._particleSoABuffer))
      << "Cells 3 not equal after loading.";

  if (useUnalignedViews) {
    atmFunctorNoHWY.SoAFunctorTriple(cell1NoHWY._particleSoABuffer.constructView(1, cell1NoHWY.size()),
                                     cell2NoHWY._particleSoABuffer.constructView(1, cell2NoHWY.size()),
                                     cell3NoHWY._particleSoABuffer.constructView(1, cell3NoHWY.size()), newton3);
    atmFunctorHWY.SoAFunctorTriple(cell1HWY._particleSoABuffer.constructView(1, cell1HWY.size()),
                                   cell2HWY._particleSoABuffer.constructView(1, cell2HWY.size()),
                                   cell3HWY._particleSoABuffer.constructView(1, cell3HWY.size()), newton3);
  } else {
    atmFunctorNoHWY.SoAFunctorTriple(cell1NoHWY._particleSoABuffer, cell2NoHWY._particleSoABuffer,
                                     cell3NoHWY._particleSoABuffer, newton3);
    atmFunctorHWY.SoAFunctorTriple(cell1HWY._particleSoABuffer, cell2HWY._particleSoABuffer,
                                   cell3HWY._particleSoABuffer, newton3);
  }
  ASSERT_TRUE(SoAParticlesEqual(cell1HWY._particleSoABuffer, cell1NoHWY._particleSoABuffer))
      << "Cells 1 not equal after applying functor.";
  ASSERT_TRUE(SoAParticlesEqual(cell2HWY._particleSoABuffer, cell2NoHWY._particleSoABuffer))
      << "Cells 2 not equal after applying functor.";
  ASSERT_TRUE(SoAParticlesEqual(cell3HWY._particleSoABuffer, cell3NoHWY._particleSoABuffer))
      << "Cells 3 not equal after applying functor.";

  atmFunctorHWY.SoAExtractor(cell1HWY, cell1HWY._particleSoABuffer, 0);
  atmFunctorHWY.SoAExtractor(cell2HWY, cell2HWY._particleSoABuffer, 0);
  atmFunctorHWY.SoAExtractor(cell3HWY, cell3HWY._particleSoABuffer, 0);
  atmFunctorHWY.SoAExtractor(cell1NoHWY, cell1NoHWY._particleSoABuffer, 0);
  atmFunctorHWY.SoAExtractor(cell2NoHWY, cell2NoHWY._particleSoABuffer, 0);
  atmFunctorHWY.SoAExtractor(cell3NoHWY, cell3NoHWY._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cell1HWY, cell1NoHWY)) << "Cells 1 not equal after extracting.";
  ASSERT_TRUE(AoSParticlesEqual(cell2HWY, cell2NoHWY)) << "Cells 2 not equal after extracting.";
  ASSERT_TRUE(AoSParticlesEqual(cell3HWY, cell3NoHWY)) << "Cells 3 not equal after extracting.";

  atmFunctorHWY.endTraversal(newton3);
  atmFunctorNoHWY.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(atmFunctorHWY.getPotentialEnergy(), atmFunctorNoHWY.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(atmFunctorHWY.getVirial(), atmFunctorNoHWY.getVirial(), tolerance) << "global virial";
}

template <bool mixing>
void ATMFunctorHWYTest::testATMFunctorVSATMFunctorHWYTwoCells(bool newton3, bool doDeleteSomeParticles,
                                                              bool useUnalignedViews) {
  FMCell cell1HWY;
  FMCell cell2HWY;

  size_t numParticles = 7;

  ParticlePropertiesLibrary<double, size_t> PPL{_cutoff};
  if constexpr (mixing) {
    PPL.addSiteType(0, 1.);
    PPL.addATMParametersToSite(0, 1.);
    PPL.addSiteType(1, 1.5);
    PPL.addATMParametersToSite(1, 1.);
    PPL.addSiteType(2, 2.);
    PPL.addATMParametersToSite(2, 1.);
    PPL.addSiteType(3, 2.5);
    PPL.addATMParametersToSite(3, 1.);
    PPL.addSiteType(4, 3.);
    PPL.addATMParametersToSite(4, 1.);
    PPL.calculateMixingCoefficients();
  }

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
      return mdLib::AxilrodTellerMutoFunctor<Molecule, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::AxilrodTellerMutoFunctor<Molecule, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto atmFunctorHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::AxilrodTellerMutoFunctorHWY<Molecule, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
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

  ASSERT_TRUE(AoSParticlesEqual(cell1HWY, cell1NoHWY)) << "Cells 1 not equal after copy initialization.";
  ASSERT_TRUE(AoSParticlesEqual(cell2HWY, cell2NoHWY)) << "Cells 2 not equal after copy initialization.";

  atmFunctorNoHWY.SoALoader(cell1NoHWY, cell1NoHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctorNoHWY.SoALoader(cell2NoHWY, cell2NoHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctorHWY.SoALoader(cell1HWY, cell1HWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctorHWY.SoALoader(cell2HWY, cell2HWY._particleSoABuffer, 0, /*skipSoAResize*/ false);

  ASSERT_TRUE(SoAParticlesEqual(cell1HWY._particleSoABuffer, cell1NoHWY._particleSoABuffer))
      << "Cells 1 not equal after loading.";
  ASSERT_TRUE(SoAParticlesEqual(cell2HWY._particleSoABuffer, cell2NoHWY._particleSoABuffer))
      << "Cells 2 not equal after loading.";

  if (useUnalignedViews) {
    atmFunctorNoHWY.SoAFunctorPair(cell1NoHWY._particleSoABuffer.constructView(1, cell1NoHWY.size()),
                                   cell2NoHWY._particleSoABuffer.constructView(1, cell2NoHWY.size()), newton3);
    atmFunctorHWY.SoAFunctorPair(cell1HWY._particleSoABuffer.constructView(1, cell1HWY.size()),
                                 cell2HWY._particleSoABuffer.constructView(1, cell2HWY.size()), newton3);
  } else {
    atmFunctorNoHWY.SoAFunctorPair(cell1NoHWY._particleSoABuffer, cell2NoHWY._particleSoABuffer, newton3);
    atmFunctorHWY.SoAFunctorPair(cell1HWY._particleSoABuffer, cell2HWY._particleSoABuffer, newton3);
  }
  ASSERT_TRUE(SoAParticlesEqual(cell1HWY._particleSoABuffer, cell1NoHWY._particleSoABuffer))
      << "Cells 1 not equal after applying functor.";
  ASSERT_TRUE(SoAParticlesEqual(cell2HWY._particleSoABuffer, cell2NoHWY._particleSoABuffer))
      << "Cells 2 not equal after applying functor.";

  atmFunctorHWY.SoAExtractor(cell1HWY, cell1HWY._particleSoABuffer, 0);
  atmFunctorHWY.SoAExtractor(cell2HWY, cell2HWY._particleSoABuffer, 0);
  atmFunctorHWY.SoAExtractor(cell1NoHWY, cell1NoHWY._particleSoABuffer, 0);
  atmFunctorHWY.SoAExtractor(cell2NoHWY, cell2NoHWY._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cell1HWY, cell1NoHWY)) << "Cells 1 not equal after extracting.";
  ASSERT_TRUE(AoSParticlesEqual(cell2HWY, cell2NoHWY)) << "Cells 2 not equal after extracting.";

  atmFunctorHWY.endTraversal(newton3);
  atmFunctorNoHWY.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(atmFunctorHWY.getPotentialEnergy(), atmFunctorNoHWY.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(atmFunctorHWY.getVirial(), atmFunctorNoHWY.getVirial(), tolerance) << "global virial";
}

template <bool mixing>
void ATMFunctorHWYTest::testATMFunctorVSATMFunctorHWYOneCell(bool newton3, bool doDeleteSomeParticles,
                                                             bool useUnalignedViews) {
  FMCell cellHWY;

  size_t numParticles = 7;

  ParticlePropertiesLibrary<double, size_t> PPL{_cutoff};
  if constexpr (mixing) {
    PPL.addSiteType(0, 1.);
    PPL.addATMParametersToSite(0, 1.);
    PPL.addSiteType(1, 1.5);
    PPL.addATMParametersToSite(1, 1.);
    PPL.addSiteType(2, 2.);
    PPL.addATMParametersToSite(2, 1.);
    PPL.addSiteType(3, 2.5);
    PPL.addATMParametersToSite(3, 1.);
    PPL.addSiteType(4, 3.);
    PPL.addATMParametersToSite(4, 1.);
    PPL.calculateMixingCoefficients();
  }

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
      return mdLib::AxilrodTellerMutoFunctor<Molecule, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::AxilrodTellerMutoFunctor<Molecule, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto atmFunctorHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::AxilrodTellerMutoFunctorHWY<Molecule, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
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

  if (useUnalignedViews) {
    atmFunctorNoHWY.SoAFunctorSingle(cellNoHWY._particleSoABuffer.constructView(1, cellNoHWY.size()), newton3);
    atmFunctorHWY.SoAFunctorSingle(cellHWY._particleSoABuffer.constructView(1, cellHWY.size()), newton3);
  } else {
    atmFunctorNoHWY.SoAFunctorSingle(cellNoHWY._particleSoABuffer, newton3);
    atmFunctorHWY.SoAFunctorSingle(cellHWY._particleSoABuffer, newton3);
  }
  ASSERT_TRUE(SoAParticlesEqual(cellHWY._particleSoABuffer, cellNoHWY._particleSoABuffer))
      << "Cells not equal after applying functor.";

  atmFunctorHWY.SoAExtractor(cellHWY, cellHWY._particleSoABuffer, 0);
  atmFunctorHWY.SoAExtractor(cellNoHWY, cellNoHWY._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cellHWY, cellNoHWY)) << "Cells 1 not equal after extracting.";

  atmFunctorHWY.endTraversal(newton3);
  atmFunctorNoHWY.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(atmFunctorHWY.getPotentialEnergy(), atmFunctorNoHWY.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(atmFunctorHWY.getVirial(), atmFunctorNoHWY.getVirial(), tolerance) << "global virial";
}

template <bool mixing>
void ATMFunctorHWYTest::testATMFunctorVSATMFunctorHWYAoS(bool newton3, bool doDeleteSomeParticles) {
  FMCell cellHWY;

  constexpr size_t numParticles = 7;

  ParticlePropertiesLibrary<double, size_t> PPL{_cutoff};
  if constexpr (mixing) {
    PPL.addSiteType(0, 1.);
    PPL.addATMParametersToSite(0, 1.);
    PPL.addSiteType(1, 1.5);
    PPL.addATMParametersToSite(1, 1.);
    PPL.addSiteType(2, 2.);
    PPL.addATMParametersToSite(2, 1.);
    PPL.addSiteType(3, 2.5);
    PPL.addATMParametersToSite(3, 1.);
    PPL.addSiteType(4, 3.);
    PPL.addATMParametersToSite(4, 1.);
    PPL.calculateMixingCoefficients();
  }

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
      return mdLib::AxilrodTellerMutoFunctor<Molecule, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::AxilrodTellerMutoFunctor<Molecule, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto atmFunctorHWY = [&]() {
    if constexpr (mixing) {
      return mdLib::AxilrodTellerMutoFunctorHWY<Molecule, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
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

  double tolerance = 1e-8;
  EXPECT_NEAR(atmFunctorHWY.getPotentialEnergy(), atmFunctorNoHWY.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(atmFunctorHWY.getVirial(), atmFunctorNoHWY.getVirial(), tolerance) << "global virial";
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
    testATMFunctorVSATMFunctorHWYOneCell<true>(newton3, doDeleteSomeParticle, false);
  } else {
    testATMFunctorVSATMFunctorHWYOneCell<false>(newton3, doDeleteSomeParticle, false);
  }
}

TEST_P(ATMFunctorHWYTest, testATMFunctorVSATMFunctorHWYOneCellUseUnalignedViews) {
  const auto [mixing, newton3, doDeleteSomeParticle] = GetParam();
  if (mixing) {
    testATMFunctorVSATMFunctorHWYOneCell<true>(newton3, doDeleteSomeParticle, true);
  } else {
    testATMFunctorVSATMFunctorHWYOneCell<false>(newton3, doDeleteSomeParticle, true);
  }
}

TEST_P(ATMFunctorHWYTest, testATMFunctorVSATMFunctorHWYTwoCellsAlignedAccess) {
  const auto [mixing, newton3, doDeleteSomeParticle] = GetParam();
  if (mixing) {
    testATMFunctorVSATMFunctorHWYTwoCells<true>(newton3, doDeleteSomeParticle, false);
  } else {
    testATMFunctorVSATMFunctorHWYTwoCells<false>(newton3, doDeleteSomeParticle, false);
  }
}

TEST_P(ATMFunctorHWYTest, testATMFunctorVSATMFunctorHWYTwoCellsUseUnalignedViews) {
  const auto [mixing, newton3, doDeleteSomeParticle] = GetParam();
  if (mixing) {
    testATMFunctorVSATMFunctorHWYTwoCells<true>(newton3, doDeleteSomeParticle, true);
  } else {
    testATMFunctorVSATMFunctorHWYTwoCells<false>(newton3, doDeleteSomeParticle, true);
  }
}

TEST_P(ATMFunctorHWYTest, testATMFunctorVSATMFunctorHWYThreeCellsAlignedAccess) {
  const auto [mixing, newton3, doDeleteSomeParticle] = GetParam();
  if (mixing) {
    testATMFunctorVSATMFunctorHWYThreeCells<true>(newton3, doDeleteSomeParticle, false);
  } else {
    testATMFunctorVSATMFunctorHWYThreeCells<false>(newton3, doDeleteSomeParticle, false);
  }
}

TEST_P(ATMFunctorHWYTest, testATMFunctorVSATMFunctorHWYThreeCellsUseUnalignedViews) {
  const auto [mixing, newton3, doDeleteSomeParticle] = GetParam();
  if (mixing) {
    testATMFunctorVSATMFunctorHWYThreeCells<true>(newton3, doDeleteSomeParticle, true);
  } else {
    testATMFunctorVSATMFunctorHWYThreeCells<false>(newton3, doDeleteSomeParticle, true);
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