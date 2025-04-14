

#include "LJFunctorXSIMDTest.h"

#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/Particle.h"
#include "autopasTools/generators/UniformGenerator.h"
#include "molecularDynamics/molecularDynamicsLibrary/LJFunctorAVX.h"

template <class SoAType>
bool LJFunctorXSIMDTest::SoAParticlesEqual(autopas::SoA<SoAType> &soa1, autopas::SoA<SoAType> &soa2) {
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

bool LJFunctorXSIMDTest::particleEqual(Particle &p1, Particle &p2) {
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

bool LJFunctorXSIMDTest::AoSParticlesEqual(FMCell &cell1, FMCell &cell2) {
  EXPECT_GT(cell1._particles.size(), 0);
  EXPECT_EQ(cell1._particles.size(), cell2._particles.size());

  bool ret = true;
  for (size_t i = 0; i < cell1._particles.size(); ++i) {
    ret = ret and particleEqual(cell1._particles[i], cell2._particles[i]);
  }

  return ret;
}

template <bool mixing>
void LJFunctorXSIMDTest::testLJFunctorVSLJFunctorXSIMDTwoCells(bool newton3, bool doDeleteSomeParticles,
                                                               bool useUnalignedViews) {
  FMCell cell1XSIMD;
  FMCell cell2XSIMD;

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
      cell1XSIMD, defaultParticle, _lowCorner, {_highCorner[0] / 2, _highCorner[1], _highCorner[2]}, numParticles);
  autopasTools::generators::UniformGenerator::fillWithParticles(
      cell2XSIMD, defaultParticle, {_highCorner[0] / 2, _lowCorner[1], _lowCorner[2]}, _highCorner, numParticles);

  for (auto &particle : cell1XSIMD) {
    if (doDeleteSomeParticles) {
      if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
    }
    if constexpr (mixing) {
      particle.setTypeId(particle.getID() % 5);
    }
  }

  for (auto &particle : cell2XSIMD) {
    if (doDeleteSomeParticles) {
      if (particle.getID() == 4) autopas::internal::markParticleAsDeleted(particle);
    }
    if constexpr (mixing) {
      particle.setTypeId(particle.getID() % 5);
    }
  }

  // copy cells
  FMCell cell1NoXSIMD(cell1XSIMD);
  FMCell cell2NoXSIMD(cell2XSIMD);

  constexpr bool shifting = true;

  auto ljFunctorNoXSIMD = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorAVX<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctorAVX<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto ljFunctorXSIMD = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorXSIMD<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctorXSIMD<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  if constexpr (not mixing) {
    ljFunctorNoXSIMD.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    ljFunctorXSIMD.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  }

  ljFunctorXSIMD.initTraversal();
  ljFunctorNoXSIMD.initTraversal();

  ASSERT_TRUE(AoSParticlesEqual(cell1XSIMD, cell1NoXSIMD)) << "Cells 1 not equal after copy initialization.";
  ASSERT_TRUE(AoSParticlesEqual(cell2XSIMD, cell2NoXSIMD)) << "Cells 2 not equal after copy initialization.";

  ljFunctorNoXSIMD.SoALoader(cell1NoXSIMD, cell1NoXSIMD._particleSoABuffer, 0, false);
  ljFunctorNoXSIMD.SoALoader(cell2NoXSIMD, cell2NoXSIMD._particleSoABuffer, 0, false);
  ljFunctorXSIMD.SoALoader(cell1XSIMD, cell1XSIMD._particleSoABuffer, 0, false);
  ljFunctorXSIMD.SoALoader(cell2XSIMD, cell2XSIMD._particleSoABuffer, 0, false);

  ASSERT_TRUE(SoAParticlesEqual(cell1XSIMD._particleSoABuffer, cell1NoXSIMD._particleSoABuffer))
      << "Cells 1 not equal after loading.";
  ASSERT_TRUE(SoAParticlesEqual(cell2XSIMD._particleSoABuffer, cell2NoXSIMD._particleSoABuffer))
      << "Cells 2 not equal after loading.";

  if (useUnalignedViews) {
    ljFunctorNoXSIMD.SoAFunctorPair(cell1NoXSIMD._particleSoABuffer.constructView(1, cell1NoXSIMD.size()),
                                    cell2NoXSIMD._particleSoABuffer.constructView(1, cell2NoXSIMD.size()), newton3);
    ljFunctorXSIMD.SoAFunctorPair(cell1XSIMD._particleSoABuffer.constructView(1, cell1XSIMD.size()),
                                  cell2XSIMD._particleSoABuffer.constructView(1, cell2XSIMD.size()), newton3);
  } else {
    ljFunctorNoXSIMD.SoAFunctorPair(cell1NoXSIMD._particleSoABuffer, cell2NoXSIMD._particleSoABuffer, newton3);
    ljFunctorXSIMD.SoAFunctorPair(cell1XSIMD._particleSoABuffer, cell2XSIMD._particleSoABuffer, newton3);
  }
  ASSERT_TRUE(SoAParticlesEqual(cell1XSIMD._particleSoABuffer, cell1NoXSIMD._particleSoABuffer))
      << "Cells 1 not equal after applying functor.";
  ASSERT_TRUE(SoAParticlesEqual(cell2XSIMD._particleSoABuffer, cell2NoXSIMD._particleSoABuffer))
      << "Cells 2 not equal after applying functor.";

  ljFunctorXSIMD.SoAExtractor(cell1XSIMD, cell1XSIMD._particleSoABuffer, 0);
  ljFunctorXSIMD.SoAExtractor(cell2XSIMD, cell2XSIMD._particleSoABuffer, 0);
  ljFunctorXSIMD.SoAExtractor(cell1NoXSIMD, cell1NoXSIMD._particleSoABuffer, 0);
  ljFunctorXSIMD.SoAExtractor(cell2NoXSIMD, cell2NoXSIMD._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cell1XSIMD, cell1NoXSIMD)) << "Cells 1 not equal after extracting.";
  ASSERT_TRUE(AoSParticlesEqual(cell2XSIMD, cell2NoXSIMD)) << "Cells 2 not equal after extracting.";

  ljFunctorXSIMD.endTraversal(newton3);
  ljFunctorNoXSIMD.endTraversal(newton3);

  double tolerance = 1e-8;
  // EXPECT_NEAR(ljFunctorXSIMD.getUpot(), ljFunctorNoXSIMD.getUpot(), tolerance) << "global uPot";
  // EXPECT_NEAR(ljFunctorXSIMD.getVirial(), ljFunctorNoXSIMD.getVirial(), tolerance) << "global virial";
}

template <bool mixing>
void LJFunctorXSIMDTest::testLJFunctorVSLJFunctorXSIMDOneCell(bool newton3, bool doDeleteSomeParticles,
                                                              bool useUnalignedViews) {
  FMCell cellXSIMD;

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
  autopasTools::generators::UniformGenerator::fillWithParticles(cellXSIMD, defaultParticle, _lowCorner, _highCorner,
                                                                numParticles);

  for (auto &particle : cellXSIMD) {
    if (doDeleteSomeParticles) {
      if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
    }
    if constexpr (mixing) {
      particle.setTypeId(particle.getID() % 5);
    }
  }

  // copy cells
  FMCell cellNoXSIMD(cellXSIMD);
  constexpr bool shifting = true;

  auto ljFunctorNoXSIMD = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorAVX<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctorAVX<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto ljFunctorXSIMD = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorXSIMD<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctorXSIMD<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  if constexpr (not mixing) {
    ljFunctorNoXSIMD.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    ljFunctorXSIMD.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  }

  ASSERT_TRUE(AoSParticlesEqual(cellXSIMD, cellNoXSIMD)) << "Cells not equal after copy initialization.";

  ljFunctorXSIMD.initTraversal();
  ljFunctorNoXSIMD.initTraversal();

  ljFunctorNoXSIMD.SoALoader(cellNoXSIMD, cellNoXSIMD._particleSoABuffer, 0, false);
  ljFunctorXSIMD.SoALoader(cellXSIMD, cellXSIMD._particleSoABuffer, 0, false);

  ASSERT_TRUE(SoAParticlesEqual(cellXSIMD._particleSoABuffer, cellNoXSIMD._particleSoABuffer))
      << "Cells not equal after loading.";

  if (useUnalignedViews) {
    ljFunctorNoXSIMD.SoAFunctorSingle(cellNoXSIMD._particleSoABuffer.constructView(1, cellNoXSIMD.size()), newton3);
    ljFunctorXSIMD.SoAFunctorSingle(cellXSIMD._particleSoABuffer.constructView(1, cellXSIMD.size()), newton3);
  } else {
    ljFunctorNoXSIMD.SoAFunctorSingle(cellNoXSIMD._particleSoABuffer, newton3);
    ljFunctorXSIMD.SoAFunctorSingle(cellXSIMD._particleSoABuffer, newton3);
  }
  ASSERT_TRUE(SoAParticlesEqual(cellXSIMD._particleSoABuffer, cellNoXSIMD._particleSoABuffer))
      << "Cells not equal after applying functor.";

  ljFunctorXSIMD.SoAExtractor(cellXSIMD, cellXSIMD._particleSoABuffer, 0);
  ljFunctorXSIMD.SoAExtractor(cellNoXSIMD, cellNoXSIMD._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cellXSIMD, cellNoXSIMD)) << "Cells 1 not equal after extracting.";

  ljFunctorXSIMD.endTraversal(newton3);
  ljFunctorNoXSIMD.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorXSIMD.getUpot(), ljFunctorNoXSIMD.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorXSIMD.getVirial(), ljFunctorNoXSIMD.getVirial(), tolerance) << "global virial";
}

template <bool mixing>
void LJFunctorXSIMDTest::testLJFunctorVSLJFunctorXSIMDVerlet(bool newton3, bool doDeleteSomeParticles) {
  using namespace autopas::utils::ArrayMath::literals;

  FMCell cellXSIMD;

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
  autopasTools::generators::UniformGenerator::fillWithParticles(cellXSIMD, defaultParticle, _lowCorner, _highCorner,
                                                                numParticles);

  for (auto &particle : cellXSIMD) {
    if (doDeleteSomeParticles) {
      if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
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
      auto dr = cellXSIMD[i].getR() - cellXSIMD[j].getR();
      double dr2 = autopas::utils::ArrayMath::dot(dr, dr);
      if (dr2 <= _interactionLengthSquare) {
        neighborLists[i].push_back(j);
      }
    }
  }

  // copy cells
  FMCell cellNoXSIMD(cellXSIMD);
  constexpr bool shifting = true;
  constexpr bool calculateGlobals = true;

  auto ljFunctorNoXSIMD = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorAVX<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctorAVX<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto ljFunctorXSIMD = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorXSIMD<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctorXSIMD<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  if constexpr (not mixing) {
    ljFunctorNoXSIMD.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    ljFunctorXSIMD.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  }

  ASSERT_TRUE(AoSParticlesEqual(cellXSIMD, cellNoXSIMD)) << "Cells not equal after copy initialization.";

  ljFunctorXSIMD.initTraversal();
  ljFunctorNoXSIMD.initTraversal();

  ljFunctorNoXSIMD.SoALoader(cellNoXSIMD, cellNoXSIMD._particleSoABuffer, 0, false);
  ljFunctorXSIMD.SoALoader(cellXSIMD, cellXSIMD._particleSoABuffer, 0, false);

  ASSERT_TRUE(SoAParticlesEqual(cellXSIMD._particleSoABuffer, cellNoXSIMD._particleSoABuffer))
      << "Cells not equal after loading.";

  for (size_t i = 0; i < numParticles; ++i) {
    ljFunctorNoXSIMD.SoAFunctorVerlet(cellNoXSIMD._particleSoABuffer, i, neighborLists[i], newton3);
    ljFunctorXSIMD.SoAFunctorVerlet(cellXSIMD._particleSoABuffer, i, neighborLists[i], newton3);
  }

  ASSERT_TRUE(SoAParticlesEqual(cellXSIMD._particleSoABuffer, cellNoXSIMD._particleSoABuffer))
      << "Cells not equal after applying functor.";

  ljFunctorXSIMD.SoAExtractor(cellXSIMD, cellXSIMD._particleSoABuffer, 0);
  ljFunctorXSIMD.SoAExtractor(cellNoXSIMD, cellNoXSIMD._particleSoABuffer, 0);

  ASSERT_TRUE(AoSParticlesEqual(cellXSIMD, cellNoXSIMD)) << "Cells not equal after extracting.";

  ljFunctorXSIMD.endTraversal(newton3);
  ljFunctorNoXSIMD.endTraversal(newton3);

  double tolerance = 1e-8;
  // EXPECT_NEAR(ljFunctorXSIMD.getUpot(), ljFunctorNoXSIMD.getUpot(), tolerance) << "global uPot";
  // EXPECT_NEAR(ljFunctorXSIMD.getVirial(), ljFunctorNoXSIMD.getVirial(), tolerance) << "global virial";
}

template <bool mixing>
void LJFunctorXSIMDTest::testLJFunctorVSLJFunctorXSIMDAoS(bool newton3, bool doDeleteSomeParticles) {
  FMCell cellXSIMD;

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
  autopasTools::generators::UniformGenerator::fillWithParticles(cellXSIMD, defaultParticle, _lowCorner, _highCorner,
                                                                numParticles);

  for (auto &particle : cellXSIMD) {
    if (doDeleteSomeParticles) {
      if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
    }
    if constexpr (mixing) {
      particle.setTypeId(particle.getID() % 5);
    }
  }

  // copy cells
  FMCell cellNoXSIMD(cellXSIMD);
  constexpr bool shifting = true;

  auto ljFunctorNoXSIMD = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorAVX<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctorAVX<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  auto ljFunctorXSIMD = [&]() {
    if constexpr (mixing) {
      return mdLib::LJFunctorXSIMD<Molecule, shifting, true, autopas::FunctorN3Modes::Both, true>(_cutoff, PPL);
    } else {
      return mdLib::LJFunctorXSIMD<Molecule, shifting, false, autopas::FunctorN3Modes::Both, true>(_cutoff);
    }
  }();

  if constexpr (not mixing) {
    ljFunctorNoXSIMD.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
    ljFunctorXSIMD.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
  }

  ASSERT_TRUE(AoSParticlesEqual(cellXSIMD, cellNoXSIMD)) << "Cells not equal after copy initialization.";

  ljFunctorXSIMD.initTraversal();
  ljFunctorNoXSIMD.initTraversal();

  for (size_t i = 0; i < numParticles; ++i) {
    for (size_t j = newton3 ? i + 1 : 0; j < numParticles; ++j) {
      if (i == j) {
        continue;
      }
      ljFunctorNoXSIMD.AoSFunctor(cellNoXSIMD[i], cellNoXSIMD[j], newton3);
      ljFunctorXSIMD.AoSFunctor(cellXSIMD[i], cellXSIMD[j], newton3);
    }
  }

  ASSERT_TRUE(AoSParticlesEqual(cellXSIMD, cellNoXSIMD)) << "Cells not equal after applying AoSfunctor.";

  ljFunctorXSIMD.endTraversal(newton3);
  ljFunctorNoXSIMD.endTraversal(newton3);

  double tolerance = 1e-8;
  EXPECT_NEAR(ljFunctorXSIMD.getUpot(), ljFunctorNoXSIMD.getPotentialEnergy(), tolerance) << "global uPot";
  EXPECT_NEAR(ljFunctorXSIMD.getVirial(), ljFunctorNoXSIMD.getVirial(), tolerance) << "global virial";
}

TEST_P(LJFunctorXSIMDTest, testLJFunctorVSLJFunctorXSIMDAoS) {
  auto [mixing, newton3, doDeleteSomeParticle] = GetParam();
  if (mixing) {
    testLJFunctorVSLJFunctorXSIMDAoS<true>(newton3, doDeleteSomeParticle);
  } else {
    testLJFunctorVSLJFunctorXSIMDAoS<false>(newton3, doDeleteSomeParticle);
  }
}

TEST_P(LJFunctorXSIMDTest, testLJFunctorVSLJFunctorXSIMDVerlet) {
  auto [mixing, newton3, doDeleteSomeParticle] = GetParam();
  if (mixing) {
    testLJFunctorVSLJFunctorXSIMDVerlet<true>(newton3, doDeleteSomeParticle);
  } else {
    testLJFunctorVSLJFunctorXSIMDVerlet<false>(newton3, doDeleteSomeParticle);
  }
}

TEST_P(LJFunctorXSIMDTest, testLJFunctorVSLJFunctorXSIMDOneCellAlignedAccess) {
  auto [mixing, newton3, doDeleteSomeParticle] = GetParam();
  if (mixing) {
    testLJFunctorVSLJFunctorXSIMDOneCell<true>(newton3, doDeleteSomeParticle, false);
  } else {
    testLJFunctorVSLJFunctorXSIMDOneCell<false>(newton3, doDeleteSomeParticle, false);
  }
}

TEST_P(LJFunctorXSIMDTest, testLJFunctorVSLJFunctorXSIMDOneCellUseUnalignedViews) {
  auto [mixing, newton3, doDeleteSomeParticle] = GetParam();
  if (mixing) {
    testLJFunctorVSLJFunctorXSIMDOneCell<true>(newton3, doDeleteSomeParticle, true);
  } else {
    testLJFunctorVSLJFunctorXSIMDOneCell<false>(newton3, doDeleteSomeParticle, true);
  }
}

TEST_P(LJFunctorXSIMDTest, testLJFunctorVSLJFunctorXSIMDTwoCellsAlignedAccess) {
  auto [mixing, newton3, doDeleteSomeParticle] = GetParam();
  if (mixing) {
    testLJFunctorVSLJFunctorXSIMDTwoCells<true>(newton3, doDeleteSomeParticle, false);
  } else {
    testLJFunctorVSLJFunctorXSIMDTwoCells<false>(newton3, doDeleteSomeParticle, false);
  }
}

TEST_P(LJFunctorXSIMDTest, testLJFunctorVSLJFunctorXSIMDTwoCellsUseUnalignedViews) {
  auto [mixing, newton3, doDeleteSomeParticle] = GetParam();
  if (mixing) {
    testLJFunctorVSLJFunctorXSIMDTwoCells<true>(newton3, doDeleteSomeParticle, true);
  } else {
    testLJFunctorVSLJFunctorXSIMDTwoCells<false>(newton3, doDeleteSomeParticle, true);
  }
}

/**
 * Lambda to generate a readable string out of the parameters of this test.
 */
static auto toString = [](const auto &info) {
  auto [mixing, newton3, doDeleteSomeParticle] = info.param;
  std::stringstream resStream;
  resStream << (mixing ? "mixing" : "NoMixing") << (newton3 ? "N3" : "noN3") << "_"
            << (doDeleteSomeParticle ? "withDeletions" : "noDeletions");
  std::string res = resStream.str();
  std::replace(res.begin(), res.end(), '-', '_');
  std::replace(res.begin(), res.end(), '.', '_');
  return res;
};

INSTANTIATE_TEST_SUITE_P(Generated, LJFunctorXSIMDTest,
                         ::testing::Combine(::testing::Bool(), ::testing::Bool(), ::testing::Bool()), toString);
