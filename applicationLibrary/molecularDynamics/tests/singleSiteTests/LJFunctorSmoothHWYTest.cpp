/**
* @file LJFunctorSmoothHWYTest.h
* @author Luis Gall, modified by Ivander Alson Tanjaya
* @date 04/23/24, 09/08/24
 */

#include "LJFunctorSmoothHWYTest.h"

#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/Particle.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "molecularDynamicsLibrary/LJFunctorAVX.h"
#include "molecularDynamicsLibrary/LJFunctorSmoothHWY.h"
#include "molecularDynamicsLibrary/LJFunctorSmooth.h"
#include "molecularDynamicsLibrary/LJFunctorSmoothHWYGS.h"

template <class SoAType>
bool LJFunctorSmoothHWYTest::SoAParticlesEqual(autopas::SoA<SoAType> &soa1, autopas::SoA<SoAType> &soa2) {

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

bool LJFunctorSmoothHWYTest::particleEqual(Particle &p1, Particle &p2) {
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

bool LJFunctorSmoothHWYTest::AoSParticlesEqual(FMCell &cell1, FMCell &cell2) {
 EXPECT_GT(cell1.size(), 0);
 EXPECT_EQ(cell1.size(), cell2.size());

 bool ret = true;

 int counter {0};

 for (size_t i = 0; i < cell1.size(); ++i) {

   if (!particleEqual(cell1._particles[i], cell2._particles[i]))
     ++counter;
 }

 return counter == 0;

 return ret;
}

template <VectorizationPattern vecPattern>
void LJFunctorSmoothHWYTest::testLJFunctorAVXvsLJFunctorHWYTwoCells(bool newton3, bool doDeleteSomeParticles, bool useUnalignedViews) {
 FMCell cell1HWY;
 FMCell cell2HWY;

 size_t numParticles = 23;

 ParticlePropertiesLibrary<double, size_t> PPL{_cutoff};
 PPL.addSiteType(0,1.,1.,1.);
 PPL.addSiteType(1,1.5,2.,1.);
 PPL.addSiteType(2,2.,1.,1.);
 PPL.addSiteType(3,2.5,2.,1.);
 PPL.addSiteType(4,3.,1.,1.);
 PPL.calculateMixingCoefficients();

 Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
 autopasTools::generators::RandomGenerator::fillWithParticles(
     cell1HWY, defaultParticle, _lowCorner, {_highCorner[0] / 2, _highCorner[1], _highCorner[2]}, numParticles);
 autopasTools::generators::RandomGenerator::fillWithParticles(
     cell2HWY, defaultParticle, {_highCorner[0] / 2, _lowCorner[1], _lowCorner[2]}, _highCorner, numParticles);

 if (doDeleteSomeParticles) {
   for (auto &particle : cell1HWY) {
     if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
     if (particle.getID() == 11) autopas::internal::markParticleAsDeleted(particle);
     if (particle.getID() == 12) autopas::internal::markParticleAsDeleted(particle);
   }
   for (auto &particle : cell2HWY) {
     if (particle.getID() == 4) autopas::internal::markParticleAsDeleted(particle);
     if (particle.getID() == 20) autopas::internal::markParticleAsDeleted(particle);
     if (particle.getID() == 17) autopas::internal::markParticleAsDeleted(particle);
   }
 }

 // copy cells
 FMCell cell1NoHWY(cell1HWY);
 FMCell cell2NoHWY(cell2HWY);

 constexpr bool shifting = true;
 constexpr bool mixing = false;
 mdLib::LJFunctorSmooth<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorSmooth(_cutoff,_cutoff*0.54321/*, PPL*/);
 ljFunctorSmooth.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
 mdLib::LJFunctorSmoothHWYGS<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true, false, vecPattern> ljFunctorSmoothHWYGS(_cutoff,_cutoff*0.54321/*, PPL*/);
 ljFunctorSmoothHWYGS.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

 ljFunctorSmoothHWYGS.initTraversal();
 ljFunctorSmooth.initTraversal();

 ASSERT_TRUE(AoSParticlesEqual(cell1HWY, cell1NoHWY)) << "Cells 1 not equal after copy initialization.";
 ASSERT_TRUE(AoSParticlesEqual(cell2HWY, cell2NoHWY)) << "Cells 2 not equal after copy initialization.";

 ljFunctorSmooth.SoALoader(cell1NoHWY, cell1NoHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
 ljFunctorSmooth.SoALoader(cell2NoHWY, cell2NoHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
 ljFunctorSmoothHWYGS.SoALoader(cell1HWY, cell1HWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
 ljFunctorSmoothHWYGS.SoALoader(cell2HWY, cell2HWY._particleSoABuffer, 0, /*skipSoAResize*/ false);

 ASSERT_TRUE(SoAParticlesEqual(cell1HWY._particleSoABuffer, cell1NoHWY._particleSoABuffer))
     << "Cells 1 not equal after loading.";
 ASSERT_TRUE(SoAParticlesEqual(cell2HWY._particleSoABuffer, cell2NoHWY._particleSoABuffer))
     << "Cells 2 not equal after loading.";

 if (useUnalignedViews) {
   ljFunctorSmooth.SoAFunctorPair(cell1NoHWY._particleSoABuffer.constructView(1, cell1NoHWY.size()),
                                  cell2NoHWY._particleSoABuffer.constructView(1, cell2NoHWY.size()), newton3);
   ljFunctorSmoothHWYGS.SoAFunctorPair(cell1HWY._particleSoABuffer.constructView(1, cell1HWY.size()),
                                       cell2HWY._particleSoABuffer.constructView(1, cell2HWY.size()), newton3);
 } else {
   ljFunctorSmooth.SoAFunctorPair(cell1NoHWY._particleSoABuffer, cell2NoHWY._particleSoABuffer, newton3);
   ljFunctorSmoothHWYGS.SoAFunctorPair(cell1HWY._particleSoABuffer, cell2HWY._particleSoABuffer, newton3);
 }
 ASSERT_TRUE(SoAParticlesEqual(cell1HWY._particleSoABuffer, cell1NoHWY._particleSoABuffer))
     << "Cells 1 not equal after applying functor.";
 ASSERT_TRUE(SoAParticlesEqual(cell2HWY._particleSoABuffer, cell2NoHWY._particleSoABuffer))
     << "Cells 2 not equal after applying functor.";

 ljFunctorSmoothHWYGS.SoAExtractor(cell1HWY, cell1HWY._particleSoABuffer, 0);
 ljFunctorSmoothHWYGS.SoAExtractor(cell2HWY, cell2HWY._particleSoABuffer, 0);
 ljFunctorSmoothHWYGS.SoAExtractor(cell1NoHWY, cell1NoHWY._particleSoABuffer, 0);
 ljFunctorSmoothHWYGS.SoAExtractor(cell2NoHWY, cell2NoHWY._particleSoABuffer, 0);

 ASSERT_TRUE(AoSParticlesEqual(cell1HWY, cell1NoHWY)) << "Cells 1 not equal after extracting.";
 ASSERT_TRUE(AoSParticlesEqual(cell2HWY, cell2NoHWY)) << "Cells 2 not equal after extracting.";

 ljFunctorSmoothHWYGS.endTraversal(newton3);
 ljFunctorSmooth.endTraversal(newton3);

 double tolerance = 1e-8;
 // TODO : uncomment before release
 // EXPECT_NEAR(ljFunctorAVX.getPotentialEnergy(), ljFunctorHWY.getUpot(), tolerance) << "global uPot";
 // EXPECT_NEAR(ljFunctorAVX.getVirial(), ljFunctorHWY.getVirial(), tolerance) << "global virial";
}

template <VectorizationPattern vecPattern>
void LJFunctorSmoothHWYTest::testLJFunctorAVXvsLJFunctorHWYOneCell(bool newton3, bool doDeleteSomeParticles, bool useUnalignedViews) {
 FMCell cellHWY;

 size_t numParticles = 23;

 ParticlePropertiesLibrary<double, size_t> PPL{_cutoff};
 PPL.addSiteType(0,1.,1.,1.);
 PPL.addSiteType(1,1.5,2.,1.);
 PPL.addSiteType(2,2.,1.,1.);
 PPL.addSiteType(3,2.5,2.,1.);
 PPL.addSiteType(4,3.,1.,1.);
 PPL.calculateMixingCoefficients();

 Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
 autopasTools::generators::RandomGenerator::fillWithParticles(cellHWY, defaultParticle, _lowCorner, _highCorner,
                                                              numParticles);

 if (true) {
   for (auto &particle : cellHWY) {
     if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
     if (particle.getID() == 11) autopas::internal::markParticleAsDeleted(particle);
     if (particle.getID() == 12) autopas::internal::markParticleAsDeleted(particle);
   }
 }

 // copy cells
 FMCell cellNoHWY(cellHWY);
 constexpr bool shifting = true;
 constexpr bool mixing = false;
 mdLib::LJFunctorSmooth<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorSmooth(_cutoff,_cutoff*0.54321/*, PPL*/);
 ljFunctorSmooth.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
 mdLib::LJFunctorSmoothHWYGS<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true, false, vecPattern> ljFunctorSmoothHWYGS(_cutoff,_cutoff*0.54321/*, PPL*/);
 ljFunctorSmoothHWYGS.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

 ASSERT_TRUE(AoSParticlesEqual(cellHWY, cellNoHWY)) << "Cells not equal after copy initialization.";

 ljFunctorSmooth.initTraversal();
 ljFunctorSmoothHWYGS.initTraversal();

 ljFunctorSmooth.SoALoader(cellNoHWY, cellNoHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
 ljFunctorSmoothHWYGS.SoALoader(cellHWY, cellHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);

 ASSERT_TRUE(SoAParticlesEqual(cellHWY._particleSoABuffer, cellNoHWY._particleSoABuffer))
     << "Cells not equal after loading.";

 if (useUnalignedViews) {
   ljFunctorSmooth.SoAFunctorSingle(cellNoHWY._particleSoABuffer.constructView(1, cellNoHWY.size()), newton3);
   ljFunctorSmoothHWYGS.SoAFunctorSingle(cellHWY._particleSoABuffer.constructView(1, cellHWY.size()), newton3);
 } else {
   ljFunctorSmooth.SoAFunctorSingle(cellNoHWY._particleSoABuffer, newton3);
   ljFunctorSmoothHWYGS.SoAFunctorSingle(cellHWY._particleSoABuffer, newton3);
 }
 ASSERT_TRUE(SoAParticlesEqual(cellHWY._particleSoABuffer, cellNoHWY._particleSoABuffer))
     << "Cells not equal after applying functor.";

 ljFunctorSmoothHWYGS.SoAExtractor(cellHWY, cellHWY._particleSoABuffer, 0);
 ljFunctorSmooth.SoAExtractor(cellNoHWY, cellNoHWY._particleSoABuffer, 0);

 ASSERT_TRUE(AoSParticlesEqual(cellHWY, cellNoHWY)) << "Cells 1 not equal after extracting.";

 ljFunctorSmoothHWYGS.endTraversal(newton3);
 ljFunctorSmooth.endTraversal(newton3);

 double tolerance = 1e-8;
 // TODO : uncomment before release
 // EXPECT_NEAR(ljFunctorAVX.getPotentialEnergy(), ljFunctorHWY.getUpot(), tolerance) << "global uPot";
 // EXPECT_NEAR(ljFunctorAVX.getVirial(), ljFunctorHWY.getVirial(), tolerance) << "global virial";
}

template <VectorizationPattern vecPattern>
void LJFunctorSmoothHWYTest::testLJFunctorAVXvsLJFunctorHWYVerlet(bool newton3, bool doDeleteSomeParticles) {

 using namespace autopas::utils::ArrayMath::literals;

 FMCell cellAVX;

 constexpr size_t numParticles = 23;

 ParticlePropertiesLibrary<double, size_t> PPL{_cutoff};
 PPL.addSiteType(0,1.,1.,1.);
 PPL.addSiteType(1,1.5,2.,1.);
 PPL.addSiteType(2,2.,1.,1.);
 PPL.addSiteType(3,2.5,2.,1.);
 PPL.addSiteType(4,3.,1.,1.);
 PPL.calculateMixingCoefficients();

 Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
 autopasTools::generators::RandomGenerator::fillWithParticles(cellAVX, defaultParticle, _lowCorner, _highCorner,
                                                              numParticles);

 if (doDeleteSomeParticles) {
   for (auto &particle : cellAVX) {
     if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
     if (particle.getID() == 11) autopas::internal::markParticleAsDeleted(particle);
     if (particle.getID() == 12) autopas::internal::markParticleAsDeleted(particle);
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
 mdLib::LJFunctorSmoothHWYGS<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true, false, vecPattern> ljFunctorSmoothHWYGS(_cutoff,_cutoff*0.54321/*, PPL*/);
 ljFunctorSmoothHWYGS.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
 mdLib::LJFunctorSmooth<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, calculateGlobals> ljFunctorSmooth(_cutoff,_cutoff*0.54321/*, PPL*/);
 ljFunctorSmooth.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

 ASSERT_TRUE(AoSParticlesEqual(cellAVX, cellHWY)) << "Cells not equal after copy initialization.";

 ljFunctorSmooth.initTraversal();
 ljFunctorSmoothHWYGS.initTraversal();

 ljFunctorSmoothHWYGS.SoALoader(cellHWY, cellHWY._particleSoABuffer, 0, /*skipSoAResize*/ false);
 ljFunctorSmooth.SoALoader(cellAVX, cellAVX._particleSoABuffer, 0, /*skipSoAResize*/ false);

 ASSERT_TRUE(SoAParticlesEqual(cellAVX._particleSoABuffer, cellHWY._particleSoABuffer))
     << "Cells not equal after loading.";

 for (size_t i = 0; i < numParticles; ++i) {
   ljFunctorSmoothHWYGS.SoAFunctorVerlet(cellHWY._particleSoABuffer, i, neighborLists[i], newton3);
   ljFunctorSmooth.SoAFunctorVerlet(cellAVX._particleSoABuffer, i, neighborLists[i], newton3);
 }

 ASSERT_TRUE(SoAParticlesEqual(cellAVX._particleSoABuffer, cellHWY._particleSoABuffer))
     << "Cells not equal after applying functor.";

 ljFunctorSmooth.SoAExtractor(cellAVX, cellAVX._particleSoABuffer, 0);
 ljFunctorSmooth.SoAExtractor(cellHWY, cellHWY._particleSoABuffer, 0);

 ASSERT_TRUE(AoSParticlesEqual(cellAVX, cellHWY)) << "Cells not equal after extracting.";

 ljFunctorSmooth.endTraversal(newton3);
 ljFunctorSmoothHWYGS.endTraversal(newton3);

 double tolerance = 1e-8;
 // TODO : uncomment before release
 // EXPECT_NEAR(ljFunctorAVX.getPotentialEnergy(), ljFunctorHWY.getUpot(), tolerance) << "global uPot";
 // EXPECT_NEAR(ljFunctorAVX.getVirial(), ljFunctorHWY.getVirial(), tolerance) << "global virial";
}

void LJFunctorSmoothHWYTest::testLJFunctorAVXvsLJFunctorHWYAoS(bool newton3, bool doDeleteSomeParticles) {
 FMCell cellHWY;

 constexpr size_t numParticles = 23;

 ParticlePropertiesLibrary<double, size_t> PPL{_cutoff};
 PPL.addSiteType(0,1.,1.,1.);
 PPL.addSiteType(1,1.5,2.,1.);
 PPL.addSiteType(2,2.,1.,1.);
 PPL.addSiteType(3,2.5,2.,1.);
 PPL.addSiteType(4,3.,1.,1.);
 PPL.calculateMixingCoefficients();

 Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
 autopasTools::generators::RandomGenerator::fillWithParticles(cellHWY, defaultParticle, _lowCorner, _highCorner,
                                                              numParticles);

 if (doDeleteSomeParticles) {
   for (auto &particle : cellHWY) {
     if (particle.getID() == 3) autopas::internal::markParticleAsDeleted(particle);
     if (particle.getID() == 11) autopas::internal::markParticleAsDeleted(particle);
     if (particle.getID() == 12) autopas::internal::markParticleAsDeleted(particle);
   }
 }

 // copy cells
 FMCell cellNoHWY(cellHWY);
 constexpr bool shifting = true;
 constexpr bool mixing = false;
 mdLib::LJFunctorSmooth<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true> ljFunctorSmooth(_cutoff,_cutoff*0.54321/*, PPL*/);
 ljFunctorSmooth.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);
 mdLib::LJFunctorSmoothHWYGS<Molecule, shifting, mixing, autopas::FunctorN3Modes::Both, true, false> ljFunctorSmoothHWYGS(_cutoff,_cutoff*0.54321/*, PPL*/);
 ljFunctorSmoothHWYGS.setParticleProperties(_epsilon * 24.0, _sigma * _sigma);

 ASSERT_TRUE(AoSParticlesEqual(cellHWY, cellNoHWY)) << "Cells not equal after copy initialization.";

 ljFunctorSmoothHWYGS.initTraversal();
 ljFunctorSmooth.initTraversal();

 for (size_t i = 0; i < numParticles; ++i) {
   for (size_t j = newton3 ? i + 1 : 0; j < numParticles; ++j) {
     if (i == j) {
       continue;
     }
     ljFunctorSmooth.AoSFunctor(cellNoHWY[i], cellNoHWY[j], newton3);
     ljFunctorSmoothHWYGS.AoSFunctor(cellHWY[i], cellHWY[j], newton3);
   }
 }

 ASSERT_TRUE(AoSParticlesEqual(cellHWY, cellNoHWY)) << "Cells not equal after applying AoSfunctor.";

 ljFunctorSmoothHWYGS.endTraversal(newton3);
 ljFunctorSmooth.endTraversal(newton3);

 double tolerance = 1e-8;
 // TODO : uncomment before release
 // EXPECT_NEAR(ljFunctorHWY.getUpot(), ljFunctorAVX.getPotentialEnergy(), tolerance) << "global uPot";
 // EXPECT_NEAR(ljFunctorHWY.getVirial(), ljFunctorAVX.getVirial(), tolerance) << "global virial";
}

TEST_P(LJFunctorSmoothHWYTest, testLJFunctorVSLJFunctorHWYAoS) {
 auto [newton3, doDeleteSomeParticle, _] = GetParam();
 testLJFunctorAVXvsLJFunctorHWYAoS(newton3, doDeleteSomeParticle);
}

TEST_P(LJFunctorSmoothHWYTest, testLJFunctorVSLJFunctorHWYVerlet) {
 // different vectorization patterns are currently not supported for Verlet Functor
 auto [newton3, doDeleteSomeParticle, vecPattern] = GetParam();
 switch (vecPattern)
 {
   case VectorizationPattern::p1xVec:
     testLJFunctorAVXvsLJFunctorHWYOneCell<VectorizationPattern::p1xVec>(newton3, doDeleteSomeParticle, false);
     break;
   case VectorizationPattern::p2xVecDiv2:
   case VectorizationPattern::pVecDiv2x2:
   case VectorizationPattern::pVecx1:
   case VectorizationPattern::pVecxVec:
     return;
   default:
     throw std::runtime_error("No vectorization pattern matched");
 }
}

TEST_P(LJFunctorSmoothHWYTest, testLJFunctorVSLJFunctorHWYOneCellAlignedAccess) {
 auto [newton3, doDeleteSomeParticle, vecPattern] = GetParam();
 switch (vecPattern)
 {
   case VectorizationPattern::p1xVec:
     testLJFunctorAVXvsLJFunctorHWYOneCell<VectorizationPattern::p1xVec>(newton3, doDeleteSomeParticle, false);
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

TEST_P(LJFunctorSmoothHWYTest, testLJFunctorVSLJFunctorHWYOneCellUseUnalignedViews) {
 auto [newton3, doDeleteSomeParticle, vecPattern] = GetParam();
 switch (vecPattern)
 {
   case VectorizationPattern::p1xVec:
     testLJFunctorAVXvsLJFunctorHWYOneCell<VectorizationPattern::p1xVec>(newton3, doDeleteSomeParticle, true);
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

TEST_P(LJFunctorSmoothHWYTest, testLJFunctorVSLJFunctorHWYTwoCellsAlignedAccess) {
 auto [newton3, doDeleteSomeParticle, vecPattern] = GetParam();
 switch (vecPattern)
 {
   case VectorizationPattern::p1xVec:
     testLJFunctorAVXvsLJFunctorHWYTwoCells<VectorizationPattern::p1xVec>(newton3, doDeleteSomeParticle, false);
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

TEST_P(LJFunctorSmoothHWYTest, testLJFunctorVSLJFunctorHWYTwoCellsUseUnalignedViews) {
 auto [newton3, doDeleteSomeParticle, vecPattern] = GetParam();
 switch (vecPattern)
 {
   case VectorizationPattern::p1xVec:
     testLJFunctorAVXvsLJFunctorHWYTwoCells<VectorizationPattern::p1xVec>(newton3, doDeleteSomeParticle, true);
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

std::vector<VectorizationPattern> patternsSmooth {
   VectorizationPattern::p1xVec,
   // VectorizationPattern::pVecDiv2x2,
   // VectorizationPattern::pVecx1,
};

std::map<VectorizationPattern, std::string> patternsToStringSmooth {
   { VectorizationPattern::p1xVec, "1xVec"},
   { VectorizationPattern::pVecx1, "Vecx1" },
};

static auto toString = [](const auto &info) {
 auto [newton3, doDeleteSomeParticle, vecPattern] = info.param;
 std::stringstream resStream;
 resStream << patternsToStringSmooth[vecPattern] << (newton3 ? "N3" : "noN3") << "_" << (doDeleteSomeParticle ? "withDeletions" : "noDeletions");
 std::string res = resStream.str();
 std::replace(res.begin(), res.end(), '-', '_');
 std::replace(res.begin(), res.end(), '.', '_');
 return res;
};

INSTANTIATE_TEST_SUITE_P(Generated, LJFunctorSmoothHWYTest, ::testing::Combine(::testing::Bool(), ::testing::Bool(), ::testing::ValuesIn(patternsSmooth)), toString);