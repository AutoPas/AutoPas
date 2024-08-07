/**
 * @file SPHTest.cpp
 * @author seckler
 * @date 22.01.18
 */

#include "SPHTest.h"

#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletLists.h"

using DensityFunctorType = sphLib::SPHCalcDensityFunctor<sphLib::SPHParticle>;

using HydroForceFunctorType = sphLib::SPHCalcHydroForceFunctor<sphLib::SPHParticle>;

TEST_F(SPHTest, testW) {
  double value = sphLib::SPHKernels::W({1., 1., 1.}, 1.);
  double shouldBeValue = 0.00944773;
  EXPECT_NEAR(value, shouldBeValue, 1e-8);

  value = sphLib::SPHKernels::W({1., .5, .25}, .5);
  shouldBeValue = 0.00151727;
  EXPECT_NEAR(value, shouldBeValue, 1e-8);
}

TEST_F(SPHTest, testGradW) {
  std::array<double, 3> value = sphLib::SPHKernels::gradW({1., 1., 1.}, 1.);
  std::array<double, 3> shouldBeValue = {-0.0213086, -0.0213086, -0.0213086};
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(value[i], shouldBeValue[i], 1e-7);
  }

  value = sphLib::SPHKernels::gradW({1., .5, .25}, .5);
  shouldBeValue = {-0.038073, -0.0190365, -0.00951825};
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(value[i], shouldBeValue[i], 1e-7);
  }
}

TEST_F(SPHTest, testSPHCalcDensityFunctor) {
  sphLib::SPHParticle sphParticle1({0., 0., 0.}, {1., .5, .25}, 1, 2.5, 0.7, 0.6);
  sphLib::SPHParticle sphParticle2({.1, .2, .3}, {-1., -.3, -.5}, 2, 1.5, 1.3, 0.8);

  DensityFunctorType densityFunctor;
  densityFunctor.AoSFunctor(sphParticle1, sphParticle2);
  // densityFunctor.AoSFunctor(sphParticle2, sphParticle1);

  EXPECT_NEAR(sphParticle1.getDensity(), 0.559026, 1e-6);
  EXPECT_NEAR(sphParticle2.getDensity(), 0.172401, 1e-6);
}

TEST_F(SPHTest, testSPHCalcDensityFunctorSoALoadExtract) {
  sphLib::SPHParticle sphParticle1({0., 0., 0.}, {1., .5, .25}, 1, 2.5, 0.7, 0.6);
  sphLib::SPHParticle sphParticle2({.1, .2, .3}, {-1., -.3, -.5}, 2, 1.5, 1.3, 0.8);

  // add particles to cell
  autopas::FullParticleCell<sphLib::SPHParticle> cell;
  cell.addParticle(sphParticle1);
  cell.addParticle(sphParticle2);

  // declare + instantiate soa
  auto &soa = cell._particleSoABuffer;

  // declare functor
  DensityFunctorType densityFunctor;

  // load soa
  densityFunctor.SoALoader(cell, soa, 0, /*skipSoAResize*/ false);

  auto massptr = soa.begin<sphLib::SPHParticle::AttributeNames::mass>();
  auto densityptr = soa.begin<sphLib::SPHParticle::AttributeNames::density>();
  auto smthlngthptr = soa.begin<sphLib::SPHParticle::AttributeNames::smth>();
  auto xptr = soa.begin<sphLib::SPHParticle::AttributeNames::posX>();
  auto yptr = soa.begin<sphLib::SPHParticle::AttributeNames::posY>();
  auto zptr = soa.begin<sphLib::SPHParticle::AttributeNames::posZ>();

  // check loading
  {
    auto iterator = cell.begin();
    for (size_t i = 0; i < 2; i++, ++iterator) {
      EXPECT_EQ(massptr[i], iterator->getMass());
      EXPECT_EQ(densityptr[i], 0.);
      EXPECT_EQ(smthlngthptr[i], iterator->getSmoothingLength());
      EXPECT_EQ(xptr[i], iterator->getR()[0]);
      EXPECT_EQ(yptr[i], iterator->getR()[1]);
      EXPECT_EQ(zptr[i], iterator->getR()[2]);
    }
  }

  // simulate functor call by changing density
  densityptr[0] = 3.;
  densityptr[1] = 2.;

  // extract soa
  densityFunctor.SoAExtractor(cell, soa, 0);

  // check extraction
  {
    std::array<double, 2> expectedDensity{3., 2.};
    auto iterator = cell.begin();
    for (size_t i = 0; i < 2; i++, ++iterator) {
      EXPECT_EQ(massptr[i], iterator->getMass());
      EXPECT_EQ(densityptr[i], iterator->getDensity());
      EXPECT_EQ(expectedDensity[i], iterator->getDensity());
      EXPECT_EQ(smthlngthptr[i], iterator->getSmoothingLength());
      EXPECT_EQ(xptr[i], iterator->getR()[0]);
      EXPECT_EQ(yptr[i], iterator->getR()[1]);
      EXPECT_EQ(zptr[i], iterator->getR()[2]);
    }
  }
}

TEST_F(SPHTest, testSPHCalcHydroForceFunctorSoALoadExtract) {
  sphLib::SPHParticle sphParticle1({0., 0., 0.}, {1., .5, .25}, 1, 2.5, 0.7, 0.6);
  sphLib::SPHParticle sphParticle2({.1, .2, .3}, {-1., -.3, -.5}, 2, 1.5, 1.3, 0.8);

  // simulate density functor call
  sphParticle1.setDensity(3.);
  sphParticle2.setDensity(2.);

  // add particles to cell
  autopas::FullParticleCell<sphLib::SPHParticle> cell;
  cell.addParticle(sphParticle1);
  cell.addParticle(sphParticle2);

  // declare + instantiate soa
  auto &soa = cell._particleSoABuffer;

  // declare functor
  HydroForceFunctorType hydroForceFunctor;

  // load soa
  hydroForceFunctor.SoALoader(cell, soa, 0, /*skipSoAResize*/ false);

  auto *const massptr = soa.begin<sphLib::SPHParticle::AttributeNames::mass>();
  auto *const densityptr = soa.begin<sphLib::SPHParticle::AttributeNames::density>();
  auto *const smthlngthptr = soa.begin<sphLib::SPHParticle::AttributeNames::smth>();
  auto *const soundSpeedptr = soa.begin<sphLib::SPHParticle::AttributeNames::soundSpeed>();
  auto *const pressureptr = soa.begin<sphLib::SPHParticle::AttributeNames::pressure>();
  auto *const vsigmaxptr = soa.begin<sphLib::SPHParticle::AttributeNames::vsigmax>();
  auto *const engDotptr = soa.begin<sphLib::SPHParticle::AttributeNames::engDot>();
  auto *const xptr = soa.begin<sphLib::SPHParticle::AttributeNames::posX>();
  auto *const yptr = soa.begin<sphLib::SPHParticle::AttributeNames::posY>();
  auto *const zptr = soa.begin<sphLib::SPHParticle::AttributeNames::posZ>();
  auto *const velXptr = soa.begin<sphLib::SPHParticle::AttributeNames::velX>();
  auto *const velYptr = soa.begin<sphLib::SPHParticle::AttributeNames::velY>();
  auto *const velZptr = soa.begin<sphLib::SPHParticle::AttributeNames::velZ>();
  auto *const accXptr = soa.begin<sphLib::SPHParticle::AttributeNames::accX>();
  auto *const accYptr = soa.begin<sphLib::SPHParticle::AttributeNames::accY>();
  auto *const accZptr = soa.begin<sphLib::SPHParticle::AttributeNames::accZ>();

  // check loading
  {
    auto iterator = cell.begin();
    for (size_t i = 0; i < 2; i++, ++iterator) {
      EXPECT_EQ(massptr[i], iterator->getMass());
      EXPECT_EQ(densityptr[i], iterator->getDensity());
      EXPECT_EQ(smthlngthptr[i], iterator->getSmoothingLength());
      EXPECT_EQ(soundSpeedptr[i], iterator->getSoundSpeed());
      EXPECT_EQ(pressureptr[i], iterator->getPressure());
      EXPECT_EQ(vsigmaxptr[i], iterator->getVSigMax());
      EXPECT_EQ(engDotptr[i], iterator->getEngDot());
      EXPECT_EQ(xptr[i], iterator->getR()[0]);
      EXPECT_EQ(yptr[i], iterator->getR()[1]);
      EXPECT_EQ(zptr[i], iterator->getR()[2]);
      EXPECT_EQ(velXptr[i], iterator->getV()[0]);
      EXPECT_EQ(velYptr[i], iterator->getV()[1]);
      EXPECT_EQ(velZptr[i], iterator->getV()[2]);
      EXPECT_EQ(accXptr[i], iterator->getAcceleration()[0]);
      EXPECT_EQ(accYptr[i], iterator->getAcceleration()[1]);
      EXPECT_EQ(accZptr[i], iterator->getAcceleration()[2]);
    }
  }

  // simulate functor call by changing density
  engDotptr[0] = -3.;
  engDotptr[1] = -2.;
  vsigmaxptr[0] = 1.;
  vsigmaxptr[1] = 1.5;
  accXptr[0] = 0.1;
  accXptr[1] = -0.1;
  accYptr[0] = 0.2;
  accYptr[1] = -0.2;
  accZptr[0] = 0.3;
  accZptr[1] = -0.3;

  // extract soa
  hydroForceFunctor.SoAExtractor(cell, soa, 0);

  // check extraction
  {
    std::array<double, 2> expectedEngdot{-3., -2.};
    std::array<double, 2> expectedvsigmax{1., 1.5};
    std::array<double, 2> expectedAccX{.1, -.1};
    std::array<double, 2> expectedAccY{.2, -.2};
    std::array<double, 2> expectedAccZ{.3, -.3};
    auto iterator = cell.begin();
    for (size_t i = 0; i < 2; i++, ++iterator) {
      EXPECT_EQ(massptr[i], iterator->getMass());
      EXPECT_EQ(densityptr[i], iterator->getDensity());
      EXPECT_EQ(smthlngthptr[i], iterator->getSmoothingLength());
      EXPECT_EQ(soundSpeedptr[i], iterator->getSoundSpeed());
      EXPECT_EQ(pressureptr[i], iterator->getPressure());
      EXPECT_EQ(expectedvsigmax[i], iterator->getVSigMax());
      EXPECT_EQ(expectedEngdot[i], iterator->getEngDot());
      EXPECT_EQ(xptr[i], iterator->getR()[0]);
      EXPECT_EQ(yptr[i], iterator->getR()[1]);
      EXPECT_EQ(zptr[i], iterator->getR()[2]);
      EXPECT_EQ(velXptr[i], iterator->getV()[0]);
      EXPECT_EQ(velYptr[i], iterator->getV()[1]);
      EXPECT_EQ(velZptr[i], iterator->getV()[2]);
      EXPECT_EQ(expectedAccX[i], iterator->getAcceleration()[0]);
      EXPECT_EQ(expectedAccY[i], iterator->getAcceleration()[1]);
      EXPECT_EQ(expectedAccZ[i], iterator->getAcceleration()[2]);
    }
  }
}

TEST_F(SPHTest, testSPHCalcDensityFunctorSoAvsAoSSingleCell) {
  autopas::FullParticleCell<sphLib::SPHParticle> cellUsingSoA;
  autopas::FullParticleCell<sphLib::SPHParticle> cellUsingAoS;
  {
    sphLib::SPHParticle defaultSphParticle({0., 0., 0.}, {1., .5, .25}, 0, 2.5, 0.7, 0.6);
    autopasTools::generators::RandomGenerator::fillWithParticles(cellUsingSoA, defaultSphParticle, {0., 0., 0.},
                                                                 {1., 1., 1.}, 30);
    autopasTools::generators::RandomGenerator::fillWithParticles(cellUsingAoS, defaultSphParticle, {0., 0., 0.},
                                                                 {1., 1., 1.}, 30);
  }

  // declare functor
  DensityFunctorType densityFunctor;

  // ------------- aos ------------------ (using newton3)
  for (auto outer = cellUsingAoS.begin(); outer != cellUsingAoS.end(); ++outer) {
    auto &p1 = *outer;

    auto inner = outer;
    ++inner;
    for (; inner != cellUsingAoS.end(); ++inner) {
      auto &p2 = *inner;

      densityFunctor.AoSFunctor(p1, p2, true);
    }
  }

  // ------------- soa ------------------ (using newton3)

  // load soa
  densityFunctor.SoALoader(cellUsingSoA, cellUsingSoA._particleSoABuffer, 0, /*skipSoAResize*/ false);

  // functors (single cell)
  densityFunctor.SoAFunctorSingle(cellUsingSoA._particleSoABuffer, true);

  // extract soa
  densityFunctor.SoAExtractor(cellUsingSoA, cellUsingSoA._particleSoABuffer, 0);

  // check same densities
  {
    auto iteratorAoS = cellUsingAoS.begin();
    auto iteratorSoA = cellUsingSoA.begin();
    for (; iteratorAoS != cellUsingAoS.end(); ++iteratorAoS, ++iteratorSoA) {
      ASSERT_TRUE(iteratorSoA != cellUsingSoA.end());
      EXPECT_NEAR(iteratorAoS->getDensity(), iteratorSoA->getDensity(), 1.e-15 * fabs(iteratorAoS->getDensity()));
    }
  }
}

TEST_F(SPHTest, testSPHCalcHydroForceFunctorSoAvsAoSSingleCell) {
  autopas::FullParticleCell<sphLib::SPHParticle> cellUsingSoA;
  autopas::FullParticleCell<sphLib::SPHParticle> cellUsingAoS;
  {
    sphLib::SPHParticle defaultSphParticle({0., 0., 0.}, {1., .5, .25}, 0, 2.5, 0.7, 0.6);
    autopasTools::generators::RandomGenerator::fillWithParticles(cellUsingSoA, defaultSphParticle, {0., 0., 0.},
                                                                 {1., 1., 1.}, 30);
    autopasTools::generators::RandomGenerator::fillWithParticles(cellUsingAoS, defaultSphParticle, {0., 0., 0.},
                                                                 {1., 1., 1.}, 30);
  }

  // simulate density functor call by setting density to sth. between 0 and 1
  {
    auto iteratorAoS = cellUsingAoS.begin();
    auto iteratorSoA = cellUsingSoA.begin();
    for (; iteratorAoS != cellUsingAoS.end(); ++iteratorAoS, ++iteratorSoA) {
      ASSERT_TRUE(iteratorSoA != cellUsingSoA.end());
      double density = static_cast<double>(rand()) / RAND_MAX;
      double pressure = static_cast<double>(rand()) / RAND_MAX;
      iteratorAoS->setDensity(density);
      iteratorAoS->setPressure(pressure);
      iteratorSoA->setDensity(density);
      iteratorSoA->setPressure(pressure);
    }
  }

  // declare functor
  HydroForceFunctorType hydroForceFunctor;

  // ------------- aos ------------------ (using newton3)
  for (auto outer = cellUsingAoS.begin(); outer != cellUsingAoS.end(); ++outer) {
    auto &p1 = *outer;

    auto inner = outer;
    ++inner;
    for (; inner != cellUsingAoS.end(); ++inner) {
      auto &p2 = *inner;

      hydroForceFunctor.AoSFunctor(p1, p2, true);
    }
  }

  // ------------- soa ------------------ (using newton3)

  // load soa
  hydroForceFunctor.SoALoader(cellUsingSoA, cellUsingSoA._particleSoABuffer, 0, /*skipSoAResize*/ false);

  // functors (single cell)
  hydroForceFunctor.SoAFunctorSingle(cellUsingSoA._particleSoABuffer, true);

  // extract soa
  hydroForceFunctor.SoAExtractor(cellUsingSoA, cellUsingSoA._particleSoABuffer, 0);

  // check same values
  {
    auto iteratorAoS = cellUsingAoS.begin();
    auto iteratorSoA = cellUsingSoA.begin();
    for (; iteratorAoS != cellUsingAoS.end(); ++iteratorAoS, ++iteratorSoA) {
      ASSERT_TRUE(iteratorSoA != cellUsingSoA.end());
      EXPECT_NEAR(iteratorAoS->getVSigMax(), iteratorSoA->getVSigMax(), 1.e-15 * fabs(iteratorAoS->getVSigMax()));
      EXPECT_NEAR(iteratorAoS->getEngDot(), iteratorSoA->getEngDot(), 1.e-15 * fabs(iteratorAoS->getEngDot()));
      EXPECT_NEAR(iteratorAoS->getAcceleration()[0], iteratorSoA->getAcceleration()[0],
                  1.e-14 * fabs(iteratorAoS->getAcceleration()[0]));
      EXPECT_NEAR(iteratorAoS->getAcceleration()[1], iteratorSoA->getAcceleration()[1],
                  1.e-14 * fabs(iteratorAoS->getAcceleration()[1]));
      EXPECT_NEAR(iteratorAoS->getAcceleration()[2], iteratorSoA->getAcceleration()[2],
                  1.e-14 * fabs(iteratorAoS->getAcceleration()[2]));
    }
  }
}

TEST_F(SPHTest, testSPHCalcDensityFunctorSoAvsAoSCellPair) {
  autopas::FullParticleCell<sphLib::SPHParticle> cellUsingSoA1;
  autopas::FullParticleCell<sphLib::SPHParticle> cellUsingSoA2;
  autopas::FullParticleCell<sphLib::SPHParticle> cellUsingAoS1;
  autopas::FullParticleCell<sphLib::SPHParticle> cellUsingAoS2;
  {
    sphLib::SPHParticle defaultSphParticle({0., 0., 0.}, {1., .5, .25}, 0, 2.5, 0.7, 0.6);
    autopasTools::generators::RandomGenerator::fillWithParticles(cellUsingSoA1, defaultSphParticle, {0., 0., 0.},
                                                                 {.5, 1., 1.}, 30);
    autopasTools::generators::RandomGenerator::fillWithParticles(cellUsingAoS1, defaultSphParticle, {0., 0., 0.},
                                                                 {.5, 1., 1.}, 30);
    autopasTools::generators::RandomGenerator::fillWithParticles(cellUsingAoS2, defaultSphParticle, {0.5, 0., 0.},
                                                                 {1., 1., 1.}, 20);
    autopasTools::generators::RandomGenerator::fillWithParticles(cellUsingSoA2, defaultSphParticle, {0.5, 0., 0.},
                                                                 {1., 1., 1.}, 20);
  }

  // declare functor
  DensityFunctorType densityFunctor;

  // ------------- aos ------------------ (using newton3)
  for (auto outer = cellUsingAoS1.begin(); outer != cellUsingAoS1.end(); ++outer) {
    auto &p1 = *outer;

    for (auto inner = cellUsingAoS2.begin(); inner != cellUsingAoS2.end(); ++inner) {
      auto &p2 = *inner;

      densityFunctor.AoSFunctor(p1, p2, true);
    }
  }

  // ------------- soa ------------------ (using newton3)

  // load soa
  densityFunctor.SoALoader(cellUsingSoA1, cellUsingSoA1._particleSoABuffer, 0, /*skipSoAResize*/ false);
  densityFunctor.SoALoader(cellUsingSoA2, cellUsingSoA2._particleSoABuffer, 0, /*skipSoAResize*/ false);

  // functors (single cell)
  densityFunctor.SoAFunctorPair(cellUsingSoA1._particleSoABuffer, cellUsingSoA2._particleSoABuffer, true);

  // extract soa
  densityFunctor.SoAExtractor(cellUsingSoA1, cellUsingSoA1._particleSoABuffer, 0);
  densityFunctor.SoAExtractor(cellUsingSoA2, cellUsingSoA2._particleSoABuffer, 0);

  // check same densities
  {
    auto iteratorAoS = cellUsingAoS1.begin();
    auto iteratorSoA = cellUsingSoA1.begin();
    for (; iteratorAoS != cellUsingAoS1.end(); ++iteratorAoS, ++iteratorSoA) {
      ASSERT_TRUE(iteratorSoA != cellUsingSoA1.end());
      EXPECT_NEAR(iteratorAoS->getDensity(), iteratorSoA->getDensity(), 1.e-15 * fabs(iteratorAoS->getDensity()));
    }
  }
  {
    auto iteratorAoS = cellUsingAoS2.begin();
    auto iteratorSoA = cellUsingSoA2.begin();
    for (; iteratorAoS != cellUsingAoS2.end(); ++iteratorAoS, ++iteratorSoA) {
      ASSERT_TRUE(iteratorSoA != cellUsingSoA2.end());
      EXPECT_NEAR(iteratorAoS->getDensity(), iteratorSoA->getDensity(), 1.e-15 * fabs(iteratorAoS->getDensity()));
    }
  }
}

TEST_F(SPHTest, testSPHCalcHydroForceFunctorSoAvsAoSCellPair) {
  autopas::FullParticleCell<sphLib::SPHParticle> cellUsingSoA1;
  autopas::FullParticleCell<sphLib::SPHParticle> cellUsingSoA2;
  autopas::FullParticleCell<sphLib::SPHParticle> cellUsingAoS1;
  autopas::FullParticleCell<sphLib::SPHParticle> cellUsingAoS2;
  {
    sphLib::SPHParticle defaultSphParticle({0., 0., 0.}, {1., .5, .25}, 0, 2.5, 0.7, 0.6);
    autopasTools::generators::RandomGenerator::fillWithParticles(cellUsingSoA1, defaultSphParticle, {0., 0., 0.},
                                                                 {.5, 1., 1.}, 30);
    autopasTools::generators::RandomGenerator::fillWithParticles(cellUsingAoS1, defaultSphParticle, {0., 0., 0.},
                                                                 {.5, 1., 1.}, 30);
    autopasTools::generators::RandomGenerator::fillWithParticles(cellUsingSoA2, defaultSphParticle, {0.5, 0., 0.},
                                                                 {1., 1., 1.}, 20);
    autopasTools::generators::RandomGenerator::fillWithParticles(cellUsingAoS2, defaultSphParticle, {0.5, 0., 0.},
                                                                 {1., 1., 1.}, 20);
  }

  // simulate density functor call by setting density to sth. between 0 and 1
  {  // cell 1
    auto iteratorAoS = cellUsingAoS1.begin();
    auto iteratorSoA = cellUsingSoA1.begin();
    for (; iteratorAoS != cellUsingAoS1.end(); ++iteratorAoS, ++iteratorSoA) {
      ASSERT_TRUE(iteratorSoA != cellUsingSoA1.end());
      double density = static_cast<double>(rand()) / RAND_MAX;
      double pressure = static_cast<double>(rand()) / RAND_MAX;
      iteratorAoS->setDensity(density);
      iteratorAoS->setPressure(pressure);
      iteratorSoA->setDensity(density);
      iteratorSoA->setPressure(pressure);
    }
  }
  {  // cell 2
    auto iteratorAoS = cellUsingAoS2.begin();
    auto iteratorSoA = cellUsingSoA2.begin();
    for (; iteratorAoS != cellUsingAoS2.end(); ++iteratorAoS, ++iteratorSoA) {
      ASSERT_TRUE(iteratorSoA != cellUsingSoA2.end());
      double density = static_cast<double>(rand()) / RAND_MAX;
      double pressure = static_cast<double>(rand()) / RAND_MAX;
      iteratorAoS->setDensity(density);
      iteratorAoS->setPressure(pressure);
      iteratorSoA->setDensity(density);
      iteratorSoA->setPressure(pressure);
    }
  }

  // declare functor
  HydroForceFunctorType hydroForceFunctor;

  // ------------- aos ------------------ (using newton3)
  for (auto outer = cellUsingAoS1.begin(); outer != cellUsingAoS1.end(); ++outer) {
    auto &p1 = *outer;

    for (auto inner = cellUsingAoS2.begin(); inner != cellUsingAoS2.end(); ++inner) {
      auto &p2 = *inner;

      hydroForceFunctor.AoSFunctor(p1, p2, true);
    }
  }

  // ------------- soa ------------------ (using newton3)

  // load soa
  hydroForceFunctor.SoALoader(cellUsingSoA1, cellUsingSoA1._particleSoABuffer, 0, /*skipSoAResize*/ false);
  hydroForceFunctor.SoALoader(cellUsingSoA2, cellUsingSoA2._particleSoABuffer, 0, /*skipSoAResize*/ false);

  // functors (single cell)
  hydroForceFunctor.SoAFunctorPair(cellUsingSoA1._particleSoABuffer, cellUsingSoA2._particleSoABuffer, true);

  // extract soa
  hydroForceFunctor.SoAExtractor(cellUsingSoA1, cellUsingSoA1._particleSoABuffer, 0);
  hydroForceFunctor.SoAExtractor(cellUsingSoA2, cellUsingSoA2._particleSoABuffer, 0);

  // check same results properties
  {
    auto iteratorAoS = cellUsingAoS1.begin();
    auto iteratorSoA = cellUsingSoA1.begin();
    for (; iteratorAoS != cellUsingAoS1.end(); ++iteratorAoS, ++iteratorSoA) {
      ASSERT_TRUE(iteratorSoA != cellUsingSoA1.end());
      EXPECT_NEAR(iteratorAoS->getVSigMax(), iteratorSoA->getVSigMax(), 1.e-14 * fabs(iteratorAoS->getVSigMax()));
      EXPECT_NEAR(iteratorAoS->getEngDot(), iteratorSoA->getEngDot(), 1.e-14 * fabs(iteratorAoS->getEngDot()));
      EXPECT_NEAR(iteratorAoS->getAcceleration()[0], iteratorSoA->getAcceleration()[0],
                  1.e-14 * fabs(iteratorAoS->getAcceleration()[0]));
      EXPECT_NEAR(iteratorAoS->getAcceleration()[1], iteratorSoA->getAcceleration()[1],
                  1.e-14 * fabs(iteratorAoS->getAcceleration()[1]));
      EXPECT_NEAR(iteratorAoS->getAcceleration()[2], iteratorSoA->getAcceleration()[2],
                  1.e-14 * fabs(iteratorAoS->getAcceleration()[2]));
    }
  }
  {
    auto iteratorAoS = cellUsingAoS2.begin();
    auto iteratorSoA = cellUsingSoA2.begin();
    for (; iteratorAoS != cellUsingAoS2.end(); ++iteratorAoS, ++iteratorSoA) {
      ASSERT_TRUE(iteratorSoA != cellUsingSoA2.end());
      EXPECT_NEAR(iteratorAoS->getVSigMax(), iteratorSoA->getVSigMax(), 1.e-14 * fabs(iteratorAoS->getVSigMax()));
      EXPECT_NEAR(iteratorAoS->getEngDot(), iteratorSoA->getEngDot(), 1.e-14 * fabs(iteratorAoS->getEngDot()));
      EXPECT_NEAR(iteratorAoS->getAcceleration()[0], iteratorSoA->getAcceleration()[0],
                  1.e-14 * fabs(iteratorAoS->getAcceleration()[0]));
      EXPECT_NEAR(iteratorAoS->getAcceleration()[1], iteratorSoA->getAcceleration()[1],
                  1.e-14 * fabs(iteratorAoS->getAcceleration()[1]));
      EXPECT_NEAR(iteratorAoS->getAcceleration()[2], iteratorSoA->getAcceleration()[2],
                  1.e-14 * fabs(iteratorAoS->getAcceleration()[2]));
    }
  }
}

TEST_F(SPHTest, testSPHCalcPressure) {
  sphLib::SPHParticle sphParticle1({0., 0., 0.}, {1., .5, .25}, 1, 2.5, 0.7, 0.6);
  sphLib::SPHParticle sphParticle2({.1, .2, .3}, {-1., -.3, -.5}, 2, 1.5, 1.3, 0.8);

  // simulate density functor call:
  sphParticle1.addDensity(0.559026);
  sphParticle2.addDensity(0.172401);

  // set pressure:
  sphParticle1.setEnergy(2.5);
  sphParticle2.setEnergy(2.5);
  sphParticle1.calcPressure();
  sphParticle2.calcPressure();

  EXPECT_NEAR(sphParticle1.getPressure(), 0.559026, 1e-6);
  EXPECT_NEAR(sphParticle2.getPressure(), 0.172401, 1e-6);
  EXPECT_NEAR(sphParticle1.getSoundSpeed(), 1.18322, 1e-5);
  EXPECT_NEAR(sphParticle2.getSoundSpeed(), 1.18322, 1e-5);
}

TEST_F(SPHTest, testSPHCalcHydroForceFunctor) {
  sphLib::SPHParticle sphParticle1({0., 0., 0.}, {1., .5, .25}, 1, 2.5, 0.7, 0.6);
  sphLib::SPHParticle sphParticle2({.1, .2, .3}, {-1., -.3, -.5}, 2, 1.5, 1.3, 0.8);

  // simulate density functor call:
  sphParticle1.addDensity(0.559026);
  sphParticle2.addDensity(0.172401);

  // set pressure:
  sphParticle1.setEnergy(2.5);
  sphParticle2.setEnergy(2.5);
  sphParticle1.setPressure(0.559026);
  sphParticle2.setPressure(0.172401);
  sphParticle1.setSoundSpeed(1.18322);
  sphParticle2.setSoundSpeed(1.18322);

  HydroForceFunctorType hydroForceFunctor;
  hydroForceFunctor.AoSFunctor(sphParticle1, sphParticle2);
  // hydroForceFunctor.AoSFunctor(sphParticle2, sphParticle1);

  EXPECT_NEAR(sphParticle1.getAcceleration().at(0), -2.26921, 1e-5);
  EXPECT_NEAR(sphParticle1.getAcceleration().at(1), -4.53843, 1e-5);
  EXPECT_NEAR(sphParticle1.getAcceleration().at(2), -6.80764, 1e-5);
  EXPECT_NEAR(sphParticle2.getAcceleration().at(0), 3.78202, 1e-5);
  EXPECT_NEAR(sphParticle2.getAcceleration().at(1), 7.56405, 1e-5);
  EXPECT_NEAR(sphParticle2.getAcceleration().at(2), 11.3461, 1e-4);
  EXPECT_NEAR(sphParticle1.getEngDot(), 5.46311, 1e-5);
  EXPECT_NEAR(sphParticle2.getEngDot(), 13.0197, 1e-4);

  sphParticle1.calcDt();
  sphParticle2.calcDt();
  EXPECT_NEAR(sphParticle1.getDt(), 0.0595165, 1e-7);
  EXPECT_NEAR(sphParticle2.getDt(), 0.110531, 1e-6);
}

TEST_F(SPHTest, testSPHCalcDensityFunctorNewton3OnOff) {
  sphLib::SPHParticle sphParticle1({0., 0., 0.}, {1., .5, .25}, 1, 2.5, 0.7, 0.6);
  sphLib::SPHParticle sphParticle2({.1, .2, .3}, {-1., -.3, -.5}, 2, 1.5, 1.3, 0.8);

  DensityFunctorType densityFunctor;
  densityFunctor.AoSFunctor(sphParticle1, sphParticle2);

  double density1 = sphParticle1.getDensity();
  double density2 = sphParticle2.getDensity();

  densityFunctor.AoSFunctor(sphParticle1, sphParticle2, false);
  densityFunctor.AoSFunctor(sphParticle2, sphParticle1, false);

  EXPECT_NEAR(sphParticle1.getDensity(), 2 * density1, 1e-6);
  EXPECT_NEAR(sphParticle2.getDensity(), 2 * density2, 1e-6);
}

TEST_F(SPHTest, testSPHCalcHydroForceFunctorNewton3OnOff) {
  sphLib::SPHParticle sphParticle1({0., 0., 0.}, {1., .5, .25}, 1, 2.5, 0.7, 0.6);
  sphLib::SPHParticle sphParticle2({.1, .2, .3}, {-1., -.3, -.5}, 2, 1.5, 1.3, 0.8);

  // simulate density functor call:
  sphParticle1.addDensity(0.559026);
  sphParticle2.addDensity(0.172401);

  // set pressure:
  sphParticle1.setEnergy(2.5);
  sphParticle2.setEnergy(2.5);
  sphParticle1.setPressure(0.559026);
  sphParticle2.setPressure(0.172401);
  sphParticle1.setSoundSpeed(1.18322);
  sphParticle2.setSoundSpeed(1.18322);

  // make copies of the particles
  sphLib::SPHParticle sphParticle3 = sphParticle1;
  sphLib::SPHParticle sphParticle4 = sphParticle2;

  HydroForceFunctorType hydroForceFunctor;
  // using newton 3
  hydroForceFunctor.AoSFunctor(sphParticle1, sphParticle2);
  // without newton 3

  // particle 3:
  hydroForceFunctor.AoSFunctor(sphParticle3, sphParticle4, false);

  EXPECT_NEAR(sphParticle3.getAcceleration().at(0), sphParticle1.getAcceleration().at(0), 1e-10);
  EXPECT_NEAR(sphParticle3.getAcceleration().at(1), sphParticle1.getAcceleration().at(1), 1e-10);
  EXPECT_NEAR(sphParticle3.getAcceleration().at(2), sphParticle1.getAcceleration().at(2), 1e-10);
  EXPECT_NEAR(sphParticle3.getEngDot(), sphParticle1.getEngDot(), 1e-10);
  sphParticle1.calcDt();
  sphParticle3.calcDt();

  EXPECT_NEAR(sphParticle3.getDt(), sphParticle1.getDt(), 1e-10);

  EXPECT_NEAR(sphParticle4.getAcceleration().at(0), 0., 1e-10);
  EXPECT_NEAR(sphParticle4.getAcceleration().at(1), 0., 1e-10);
  EXPECT_NEAR(sphParticle4.getAcceleration().at(2), 0., 1e-10);
  EXPECT_NEAR(sphParticle4.getEngDot(), 0., 1e-10);

  // particle 4:
  hydroForceFunctor.AoSFunctor(sphParticle4, sphParticle3, false);

  EXPECT_NEAR(sphParticle4.getAcceleration().at(0), sphParticle2.getAcceleration().at(0), 1e-10);
  EXPECT_NEAR(sphParticle4.getAcceleration().at(1), sphParticle2.getAcceleration().at(1), 1e-10);
  EXPECT_NEAR(sphParticle4.getAcceleration().at(2), sphParticle2.getAcceleration().at(2), 1e-10);
  EXPECT_NEAR(sphParticle4.getEngDot(), sphParticle2.getEngDot(), 1e-10);

  sphParticle2.calcDt();
  sphParticle4.calcDt();
  EXPECT_NEAR(sphParticle4.getDt(), sphParticle2.getDt(), 1e-10);
}

void densityCheck(autopas::VerletLists<sphLib::SPHParticle> &verletLists,
                  autopas::LinkedCells<sphLib::SPHParticle> &linkedCells, size_t numMolecules, double relErrTolerance) {
  std::vector<double> densityVerlet(numMolecules), densityLinked(numMolecules);
  /* get and sort by id, the */
  for (auto it = verletLists.begin(); it.isValid(); ++it) {
    sphLib::SPHParticle &m = *it;
    densityVerlet.at(m.getID()) = m.getDensity();
  }

  for (auto it = linkedCells.begin(); it.isValid(); ++it) {
    sphLib::SPHParticle &m = *it;
    densityLinked.at(m.getID()) = m.getDensity();
  }

  for (unsigned long i = 0; i < numMolecules; ++i) {
    double d1 = densityVerlet[i];
    double d2 = densityLinked[i];
    EXPECT_NEAR(d1, d2, std::fabs(d1 * relErrTolerance));
  }
}

void hydroInit(autopas::VerletLists<sphLib::SPHParticle> &verletLists) {
  for (auto itVerlet = verletLists.begin(); itVerlet.isValid(); ++itVerlet) {
    double density = static_cast<double>(rand()) / RAND_MAX;
    double pressure = static_cast<double>(rand()) / RAND_MAX;
    itVerlet->setDensity(density);
    itVerlet->setPressure(pressure);
  }
}

void hydroCheck(autopas::VerletLists<sphLib::SPHParticle> &verletLists,
                autopas::LinkedCells<sphLib::SPHParticle> &linkedCells, size_t numMolecules, double relErrTolerance) {
  std::vector<double> vsigmaxVerlet(numMolecules), vsigmaxLinked(numMolecules);
  std::vector<double> engdotVerlet(numMolecules), engdotLinked(numMolecules);
  std::vector<std::array<double, 3>> accVerlet(numMolecules), accLinked(numMolecules);
  /* get and sort by id, the */
  for (auto it = verletLists.begin(); it.isValid(); ++it) {
    sphLib::SPHParticle &m = *it;
    vsigmaxVerlet.at(m.getID()) = m.getVSigMax();
    engdotVerlet.at(m.getID()) = m.getEngDot();
    accVerlet.at(m.getID()) = m.getAcceleration();
  }

  for (auto it = verletLists.begin(); it.isValid(); ++it) {
    sphLib::SPHParticle &m = *it;
    vsigmaxLinked.at(m.getID()) = m.getVSigMax();
    engdotLinked.at(m.getID()) = m.getEngDot();
    accLinked.at(m.getID()) = m.getAcceleration();
  }

  for (unsigned long i = 0; i < numMolecules; ++i) {
    EXPECT_NEAR(vsigmaxVerlet[i], vsigmaxLinked[i], relErrTolerance * fabs(vsigmaxLinked[i]));
    EXPECT_NEAR(engdotVerlet[i], engdotLinked[i], relErrTolerance * fabs(engdotLinked[i]));
    EXPECT_NEAR(accVerlet[i][0], accLinked[i][0], relErrTolerance * fabs(accLinked[i][0]));
    EXPECT_NEAR(accVerlet[i][1], accLinked[i][1], relErrTolerance * fabs(accLinked[i][1]));
    EXPECT_NEAR(accVerlet[i][2], accLinked[i][2], relErrTolerance * fabs(accLinked[i][2]));
  }
}

template <typename FunctorType, typename InitType, typename CheckType>
void testVerLetVsLC(FunctorType &fnctr, InitType init, CheckType check, autopas::DataLayoutOption dataLayoutOption) {
  unsigned int numMolecules = 50;
  double relErrTolerance = 1e-10;
  double cutoff = 1.;
  using sphLib::SPHParticle;

  autopas::VerletLists<SPHParticle> verletLists({0., 0., 0.}, {5., 5., 5.}, cutoff, 0.5, 1);
  autopas::LinkedCells<sphLib::SPHParticle> linkedCells({0., 0., 0.}, {5., 5., 5.}, cutoff, 0.5, 1., 1);

  sphLib::SPHParticle defaultSPHParticle({0., 0., 0.}, {1., .5, .25}, 0, 2.5,
                                         cutoff / sphLib::SPHKernels::getKernelSupportRadius(), 0.6);
  autopasTools::generators::RandomGenerator::fillWithParticles(verletLists, defaultSPHParticle, verletLists.getBoxMin(),
                                                               verletLists.getBoxMax(), numMolecules);

  // init particles in verlet list container
  init(verletLists);

  // now fill second container with the molecules from the first one, because otherwise we generate new particles
  for (auto it = verletLists.begin(); it.isValid(); ++it) {
    linkedCells.addParticle(*it);
  }

  if (dataLayoutOption == autopas::DataLayoutOption::aos) {
    autopas::LCC08Traversal<autopas::FullParticleCell<SPHParticle>, FunctorType> traversalLJ(
        linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &fnctr, linkedCells.getInteractionLength(),
        linkedCells.getCellBlock().getCellLength(), autopas::DataLayoutOption::aos, true);

    autopas::VLListIterationTraversal<autopas::FullParticleCell<SPHParticle>, FunctorType> traversalLJVerlet(
        &fnctr, autopas::DataLayoutOption::aos, true);

    verletLists.rebuildNeighborLists(&traversalLJVerlet);
    verletLists.computeInteractions(&traversalLJVerlet);
    linkedCells.computeInteractions(&traversalLJ);
  } else {
    autopas::LCC08Traversal<autopas::FullParticleCell<SPHParticle>, FunctorType> traversalLJ(
        linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &fnctr, linkedCells.getInteractionLength(),
        linkedCells.getCellBlock().getCellLength(), autopas::DataLayoutOption::soa, true);
    autopas::VLListIterationTraversal<autopas::FullParticleCell<SPHParticle>, FunctorType> traversalLJVerlet(
        &fnctr, autopas::DataLayoutOption::soa, true);

    verletLists.rebuildNeighborLists(&traversalLJVerlet);
    verletLists.computeInteractions(&traversalLJVerlet);
    linkedCells.computeInteractions(&traversalLJ);
  }
  check(verletLists, linkedCells, numMolecules, relErrTolerance);
}

TEST_P(SPHTest, testVerletVsLC) {
  auto params = GetParam();
  auto [dataLayoutOption, sphFunctorType] = params;
  if (sphFunctorType == SPHFunctorType::density) {
    DensityFunctorType densityFunctor;
    testVerLetVsLC(
        densityFunctor, [](auto & /*ignored*/) {}, densityCheck, dataLayoutOption);
  } else {
    HydroForceFunctorType hydroForceFunctor;
    testVerLetVsLC(hydroForceFunctor, hydroInit, hydroCheck, dataLayoutOption);
  }
}

INSTANTIATE_TEST_SUITE_P(Generated, SPHTest,
                         ::testing::Combine(::testing::Values(autopas::DataLayoutOption::aos,
                                                              autopas::DataLayoutOption::soa),
                                            ::testing::Values(SPHFunctorType::density, SPHFunctorType::hydro)),
                         SPHTest::PrintToStringParamName());