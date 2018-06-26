//
// Created by seckler on 22.01.18.
//

#include "SPHTest.h"
#include "autopasIncludes.h"
#include "sph/autopassph.h"
#include "testingHelpers/RandomGenerator.h"

TEST_F(SPHTest, testW) {
  double value = autopas::sph::SPHKernels::W({1., 1., 1.}, 1.);
  double should_be_value = 0.00944773;
  EXPECT_NEAR(value, should_be_value, 1e-8);

  value = autopas::sph::SPHKernels::W({1., .5, .25}, .5);
  should_be_value = 0.00151727;
  EXPECT_NEAR(value, should_be_value, 1e-8);
}

TEST_F(SPHTest, testGradW) {
  std::array<double, 3> value = autopas::sph::SPHKernels::gradW({1., 1., 1.}, 1.);
  std::array<double, 3> should_be_value = {-0.0213086, -0.0213086, -0.0213086};
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(value[i], should_be_value[i], 1e-7);
  }

  value = autopas::sph::SPHKernels::gradW({1., .5, .25}, .5);
  should_be_value = {-0.038073, -0.0190365, -0.00951825};
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(value[i], should_be_value[i], 1e-7);
  }
}

TEST_F(SPHTest, testSPHCalcDensityFunctor) {
  autopas::sph::SPHParticle sphParticle1({0., 0., 0.}, {1., .5, .25}, 1, 2.5, 0.7, 0.6);
  autopas::sph::SPHParticle sphParticle2({.1, .2, .3}, {-1., -.3, -.5}, 2, 1.5, 1.3, 0.8);

  autopas::sph::SPHCalcDensityFunctor densityFunctor;
  densityFunctor.AoSFunctor(sphParticle1, sphParticle2);
  // densityFunctor.AoSFunctor(sphParticle2, sphParticle1);

  EXPECT_NEAR(sphParticle1.getDensity(), 0.559026, 1e-6);
  EXPECT_NEAR(sphParticle2.getDensity(), 0.172401, 1e-6);
}

TEST_F(SPHTest, testSPHCalcDensityFunctorSoALoadExtract) {
  autopas::sph::SPHParticle sphParticle1({0., 0., 0.}, {1., .5, .25}, 1, 2.5, 0.7, 0.6);
  autopas::sph::SPHParticle sphParticle2({.1, .2, .3}, {-1., -.3, -.5}, 2, 1.5, 1.3, 0.8);

  // add particles to cell
  autopas::FullParticleCell<autopas::sph::SPHParticle> cell;
  cell.addParticle(sphParticle1);
  cell.addParticle(sphParticle2);

  // declare + instantiate soa
  auto &soa = cell._particleSoABuffer;

  // declare functor
  autopas::sph::SPHCalcDensityFunctor densityFunctor;

  // load soa
  densityFunctor.SoALoader(cell, soa, 0);

  auto massptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::mass>();
  auto densityptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::density>();
  auto smthlngthptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::smth>();
  auto xptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::posX>();
  auto yptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::posY>();
  auto zptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::posZ>();

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
  autopas::sph::SPHParticle sphParticle1({0., 0., 0.}, {1., .5, .25}, 1, 2.5, 0.7, 0.6);
  autopas::sph::SPHParticle sphParticle2({.1, .2, .3}, {-1., -.3, -.5}, 2, 1.5, 1.3, 0.8);

  // simulate density functor call
  sphParticle1.setDensity(3.);
  sphParticle2.setDensity(2.);

  // add particles to cell
  autopas::FullParticleCell<autopas::sph::SPHParticle> cell;
  cell.addParticle(sphParticle1);
  cell.addParticle(sphParticle2);

  // declare + instantiate soa
  auto &soa = cell._particleSoABuffer;

  // declare functor
  autopas::sph::SPHCalcHydroForceFunctor hydroForceFunctor;

  // load soa
  hydroForceFunctor.SoALoader(cell, soa, 0);

  auto massptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::mass>();
  auto densityptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::density>();
  auto smthlngthptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::smth>();
  auto soundSpeedptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::soundSpeed>();
  auto pressureptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::pressure>();
  auto vsigmaxptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::vsigmax>();
  auto engDotptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::engDot>();
  auto xptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::posX>();
  auto yptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::posY>();
  auto zptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::posZ>();
  auto velXptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::velX>();
  auto velYptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::velY>();
  auto velZptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::velZ>();
  auto accXptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::accX>();
  auto accYptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::accY>();
  auto accZptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::accZ>();

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
  autopas::FullParticleCell<autopas::sph::SPHParticle> cell_using_soa;
  autopas::FullParticleCell<autopas::sph::SPHParticle> cell_using_aos;
  {
    autopas::sph::SPHParticle defaultSphParticle({0., 0., 0.}, {1., .5, .25}, 1, 2.5, 0.7, 0.6);
    RandomGenerator::fillWithParticles(cell_using_soa, defaultSphParticle, {0., 0., 0.}, {1., 1., 1.}, 30);
    RandomGenerator::fillWithParticles(cell_using_aos, defaultSphParticle, {0., 0., 0.}, {1., 1., 1.}, 30);
  }

  // declare functor
  autopas::sph::SPHCalcDensityFunctor densityFunctor;

  // ------------- aos ------------------ (using newton3)
  for (auto outer = cell_using_aos.begin(); outer.isValid(); ++outer) {
    auto &p1 = *outer;

    auto inner = outer;
    ++inner;
    for (; inner.isValid(); ++inner) {
      auto &p2 = *inner;

      densityFunctor.AoSFunctor(p1, p2, true);
    }
  }

  // ------------- soa ------------------ (using newton3)

  // load soa
  densityFunctor.SoALoader(cell_using_soa, cell_using_soa._particleSoABuffer, 0);

  // functors (single cell)
  densityFunctor.SoAFunctor(cell_using_soa._particleSoABuffer, true);

  // extract soa
  densityFunctor.SoAExtractor(cell_using_soa, cell_using_soa._particleSoABuffer);

  // check same densities
  {
    auto iteratoraos = cell_using_aos.begin();
    auto iteratorsoa = cell_using_soa.begin();
    for (; iteratoraos.isValid(); ++iteratoraos, ++iteratorsoa) {
      ASSERT_TRUE(iteratorsoa.isValid());
      EXPECT_NEAR(iteratoraos->getDensity(), iteratorsoa->getDensity(), 1.e-15 * fabs(iteratoraos->getDensity()));
    }
  }
}

TEST_F(SPHTest, testSPHCalcHydroForceFunctorSoAvsAoSSingleCell) {
  autopas::FullParticleCell<autopas::sph::SPHParticle> cell_using_soa;
  autopas::FullParticleCell<autopas::sph::SPHParticle> cell_using_aos;
  {
    autopas::sph::SPHParticle defaultSphParticle({0., 0., 0.}, {1., .5, .25}, 1, 2.5, 0.7, 0.6);
    RandomGenerator::fillWithParticles(cell_using_soa, defaultSphParticle, {0., 0., 0.}, {1., 1., 1.}, 30);
    RandomGenerator::fillWithParticles(cell_using_aos, defaultSphParticle, {0., 0., 0.}, {1., 1., 1.}, 30);
  }

  // simulate density functor call by setting density to sth. between 0 and 1
  {
    auto iteratoraos = cell_using_aos.begin();
    auto iteratorsoa = cell_using_soa.begin();
    for (; iteratoraos.isValid(); ++iteratoraos, ++iteratorsoa) {
      ASSERT_TRUE(iteratorsoa.isValid());
      double density = static_cast<double>(rand()) / RAND_MAX;
      double pressure = static_cast<double>(rand()) / RAND_MAX;
      iteratoraos->setDensity(density);
      iteratoraos->setPressure(pressure);
      iteratorsoa->setDensity(density);
      iteratorsoa->setPressure(pressure);
    }
  }

  // declare functor
  autopas::sph::SPHCalcHydroForceFunctor hydroForceFunctor;

  // ------------- aos ------------------ (using newton3)
  for (auto outer = cell_using_aos.begin(); outer.isValid(); ++outer) {
    auto &p1 = *outer;

    auto inner = outer;
    ++inner;
    for (; inner.isValid(); ++inner) {
      auto &p2 = *inner;

      hydroForceFunctor.AoSFunctor(p1, p2, true);
    }
  }

  // ------------- soa ------------------ (using newton3)

  // load soa
  hydroForceFunctor.SoALoader(cell_using_soa, cell_using_soa._particleSoABuffer, 0);

  // functors (single cell)
  hydroForceFunctor.SoAFunctor(cell_using_soa._particleSoABuffer, true);

  // extract soa
  hydroForceFunctor.SoAExtractor(cell_using_soa, cell_using_soa._particleSoABuffer);

  // check same densities
  {
    auto iteratoraos = cell_using_aos.begin();
    auto iteratorsoa = cell_using_soa.begin();
    for (; iteratoraos.isValid(); ++iteratoraos, ++iteratorsoa) {
      ASSERT_TRUE(iteratorsoa.isValid());
      EXPECT_NEAR(iteratoraos->getVSigMax(), iteratorsoa->getVSigMax(), 1.e-15 * fabs(iteratoraos->getVSigMax()));
      EXPECT_NEAR(iteratoraos->getEngDot(), iteratorsoa->getEngDot(), 1.e-15 * fabs(iteratoraos->getEngDot()));
      EXPECT_NEAR(iteratoraos->getAcceleration()[0], iteratorsoa->getAcceleration()[0],
                  1.e-14 * fabs(iteratoraos->getAcceleration()[0]));
      EXPECT_NEAR(iteratoraos->getAcceleration()[1], iteratorsoa->getAcceleration()[1],
                  1.e-14 * fabs(iteratoraos->getAcceleration()[1]));
      EXPECT_NEAR(iteratoraos->getAcceleration()[2], iteratorsoa->getAcceleration()[2],
                  1.e-14 * fabs(iteratoraos->getAcceleration()[2]));
    }
  }
}

TEST_F(SPHTest, testSPHCalcDensityFunctorSoAvsAoSCellPair) {
  autopas::FullParticleCell<autopas::sph::SPHParticle> cell_using_soa1;
  autopas::FullParticleCell<autopas::sph::SPHParticle> cell_using_soa2;
  autopas::FullParticleCell<autopas::sph::SPHParticle> cell_using_aos1;
  autopas::FullParticleCell<autopas::sph::SPHParticle> cell_using_aos2;
  {
    autopas::sph::SPHParticle defaultSphParticle({0., 0., 0.}, {1., .5, .25}, 1, 2.5, 0.7, 0.6);
    RandomGenerator::fillWithParticles(cell_using_soa1, defaultSphParticle, {0., 0., 0.}, {.5, 1., 1.}, 30);
    RandomGenerator::fillWithParticles(cell_using_aos1, defaultSphParticle, {0., 0., 0.}, {.5, 1., 1.}, 30);
    RandomGenerator::fillWithParticles(cell_using_soa2, defaultSphParticle, {0.5, 0., 0.}, {1., 1., 1.}, 20);
    RandomGenerator::fillWithParticles(cell_using_aos2, defaultSphParticle, {0.5, 0., 0.}, {1., 1., 1.}, 20);
  }

  // declare functor
  autopas::sph::SPHCalcDensityFunctor densityFunctor;

  // ------------- aos ------------------ (using newton3)
  for (auto outer = cell_using_aos1.begin(); outer.isValid(); ++outer) {
    auto &p1 = *outer;

    for (auto inner = cell_using_aos2.begin(); inner.isValid(); ++inner) {
      auto &p2 = *inner;

      densityFunctor.AoSFunctor(p1, p2, true);
    }
  }

  // ------------- soa ------------------ (using newton3)

  // load soa
  densityFunctor.SoALoader(cell_using_soa1, cell_using_soa1._particleSoABuffer, 0);
  densityFunctor.SoALoader(cell_using_soa2, cell_using_soa2._particleSoABuffer, 0);

  // functors (single cell)
  densityFunctor.SoAFunctor(cell_using_soa1._particleSoABuffer, cell_using_soa2._particleSoABuffer, true);

  // extract soa
  densityFunctor.SoAExtractor(cell_using_soa1, cell_using_soa1._particleSoABuffer);
  densityFunctor.SoAExtractor(cell_using_soa2, cell_using_soa2._particleSoABuffer);

  // check same densities
  {
    auto iteratoraos = cell_using_aos1.begin();
    auto iteratorsoa = cell_using_soa1.begin();
    for (; iteratoraos.isValid(); ++iteratoraos, ++iteratorsoa) {
      ASSERT_TRUE(iteratorsoa.isValid());
      EXPECT_NEAR(iteratoraos->getDensity(), iteratorsoa->getDensity(), 1.e-15 * fabs(iteratoraos->getDensity()));
    }
  }
  {
    auto iteratoraos = cell_using_aos2.begin();
    auto iteratorsoa = cell_using_soa2.begin();
    for (; iteratoraos.isValid(); ++iteratoraos, ++iteratorsoa) {
      ASSERT_TRUE(iteratorsoa.isValid());
      EXPECT_NEAR(iteratoraos->getDensity(), iteratorsoa->getDensity(), 1.e-15 * fabs(iteratoraos->getDensity()));
    }
  }
}

TEST_F(SPHTest, testSPHCalcHydroForceFunctorSoAvsAoSCellPair) {
  autopas::FullParticleCell<autopas::sph::SPHParticle> cell_using_soa1;
  autopas::FullParticleCell<autopas::sph::SPHParticle> cell_using_soa2;
  autopas::FullParticleCell<autopas::sph::SPHParticle> cell_using_aos1;
  autopas::FullParticleCell<autopas::sph::SPHParticle> cell_using_aos2;
  {
    autopas::sph::SPHParticle defaultSphParticle({0., 0., 0.}, {1., .5, .25}, 1, 2.5, 0.7, 0.6);
    RandomGenerator::fillWithParticles(cell_using_soa1, defaultSphParticle, {0., 0., 0.}, {.5, 1., 1.}, 30);
    RandomGenerator::fillWithParticles(cell_using_aos1, defaultSphParticle, {0., 0., 0.}, {.5, 1., 1.}, 30);
    RandomGenerator::fillWithParticles(cell_using_soa2, defaultSphParticle, {0.5, 0., 0.}, {1., 1., 1.}, 20);
    RandomGenerator::fillWithParticles(cell_using_aos2, defaultSphParticle, {0.5, 0., 0.}, {1., 1., 1.}, 20);
  }

  // simulate density functor call by setting density to sth. between 0 and 1
  {  // cell 1
    auto iteratoraos = cell_using_aos1.begin();
    auto iteratorsoa = cell_using_soa1.begin();
    for (; iteratoraos.isValid(); ++iteratoraos, ++iteratorsoa) {
      ASSERT_TRUE(iteratorsoa.isValid());
      double density = static_cast<double>(rand()) / RAND_MAX;
      double pressure = static_cast<double>(rand()) / RAND_MAX;
      iteratoraos->setDensity(density);
      iteratoraos->setPressure(pressure);
      iteratorsoa->setDensity(density);
      iteratorsoa->setPressure(pressure);
    }
  }
  {  // cell 2
    auto iteratoraos = cell_using_aos2.begin();
    auto iteratorsoa = cell_using_soa2.begin();
    for (; iteratoraos.isValid(); ++iteratoraos, ++iteratorsoa) {
      ASSERT_TRUE(iteratorsoa.isValid());
      double density = static_cast<double>(rand()) / RAND_MAX;
      double pressure = static_cast<double>(rand()) / RAND_MAX;
      iteratoraos->setDensity(density);
      iteratoraos->setPressure(pressure);
      iteratorsoa->setDensity(density);
      iteratorsoa->setPressure(pressure);
    }
  }

  // declare functor
  autopas::sph::SPHCalcHydroForceFunctor hydroForceFunctor;

  // ------------- aos ------------------ (using newton3)
  for (auto outer = cell_using_aos1.begin(); outer.isValid(); ++outer) {
    auto &p1 = *outer;

    for (auto inner = cell_using_aos2.begin(); inner.isValid(); ++inner) {
      auto &p2 = *inner;

      hydroForceFunctor.AoSFunctor(p1, p2, true);
    }
  }

  // ------------- soa ------------------ (using newton3)

  // load soa
  hydroForceFunctor.SoALoader(cell_using_soa1, cell_using_soa1._particleSoABuffer, 0);
  hydroForceFunctor.SoALoader(cell_using_soa2, cell_using_soa2._particleSoABuffer, 0);

  // functors (single cell)
  hydroForceFunctor.SoAFunctor(cell_using_soa1._particleSoABuffer, cell_using_soa2._particleSoABuffer, true);

  // extract soa
  hydroForceFunctor.SoAExtractor(cell_using_soa1, cell_using_soa1._particleSoABuffer);
  hydroForceFunctor.SoAExtractor(cell_using_soa2, cell_using_soa2._particleSoABuffer);

  // check same results properties
  {
    auto iteratoraos = cell_using_aos1.begin();
    auto iteratorsoa = cell_using_soa1.begin();
    for (; iteratoraos.isValid(); ++iteratoraos, ++iteratorsoa) {
      ASSERT_TRUE(iteratorsoa.isValid());
      EXPECT_NEAR(iteratoraos->getVSigMax(), iteratorsoa->getVSigMax(), 1.e-14 * fabs(iteratoraos->getVSigMax()));
      EXPECT_NEAR(iteratoraos->getEngDot(), iteratorsoa->getEngDot(), 1.e-14 * fabs(iteratoraos->getEngDot()));
      EXPECT_NEAR(iteratoraos->getAcceleration()[0], iteratorsoa->getAcceleration()[0],
                  1.e-14 * fabs(iteratoraos->getAcceleration()[0]));
      EXPECT_NEAR(iteratoraos->getAcceleration()[1], iteratorsoa->getAcceleration()[1],
                  1.e-14 * fabs(iteratoraos->getAcceleration()[1]));
      EXPECT_NEAR(iteratoraos->getAcceleration()[2], iteratorsoa->getAcceleration()[2],
                  1.e-14 * fabs(iteratoraos->getAcceleration()[2]));
    }
  }
  {
    auto iteratoraos = cell_using_aos2.begin();
    auto iteratorsoa = cell_using_soa2.begin();
    for (; iteratoraos.isValid(); ++iteratoraos, ++iteratorsoa) {
      ASSERT_TRUE(iteratorsoa.isValid());
      EXPECT_NEAR(iteratoraos->getVSigMax(), iteratorsoa->getVSigMax(), 1.e-14 * fabs(iteratoraos->getVSigMax()));
      EXPECT_NEAR(iteratoraos->getEngDot(), iteratorsoa->getEngDot(), 1.e-14 * fabs(iteratoraos->getEngDot()));
      EXPECT_NEAR(iteratoraos->getAcceleration()[0], iteratorsoa->getAcceleration()[0],
                  1.e-14 * fabs(iteratoraos->getAcceleration()[0]));
      EXPECT_NEAR(iteratoraos->getAcceleration()[1], iteratorsoa->getAcceleration()[1],
                  1.e-14 * fabs(iteratoraos->getAcceleration()[1]));
      EXPECT_NEAR(iteratoraos->getAcceleration()[2], iteratorsoa->getAcceleration()[2],
                  1.e-14 * fabs(iteratoraos->getAcceleration()[2]));
    }
  }
}

TEST_F(SPHTest, testSPHCalcPressure) {
  autopas::sph::SPHParticle sphParticle1({0., 0., 0.}, {1., .5, .25}, 1, 2.5, 0.7, 0.6);
  autopas::sph::SPHParticle sphParticle2({.1, .2, .3}, {-1., -.3, -.5}, 2, 1.5, 1.3, 0.8);

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
  autopas::sph::SPHParticle sphParticle1({0., 0., 0.}, {1., .5, .25}, 1, 2.5, 0.7, 0.6);
  autopas::sph::SPHParticle sphParticle2({.1, .2, .3}, {-1., -.3, -.5}, 2, 1.5, 1.3, 0.8);

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

  autopas::sph::SPHCalcHydroForceFunctor hydroForceFunctor;
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
  autopas::sph::SPHParticle sphParticle1({0., 0., 0.}, {1., .5, .25}, 1, 2.5, 0.7, 0.6);
  autopas::sph::SPHParticle sphParticle2({.1, .2, .3}, {-1., -.3, -.5}, 2, 1.5, 1.3, 0.8);

  autopas::sph::SPHCalcDensityFunctor densityFunctor;
  densityFunctor.AoSFunctor(sphParticle1, sphParticle2);

  double density1 = sphParticle1.getDensity();
  double density2 = sphParticle2.getDensity();

  densityFunctor.AoSFunctor(sphParticle1, sphParticle2, false);
  densityFunctor.AoSFunctor(sphParticle2, sphParticle1, false);

  EXPECT_NEAR(sphParticle1.getDensity(), 2 * density1, 1e-6);
  EXPECT_NEAR(sphParticle2.getDensity(), 2 * density2, 1e-6);
}

TEST_F(SPHTest, testSPHCalcHydroForceFunctorNewton3OnOff) {
  autopas::sph::SPHParticle sphParticle1({0., 0., 0.}, {1., .5, .25}, 1, 2.5, 0.7, 0.6);
  autopas::sph::SPHParticle sphParticle2({.1, .2, .3}, {-1., -.3, -.5}, 2, 1.5, 1.3, 0.8);

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
  autopas::sph::SPHParticle sphParticle3 = sphParticle1;
  autopas::sph::SPHParticle sphParticle4 = sphParticle2;

  autopas::sph::SPHCalcHydroForceFunctor hydroForceFunctor;
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

TEST_F(SPHTest, testVerletVsLCAoS) {
  unsigned int numMolecules = 500;
  double rel_err_tolerance = 1e-10;
  double cutoff = 1.;
  using autopas::sph::SPHParticle;

  autopas::VerletLists<SPHParticle> _verletLists({0., 0., 0.}, {5., 5., 5.}, cutoff, 0.5, 3);
  autopas::LinkedCells<SPHParticle, autopas::FullParticleCell<SPHParticle>> _linkedCells({0., 0., 0.}, {5., 5., 5.},
                                                                                         cutoff);

  autopas::sph::SPHParticle defaultSPHParticle({0., 0., 0.}, {1., .5, .25}, 1, 2.5, 0.7, 0.6);
  RandomGenerator::fillWithParticles(_verletLists, defaultSPHParticle, numMolecules);
  // now fill second container with the molecules from the first one, because
  // otherwise we generate new particles
  for (auto it = _verletLists.begin(); it.isValid(); ++it) {
    _linkedCells.addParticle(*it);
  }

  autopas::sph::SPHCalcDensityFunctor densityFunctor;

  _verletLists.iteratePairwiseAoS(&densityFunctor);
  _linkedCells.iteratePairwiseAoS(&densityFunctor);


  auto itDirect = _verletLists.begin();
  auto itLinked = _linkedCells.begin();

  std::vector<std::array<double, 3>> forcesDirect(numMolecules), forcesLinked(numMolecules);
  // get and sort by id, the
  for (auto it = _verletLists.begin(); it.isValid(); ++it) {
    SPHParticle &m = *it;
    forcesDirect.at(m.getID()) = m.getF();
  }

  for (auto it = _linkedCells.begin(); it.isValid(); ++it) {
    SPHParticle &m = *it;
    forcesLinked.at(m.getID()) = m.getF();
  }

  for (unsigned long i = 0; i < numMolecules; ++i) {
    for (int d = 0; d < 3; ++d) {
      double f1 = forcesDirect[i][d];
      double f2 = forcesLinked[i][d];
      EXPECT_NEAR(f1, f2, std::fabs(f1 * rel_err_tolerance));
    }
  }
}

TEST_F(SPHTest, testVerletVsLCSoA) {
  unsigned int numMolecules = 500;
  double rel_err_tolerance = 1e-10;
  double cutoff = 1.;
  using autopas::sph::SPHParticle;

  autopas::VerletLists<SPHParticle> _verletLists({0., 0., 0.}, {5., 5., 5.}, cutoff, 0.5, 3);
  autopas::LinkedCells<SPHParticle, autopas::FullParticleCell<SPHParticle>> _linkedCells({0., 0., 0.}, {5., 5., 5.},
                                                                                         cutoff);

  autopas::sph::SPHParticle defaultSPHParticle({0., 0., 0.}, {1., .5, .25}, 1, 2.5, 0.7, 0.6);
  RandomGenerator::fillWithParticles(_verletLists, defaultSPHParticle, numMolecules);
  // now fill second container with the molecules from the first one, because
  // otherwise we generate new particles
  for (auto it = _verletLists.begin(); it.isValid(); ++it) {
    _linkedCells.addParticle(*it);
  }

  autopas::sph::SPHCalcDensityFunctor densityFunctor;
  _verletLists.iteratePairwiseSoA(&densityFunctor);
  _linkedCells.iteratePairwiseSoA(&densityFunctor);

  auto itDirect = _verletLists.begin();
  auto itLinked = _linkedCells.begin();

  std::vector<std::array<double, 3>> forcesDirect(numMolecules), forcesLinked(numMolecules);
  // get and sort by id, the
  for (auto it = _verletLists.begin(); it.isValid(); ++it) {
    SPHParticle &m = *it;
    forcesDirect.at(m.getID()) = m.getF();
  }

  for (auto it = _linkedCells.begin(); it.isValid(); ++it) {
    SPHParticle &m = *it;
    forcesLinked.at(m.getID()) = m.getF();
  }

  for (unsigned long i = 0; i < numMolecules; ++i) {
    for (int d = 0; d < 3; ++d) {
      double f1 = forcesDirect[i][d];
      double f2 = forcesLinked[i][d];
      EXPECT_NEAR(f1, f2, std::fabs(f1 * rel_err_tolerance));
    }
  }
}
