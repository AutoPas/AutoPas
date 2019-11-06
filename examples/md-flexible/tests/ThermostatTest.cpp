/**
 * @file ThermostatTest.cpp
 * @author N. Fottner
 * @date 28/08/19.
 */

#include "ThermostatTest.h"

void ThermostatTest::initContainer(AutoPasType &autopas, const Molecule &dummy, std::array<size_t, 3> particlesPerDim) {
  constexpr double particleSpacing = 1.;
  constexpr double cutoff = 1.;

  double minimalBoxLength = cutoff + autopas.getVerletSkin();
  std::array<double, 3> boxmax = {std::max(particlesPerDim[0] * particleSpacing, minimalBoxLength),
                                  std::max(particlesPerDim[1] * particleSpacing, minimalBoxLength),
                                  std::max(particlesPerDim[2] * particleSpacing, minimalBoxLength)};
  autopas.setBoxMax(boxmax);
  autopas.setBoxMin({0., 0., 0.});
  autopas.setCutoff(cutoff);
  autopas.init();
  // place particles grid in the middle of the domain
  GridGenerator::fillWithParticles(autopas, particlesPerDim, dummy, {particleSpacing, particleSpacing, particleSpacing},
                                   {particleSpacing / 2, particleSpacing / 2, particleSpacing / 2});
}

void ThermostatTest::testBrownianMotion(const Molecule &dummyMolecule, bool useCurrentTemp) {
  initContainer(_autopas, dummyMolecule, {2, 1, 1});

  const double tInit = 0.01;
  const double tTarget = 0.1;
  const double deltaTemp = 0.01;

  ThermostatType thermostat(tInit, tTarget, deltaTemp, _particlePropertiesLibrary);

  thermostat.addBrownianMotion(_autopas, useCurrentTemp);
  // check that velocities have actually changed
  for (auto iter = _autopas.begin(); iter.isValid(); ++iter) {
    for (int i = 0; i < 3; ++i) {
      EXPECT_THAT(iter->getV()[i], ::testing::Not(::testing::DoubleEq(dummyMolecule.getV()[i])));
    }
  }
}

TEST_F(ThermostatTest, BrownianMotionTest_useCurrentTempFalse) { testBrownianMotion(Molecule(), false); }

TEST_F(ThermostatTest, BrownianMotionTest_useCurrentTempTrue) {
  Molecule m;
  m.setV({1., 1., 1.});
  testBrownianMotion(m, true);
}

// void ThermostatTest::basicApplication(double initT, double targetT, double deltaT, bool initBM, AutoPasType &autopas)
// {
//  double nrApplications = targetT / deltaT;
//  auto thermostat = Thermostat<decltype(autopas), std::remove_reference_t<decltype(_particlePropertiesLibrary)>>(
//      initT, targetT, deltaT, _particlePropertiesLibrary);
//  if (initBM) {
//    thermostat.addBrownianMotion(autopas, false);
//  } else {
//    // initial velocity value of particles necessary otherwise zero divisions causing error
//    for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
//      iter->addV({0, 0.1, 0});
//    }
//    thermostat.addBrownianMotion(autopas, false);
//  }
//  thermostat.apply(autopas);
//  for (size_t i = 1; i <= nrApplications; i++) {
//    EXPECT_NEAR(thermostat.calcTemperature(autopas), (initT + (deltaT * i)) > targetT ? targetT : initT + (deltaT *
//    i),
//                _absDelta);
//    thermostat.apply(autopas);
//  }
//  thermostat.apply(autopas);
//  EXPECT_NEAR(thermostat.calcTemperature(autopas), targetT, _absDelta);
//}
//
// void ThermostatTest::calcTemperature(size_t particlesPerDimension) {
//  auto autopas = autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>();
//  initContainer(1. /*particleSpacing*/, 1.5 /*cutoff*/,
//                {particlesPerDimension, particlesPerDimension, particlesPerDimension} /*particlesPerdim*/);
//  auto _thermostat = Thermostat<decltype(autopas), std::remove_reference_t<decltype(_particlePropertiesLibrary)>>(
//      0.1 /*initT*/, 5 /*targetT*/, 0.01 /*deltaT*/, _particlePropertiesLibrary);
//  for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
//    iter->setV({0.1, 0.1, 0.2});
//  }
//  EXPECT_NEAR(_thermostat.calcTemperature(autopas), 0.02, _absDelta);
//}
//
// TEST_F(ThermostatTest, Application) {
//  auto _autopas = autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>();
//  initContainer(1. /*particleSpacing*/, 1.5 /*cutoff*/, {2, 2, 2} /*particlesPerdim*/);
//  basicApplication(0.01, 3.0, 0.001, true, _autopas);
//  // reset _autopas
//  initContainer(1. /*particleSpacing*/, 1.5 /*cutoff*/, {2, 2, 2} /*particlesPerdim*/);
//  basicApplication(1.0, 7.0, 0.01, true, _autopas);
//  // dont use current Temp for Brownian Motion initialization
//  initContainer(1. /*particleSpacing*/, 1.5 /*cutoff*/, {2, 2, 2} /*particlesPerdim*/);
//  basicApplication(0.01, 3.0, 0.001, false, _autopas);
//  // reset _autopas
//  initContainer(1. /*particleSpacing*/, 1.5 /*cutoff*/, {2, 2, 2} /*particlesPerdim*/);
//  basicApplication(1.0, 7.0, 0.01, false, _autopas);
//}
//
// TEST_F(ThermostatTest, calcTemperature) {
//  calcTemperature(1);
//  calcTemperature(2);
//  calcTemperature(10);
//}