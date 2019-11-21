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

//  ThermostatType thermostat(1, 1, 1, _particlePropertiesLibrary);

  Thermostat::addBrownianMotion(_autopas, _particlePropertiesLibrary, useCurrentTemp);
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

/**
 * Init a system with 2x2x2 particles. Calculate its temperature and heat it up from 0 to 2.
 */
TEST_F(ThermostatTest, ApplyAndCalcTemperatureTest) {
  Molecule m;
  initContainer(_autopas, m, {2, 2, 2});
  _particlePropertiesLibrary.addType(0, 1., 1., 1.);
//  ThermostatType thermo(init 1.,target  2., delta .5, _particlePropertiesLibrary);
  constexpr double initTemp = 1.;
  constexpr double targetTemp = 2.;
  constexpr double deltaTemp = .5;
  EXPECT_NEAR(Thermostat::calcTemperature(_autopas, _particlePropertiesLibrary), 0, 1e-12);
  // add random velocities so that we do not scale zero vectors
  Thermostat::addBrownianMotion(_autopas, _particlePropertiesLibrary, false);
  // expect temperature to have changed from zero
  EXPECT_THAT(Thermostat::calcTemperature(_autopas, _particlePropertiesLibrary), ::testing::Not(::testing::DoubleNear(0, 1e-12)));
  // set system to initial temperature
  Thermostat::apply(_autopas, _particlePropertiesLibrary, initTemp, std::numeric_limits<double>::max());
  EXPECT_NEAR(Thermostat::calcTemperature(_autopas, _particlePropertiesLibrary), 1., 1e-12);
  Thermostat::apply(_autopas, _particlePropertiesLibrary, targetTemp, deltaTemp);
  EXPECT_NEAR(Thermostat::calcTemperature(_autopas, _particlePropertiesLibrary), 1.5, 1e-12);
  Thermostat::apply(_autopas, _particlePropertiesLibrary, targetTemp, deltaTemp);
  EXPECT_NEAR(Thermostat::calcTemperature(_autopas, _particlePropertiesLibrary), 2., 1e-12);
  Thermostat::apply(_autopas, _particlePropertiesLibrary, targetTemp, deltaTemp);
  EXPECT_NEAR(Thermostat::calcTemperature(_autopas, _particlePropertiesLibrary), 2., 1e-12);
}
