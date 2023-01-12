/**
 * @file ThermostatTest.cpp
 * @author N. Fottner
 * @date 28/08/19.
 */

#include "ThermostatTest.h"

#include "autopasTools/generators/GridGenerator.h"
#include "src/Thermostat.h"

void ThermostatTest::initContainer(AutoPasType &autopas, const ParticleType &dummy, std::array<size_t, 3> particlesPerDim) {
  constexpr double particleSpacing = 1.;
  constexpr double cutoff = 1.;

  double minimalBoxLength = cutoff + autopas.getVerletSkin();
  std::array<double, 3> boxMax = {std::max((double)particlesPerDim[0] * particleSpacing, minimalBoxLength),
                                  std::max((double)particlesPerDim[1] * particleSpacing, minimalBoxLength),
                                  std::max((double)particlesPerDim[2] * particleSpacing, minimalBoxLength)};
  autopas.setBoxMax(boxMax);
  autopas.setBoxMin({0., 0., 0.});
  autopas.setCutoff(cutoff);
  autopas.init();
  // place particles grid in the middle of the domain
  autopasTools::generators::GridGenerator::fillWithParticles(
      autopas, particlesPerDim, dummy, {particleSpacing, particleSpacing, particleSpacing},
      {particleSpacing / 2, particleSpacing / 2, particleSpacing / 2});
}

void ThermostatTest::testBrownianMotion(const ParticleType &dummyMolecule, const double targetTemperature) {
  // generate a significant number of molecules so that we can make statistical deductions
  initContainer(_autopas, dummyMolecule, {30, 30, 30});

  auto initTemperature = Thermostat::calcTemperature(_autopas, _particlePropertiesLibrary);

  // this only makes sense if the initial temperature is zero
  Thermostat::addBrownianMotion(_autopas, _particlePropertiesLibrary, targetTemperature);
  // check that velocities have actually changed
  for (auto iter = _autopas.begin(); iter.isValid(); ++iter) {
    for (int i = 0; i < 3; ++i) {
      EXPECT_THAT(iter->getV()[i], ::testing::Not(::testing::DoubleEq(dummyMolecule.getV()[i])));
    }
  }
  EXPECT_NEAR(Thermostat::calcTemperature(_autopas, _particlePropertiesLibrary), targetTemperature + initTemperature,
              0.02);
}

/**
 * Tests the application of Brownian Motion for the case that the molecule being applied to has no initial velocity.
 */
TEST_F(ThermostatTest, BrownianMotionTest_ZeroInitialVel) {
  ParticleType dummyMolecule;
  testBrownianMotion(dummyMolecule, 1.5);
}

/**
 * Tests the application of Brownian Motion for the case that the molecule being applied to has non-zero initial velocity.
 */
TEST_F(ThermostatTest, BrownianMotionTest_NonZeroInitialVel) {
  ParticleType dummyMolecule;
  dummyMolecule.setV({1, 1, 1});
  testBrownianMotion(dummyMolecule, 1.5);
}

/**
 * For multi-site molecules, tests the application of Brownian Motion for the case that the molecule being applied to
 * has non-zero initial angular velocity but zero translational velocity.
 */
TEST_F(ThermostatTest, BrownianMotionTest_NonZeroInitialAngVelWithZeroTransVel) {
#ifdef MD_FLEXIBLE_USE_MULTI_SITE
  ParticleType dummyMolecule;
  dummyMolecule.setAngularVel({1, 1, 1});
  testBrownianMotion(dummyMolecule, 1.5);
#endif
}

/**
 * For multi-site molecules, tests the application of Brownian Motion for the case that the molecule being applied to
 * has non-zero initial angular velocity and non-zero translational velocity.
 */
TEST_F(ThermostatTest, BrownianMotionTest_NonZeroInitialAngVelWithNonZeroTransVel) {
#ifdef MD_FLEXIBLE_USE_MULTI_SITE
  ParticleType dummyMolecule;
  dummyMolecule.setV({1, 1, 1});
  dummyMolecule.setAngularVel({1, 1, 1});
  testBrownianMotion(dummyMolecule, 1.5);
#endif
}

/**
 * Tests the application of Brownian Motion and the Thermostat for systems with multiple molecule types.
 */
TEST_F(ThermostatTest, MultiComponentTest) {
  ParticleType dummyMolecule;
  // fill with particles of type 0 and init
  initContainer(_autopas, dummyMolecule, {25, 25, 25});
  // add some type 1 particles
  dummyMolecule.setTypeId(1);
  autopasTools::generators::GridGenerator::fillWithParticles(_autopas, {25, 25, 25}, dummyMolecule);

  // init system with brownian motion and test for the given temperature
  constexpr double targetTemperature1 = 4.2;
  Thermostat::addBrownianMotion(_autopas, _particlePropertiesLibrary, targetTemperature1);
  auto temperatureMap = Thermostat::calcTemperatureComponent(_autopas, _particlePropertiesLibrary);

  for (auto &[typeId, temperature] : temperatureMap) {
    // brownian motion is probabilistic therefore a large tolerance is needed
    EXPECT_NEAR(temperature, targetTemperature1, 0.1);
  }

  // set system to a different temperature through apply and check that the temperature matches for each component
  constexpr double targetTemperature2 = targetTemperature1 + 2.;
  Thermostat::apply(_autopas, _particlePropertiesLibrary, targetTemperature2, std::numeric_limits<double>::max());
  temperatureMap = Thermostat::calcTemperatureComponent(_autopas, _particlePropertiesLibrary);
  for (auto &[typeId, temperature] : temperatureMap) {
    EXPECT_NEAR(temperature, targetTemperature2, 1e-9);
  }
}

/**
 * Initialize a system with the given temperature and scale it in steps of delta to the target temperature.
 * @param initialTemperature
 * @param targetTemperature
 * @param deltaTemperature
 */
TEST_P(ThermostatTest, testApplyAndCalcTemperature) {
  const double initialTemperature = std::get<0>(GetParam());
  const double targetTemperature = std::get<1>(GetParam());
  const double deltaTemperature = std::get<2>(GetParam());
  ParticleType m;
  initContainer(_autopas, m, {2, 2, 2});
  EXPECT_EQ(Thermostat::calcTemperature(_autopas, _particlePropertiesLibrary), 0)
      << "initially there are no velocities -> temperature should be exactly 0";
  // add random velocities so that we do not scale zero vectors
  Thermostat::addBrownianMotion(_autopas, _particlePropertiesLibrary, initialTemperature);
  EXPECT_THAT(Thermostat::calcTemperature(_autopas, _particlePropertiesLibrary),
              ::testing::Not(::testing::DoubleNear(0, 1e-12)))
      << "After Brownian motion expect temperature to have changed from zero";
  // set system to initial temperature
  Thermostat::apply(_autopas, _particlePropertiesLibrary, initialTemperature, std::numeric_limits<double>::max());
  EXPECT_NEAR(Thermostat::calcTemperature(_autopas, _particlePropertiesLibrary), initialTemperature, 1e-12)
      << "Thermostat failed to set initial temperature correctly!";

  const auto expectedIterations = std::ceil(std::abs((targetTemperature - initialTemperature) / deltaTemperature));
  for (int i = 1; i <= expectedIterations; ++i) {
    Thermostat::apply(_autopas, _particlePropertiesLibrary, targetTemperature, deltaTemperature);
    if (i != expectedIterations) {
      EXPECT_NEAR(Thermostat::calcTemperature(_autopas, _particlePropertiesLibrary),
                  initialTemperature + i * deltaTemperature, 1e-12);
    } else {
      EXPECT_NEAR(Thermostat::calcTemperature(_autopas, _particlePropertiesLibrary), targetTemperature, 1e-12);
    }
  }

  // apply once more to check that temperature stays the same
  Thermostat::apply(_autopas, _particlePropertiesLibrary, targetTemperature, deltaTemperature);
  EXPECT_NEAR(Thermostat::calcTemperature(_autopas, _particlePropertiesLibrary), targetTemperature, 1e-12)
      << "Thermostat changed the temperature when current temperature was already on target!";
}

/**
 * Aspects to consider:
 *  - increase temperature
 *  - decrease temperature
 *  - delta temperature hits target temperature exactly
 *  - delta temperature does not hit target temperature exactly
 */
INSTANTIATE_TEST_SUITE_P(Generated, ThermostatTest,
                         // tuple (initialTemperature, targetTemperature, deltaTemperature)
                         ::testing::Values(std::tuple(1., 2., .3), std::tuple(1., 2., .5), std::tuple(2., 1., -.3),
                                           std::tuple(2., 1., -.5)));
