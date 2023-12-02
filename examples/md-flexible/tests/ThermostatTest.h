/**
 * @file ThermostatTest.h
 * @author N. Fottner
 * @date 28/08/19.
 */
#pragma once
#include "AutoPasTestBase.h"
#include "autopas/AutoPasDecl.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"
#include "src/TypeDefinitions.h"
#include "testingHelpers/commonTypedefs.h"

extern template class autopas::AutoPas<ParticleType>;

class ThermostatTest : public AutoPasTestBase,
                       public ::testing::WithParamInterface<std::tuple<double, double, double>> {
 public:
  using AutoPasType = autopas::AutoPas<ParticleType>;

  ThermostatTest() : AutoPasTestBase(), _particlePropertiesLibrary(ParticlePropertiesLibrary<double, size_t>(1.)) {
    _particlePropertiesLibrary.addSiteType(0, 1.);
    _particlePropertiesLibrary.addLJSite(0, 1., 1.);
    _particlePropertiesLibrary.addSiteType(1, 2.);
    _particlePropertiesLibrary.addLJSite(1, 1., 1.);

#if MD_FLEXIBLE_MODE == MULTISITE
    _particlePropertiesLibrary.addMolType(0, {0}, {{0., 0., 0.}}, {1., 1., 1.});
    _particlePropertiesLibrary.addMolType(1, {0, 0, 1}, {{0., -0.05, 0.}, {0.5, 0., 0.}, {0., 0.25, 0.25}},
                                          {1., 1., 1.});
#endif
    _particlePropertiesLibrary.calculateMixingCoefficients();
  }

 protected:
  /**
   * Fills an autopas container with a given grid of copies of the dummy particle and initializes the container.
   * @param autopas
   * @param dummy
   * @param particlesPerDim
   */
  void initContainer(AutoPasType &autopas, const ParticleType &dummy, std::array<size_t, 3> particlesPerDim);

  /**
   * Applies brownian motion to a system and checks that all velocities have changed.
   * @param dummyMolecule
   * @param useCurrentTemp
   */
  void testBrownianMotion(const ParticleType &dummyMolecule, double targetTemperature);

  ParticlePropertiesLibrary<double, size_t> _particlePropertiesLibrary;
  AutoPasType _autopas;
};
