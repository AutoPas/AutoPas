/**
 * @file ThermostatTest.h
 * @author N. Fottner
 * @date 28/08/19.
 */
#pragma once
#include "Generator.h"
#include "PrintableMolecule.h"
#include "Thermostat.h"
#include "TimeDiscretization.h"
#include "AutoPasTestBase.h"
#include "Objects/Objects.h"
#include "autopas/AutoPas.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
#include "autopas/utils/ArrayUtils.h"
#include "testingHelpers/commonTypedefs.h"

class ThermostatTest : public AutoPasTestBase {
 public:
  using AutoPasType = autopas::AutoPas<Molecule, autopas::FullParticleCell<Molecule>>;
  using ThermostatType = Thermostat<AutoPasType, ParticlePropertiesLibrary<double, size_t>>;

  ThermostatTest() : AutoPasTestBase(), _particlePropertiesLibrary(ParticlePropertiesLibrary<double, size_t>()) {
    _particlePropertiesLibrary.addType(0, 1., 1., 1.); /*initializing the default particlePropertiesLibrary*/
  }

  static void initContainer(AutoPasType &autopas, const Molecule &dummy, std::array<size_t, 3> particlesPerDim);

  void basicApplication(double initT, double targetT, double deltaT, bool initBM, AutoPasType &autopas);

  void calcTemperature(size_t particlesPerDimension);

 protected:

  void testBrownianMotion(const Molecule &dummyMolecule, bool useCurrentTemp);

  ParticlePropertiesLibrary<double, size_t> _particlePropertiesLibrary;
  AutoPasType _autopas;
  static constexpr double _absDelta = 1e-7;
};
