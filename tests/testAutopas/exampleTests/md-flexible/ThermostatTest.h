/**
 * @file ThermostatTest.h
 * @author N. Fottner
 * @date 28/08/19.
 */
#pragma once
#include "../../../../examples/md-flexible/Generator.h"
#include "../../../../examples/md-flexible/Objects.h"
#include "../../../../examples/md-flexible/PrintableMolecule.h"
#include "../../../../examples/md-flexible/Thermostat.h"
#include "../../../../examples/md-flexible/TimeDiscretization.h"
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
#include "autopas/utils/ArrayUtils.h"

class ThermostatTest : public AutoPasTestBase {
 public:
  ThermostatTest()
      : AutoPasTestBase(), _particlePropertiesLibrary(ParticlePropertiesLibrary<double, size_t>()), absDelta(1e-7) {
    _particlePropertiesLibrary.addType(0, 1., 1., 1.); /*initializing the default particlePropertiesLibrary*/
  }

  static void initFillWithParticles(
      std::array<unsigned long, 3> particlesPerDim, double particleSpacing, double cutoff,
      autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> &autopas);

  void basicApplication(double initT, double targetT, double deltaT, bool initBM,
                        autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> &autopas);

  void calcTemperature(size_t particlesPerDimension);

 protected:
  ParticlePropertiesLibrary<double, size_t> _particlePropertiesLibrary;
  double absDelta;
};
