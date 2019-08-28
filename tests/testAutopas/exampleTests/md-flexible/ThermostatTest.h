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
      : AutoPasTestBase(),
        _autopas(autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>()),
        _particlePropertiesLibrary(ParticlePropertiesLibrary<double, size_t>()),
        _timeDiscretization(
            TimeDiscretization<decltype(_autopas), std::remove_reference_t<decltype(_particlePropertiesLibrary)>>(
                0.002 /*delta_t*/, _particlePropertiesLibrary)) {
    _particlePropertiesLibrary.addType(0, 1., 1., 1.); /*initializing the default particlePropertiesLibrary*/
  }

  void initFillWithParticles(std::array<unsigned long, 3> particlesPerDim, double particleSpacing, double cutoff);

 protected:
  autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> _autopas;
  ParticlePropertiesLibrary<double, size_t> _particlePropertiesLibrary;
  TimeDiscretization<decltype(_autopas), std::remove_reference_t<decltype(_particlePropertiesLibrary)>>
      _timeDiscretization;
};
