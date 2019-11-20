/**
 * @file TimeDiscretizationTest.h
 * @author N. Fottner
 * @date 05/22/19.
 */

#pragma once

#include <gtest/gtest.h>
#include <vector>

#include "AutoPasTestBase.h"
#include "PrintableMolecule.h"
#include "TimeDiscretization.h"
#include "autopas/AutoPas.h"
#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
#include "autopas/utils/ArrayMath.h"
#include "testingHelpers/GridGenerator.h"
#include "testingHelpers/RandomGenerator.h"
#include "testingHelpers/commonTypedefs.h"

class TimeDiscretizationTest : public AutoPasTestBase {
 public:
  TimeDiscretizationTest()
      : AutoPasTestBase() {
    functor.setParticleProperties(1., 1.);
    _particlePropertiesLibrary.addType(0, 1, 1, 1);
  }
  void fillWithParticlesAndInit(autopas::AutoPas<PrintableMolecule,
                                                 autopas::FullParticleCell<PrintableMolecule>> &autopas);

 protected:
  double cutoff{2.};
  ParticlePropertiesLibrary<double, size_t> _particlePropertiesLibrary;
  autopas::LJFunctor<PrintableMolecule, autopas::ParticleCell<PrintableMolecule>> functor{autopas::LJFunctor<PrintableMolecule, autopas::ParticleCell<PrintableMolecule>>(cutoff, 0.0)};
};