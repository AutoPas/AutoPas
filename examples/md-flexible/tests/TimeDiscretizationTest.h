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
      : AutoPasTestBase(),
        cutoff{1.},
        boxmin{{0., 0., 0.}},
        boxmax{{5., 5., 5.}},
        _particlePropertiesLibrary(),
        functor{autopas::LJFunctor<PrintableMolecule, autopas::ParticleCell<PrintableMolecule>>(cutoff, 0.0)} {
    functor.setParticleProperties(1., 1.);
    _particlePropertiesLibrary.addType(0, 1, 1, 1);
  }
  void globalForceTest(autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> &auto1,
                       autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> &auto2,
                       int iterations);
  void initFillWithParticles(autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> &autopas,
                             std::array<unsigned long, 3> particlesPerDim);

  static std::array<double, 3> nextPosition(std::array<double, 3> position, std::array<double, 3> force,
                                            std::array<double, 3> velocity, double particle_delta_t);
  static std::array<double, 3> nextVelocity(std::array<double, 3> velocity, std::array<double, 3> force,
                                            std::array<double, 3> oldf, double particle_delta_t);

  void Pos_and_Velo_Test(autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> &autopas,
                         size_t numberOfParticles, int iterations);

 protected:
  double cutoff;
  std::array<double, 3> boxmin;
  std::array<double, 3> boxmax;
  ParticlePropertiesLibrary<double, size_t> _particlePropertiesLibrary;
  autopas::LJFunctor<PrintableMolecule, autopas::ParticleCell<PrintableMolecule>> functor;
};