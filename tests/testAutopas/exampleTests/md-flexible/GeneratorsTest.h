
#pragma once
#include <gtest/gtest.h>
#include <math.h>
#include <vector>
#include "../../../../examples/md-flexible/ParticleClassLibrary.h"
#include "../../../../examples/md-flexible/PrintableMolecule.h"
#include "../../../../examples/md-flexible/TimeDiscretization.h"
#include "../../../../src/autopas/utils/ArrayMath.h"
#include "../../testingHelpers/GaussianGenerator.h"
#include "../../testingHelpers/GridGenerator.h"
#include "../../testingHelpers/RandomGenerator.h"
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "testingHelpers/commonTypedefs.h"

class GeneratorsTest : public AutoPasTestBase {
 public:
  GeneratorsTest()
      : AutoPasTestBase(),
        epsilon{1.0},
        sigma{1.0},
        cutoff{1.},
        boxmin{{0., 0., 0.}},
        boxmax{{5., 5., 5.}},
        PCL{ParticleClassLibrary(epsilon, sigma, 1.0, 800)}
  // functor{autopas::LJFunctor<Particle, FPCell,
  //                           autopas::FunctorN3Modes::Both, true>(cutoff, PCL, 0.0)}
  {}

  void MolSimTaskGeneration(autopas::AutoPas<Particle, FPCell> &autopas);

 protected:
  double epsilon;
  double sigma;
  double cutoff;
  array<double, 3> boxmin;
  array<double, 3> boxmax;
  ParticleClassLibrary PCL;
  // autopas::LJFunctor<Particle, FPCell, autopas::FunctorN3Modes::Both, true> functor;
};