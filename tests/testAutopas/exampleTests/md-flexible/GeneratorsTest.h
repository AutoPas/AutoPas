
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
  GeneratorsTest() : AutoPasTestBase() {
    epsilon = 5.0;
    sigma = 1.0;
    cutoff = 1;
    boxmin = {0., 0., 0.};
    boxmax = {5., 5., 5.};
    map<unsigned long, double> universalMap;
    for (unsigned long i = 0; i < 1000; i++) {
      universalMap.emplace(i, 1.0);
    }
    PCL = ParticleClassLibrary(universalMap, universalMap, universalMap);
    functor = new autopas::LJFunctor<PrintableMolecule, autopas::ParticleCell<PrintableMolecule>,
                                     autopas::FunctorN3Modes::Both, true>(cutoff, PCL, 0.0);
  }
  ~GeneratorsTest() { delete functor; }
  void MolSimTaskGeneration(autopas::AutoPas<Particle, FPCell> &autopas);

 protected:
  double epsilon;
  double sigma;
  double cutoff;
  array<double, 3> boxmin;
  array<double, 3> boxmax;
  ParticleClassLibrary PCL;
  autopas::LJFunctor<PrintableMolecule, autopas::ParticleCell<PrintableMolecule>, autopas::FunctorN3Modes::Both, true>
      *functor;
};