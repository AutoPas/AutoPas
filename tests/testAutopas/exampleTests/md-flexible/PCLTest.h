#pragma once
#include <gtest/gtest.h>
#include "../../../../examples/md-flexible/ParticleClassLibrary.h"
#include "../../../../examples/md-flexible/PrintableMolecule.h"
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "testingHelpers/commonTypedefs.h"
using namespace std;

class PCLTest : public AutoPasTestBase {
 public:
  PCLTest() : AutoPasTestBase(),
    epsilon{1.0},
    sigma {1.0},
    mass {1.0},
    PCL { ParticleClassLibrary(epsilon,sigma, mass)}
    {}

  double mixingE(double e1, double e2);
  double mixingS(double s1, double s2);

 protected:
  PrintableMolecule dummyParticle;
  double epsilon;
  double sigma;
  double mass;
  ParticleClassLibrary PCL;
};