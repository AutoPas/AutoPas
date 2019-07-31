#pragma once
#include <gtest/gtest.h>
#include "../../../../examples/md-flexible/PrintableMolecule.h"
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
#include "testingHelpers/commonTypedefs.h"
using namespace std;

class PCLTest : public AutoPasTestBase {
 public:
  PCLTest() : AutoPasTestBase(), epsilon{1.0}, sigma{1.0}, epsilon2{2.0}, sigma2{2.0}, mass{1.0} {}

  static double mixingE(double e1, double e2);
  static double mixingS(double s1, double s2);

 protected:
  double epsilon;
  double sigma;
  double epsilon2;
  double sigma2;
  double mass;
  ParticlePropertiesLibrary PPL;
};