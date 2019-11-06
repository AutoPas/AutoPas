/**
 * @file PCLTest.h
 * @author N. Fottner
 * @date 7/7/19
 */
#pragma once
#include <gtest/gtest.h>
#include "PrintableMolecule.h"
#include "Simulation.h"
#include "parsing/YamlParser.h"
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
#include "testingHelpers/commonTypedefs.h"

class ParticlePropertiesLibraryTest : public AutoPasTestBase {
 public:
  ParticlePropertiesLibraryTest()
      : AutoPasTestBase(), epsilon{1.0}, sigma{1.0}, epsilon2{2.0}, sigma2{2.0}, mass{1.0} {}

  static double mixingE(double e1, double e2);
  static double mixingS(double s1, double s2);

 protected:
  double epsilon;
  double sigma;
  double epsilon2;
  double sigma2;
  double mass;
  ParticlePropertiesLibrary<double, size_t> PPL;
};