/**
* @file MulticenteredLJFunctorTest.h
* @author S. Newcome
* @date 12/05/2022
*/

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/molecularDynamics/MulticenteredMoleculeLJ.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"

struct multisiteMolecule {
  std::array<double,3> CoMPosition;
  std::array<double,4> quaternion;
  std::array<double,3> force;
  std::array<double,3> torque;
  std::vector<std::array<double,3>> sitePositions;
  std::vector<double> siteEpsilons;
  std::vector<double> siteSigmas;
};

/**
 * Test class for MulticenteredLJFunctor
 */
class MulticenteredLJFunctorTest : public AutoPasTestBase {
 public:
  /**
   * Constructor
   */
  MulticenteredLJFunctorTest() = default;

  void testForceCalculation(multisiteMolecule molA, multisiteMolecule molB, double cutoff);
};


