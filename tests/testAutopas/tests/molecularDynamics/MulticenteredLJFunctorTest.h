/**
* @file MulticenteredLJFunctorTest.h
* @author S. Newcome
* @date 12/05/2022
*/

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/AutoPasDecl.h"
#include "autopas/molecularDynamics/LJMulticenterFunctor.h"
#include "autopas/molecularDynamics/MulticenteredMoleculeLJ.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/Quaternion.h"

struct multisiteMolecule {
  std::array<double,3> CoMPosition;
  std::array<double,4> quaternion;
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

  /**
   * Tests the correctness of the AoS functor for a given molA, molB, and cutoff.
   * @tparam newton3
   * @param molA
   * @param molB
   * @param cutoff
   */
  template <bool newton3>
  void testAoSForceCalculation(multisiteMolecule molA, multisiteMolecule molB, double cutoff);
};


