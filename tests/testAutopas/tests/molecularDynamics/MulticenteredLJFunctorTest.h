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
   * Tests the correctness of the AoS functor for a given molA, molB, PPL, and cutoff.
   * @tparam newton3
   * @param molA
   * @param molB
   * @param cutoff
   */
  template <bool newton3>
  void testAoSForceCalculation(autopas::MulticenteredMoleculeLJ molA, autopas::MulticenteredMoleculeLJ molB, ParticlePropertiesLibrary<double, size_t> PPL, double cutoff);
};


