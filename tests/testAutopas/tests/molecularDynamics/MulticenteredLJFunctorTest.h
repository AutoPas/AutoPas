/**
* @file MulticenteredLJFunctorTest.h
* @author S. Newcome
* @date 12/05/2022
*/

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/AutoPasDecl.h"
#include "autopas/molecularDynamics/LJFunctor.h"
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
   * Generates a reproducible set of molecules of varying numbers of molecules and a particle property library.
   * @param molecules
   * @param PPL
   */
  void generateMoleculesAndPPL(std::vector<autopas::MulticenteredMoleculeLJ> *molecules, ParticlePropertiesLibrary<double, size_t> *PPL);

  /**
   * Tests the correctness of the AoS functor for a given molA, molB, PPL, and cutoff.
   * @tparam newton3
   * @param molA
   * @param molB
   * @param PPL
   * @param cutoff
   */
  template <bool newton3>
  void testAoSForceCalculation(autopas::MulticenteredMoleculeLJ molA, autopas::MulticenteredMoleculeLJ molB, ParticlePropertiesLibrary<double, size_t> PPL, double cutoff);

  /**
   * Compares the forces produced by the single-site functor and the multi-site functor applied to single-site molecules.
   * @tparam newton3
   * @param molA
   * @param molB
   * @param PPL
   * @param cutoff
   */
  template <bool newton3>
  void singleSiteSanityCheck(autopas::MulticenteredMoleculeLJ molA, autopas::MulticenteredMoleculeLJ molB, ParticlePropertiesLibrary<double, size_t> PPL, double cutoff);

  /**
   * Compares the correctness of the SoACellFunctor against that of the AoSFunctor.
   * @tparam newton3
   * @param molecules vector of multi-site molecules
   * @param PPL Particle Properties Library
   * @param cutoff
   */
  template <bool newton3>
  void testSoACellAgainstAoS(std::vector<autopas::MulticenteredMoleculeLJ> molecules, ParticlePropertiesLibrary<double, size_t> PPL, double cutoff);

  /**
   * Compares the correctness of the SoACellPairFunctor against that of the AoSFunctor.
   * @tparam newton3
   * @param moleculesA
   * @param moleculesB
   * @param PPL
   * @param cutoff
   */
  template <bool newton3>
  void testSoACellPairAgainstAoS(std::vector<autopas::MulticenteredMoleculeLJ> moleculesA, std::vector<autopas::MulticenteredMoleculeLJ> moleculesB, ParticlePropertiesLibrary<double, size_t> PPL, double cutoff);
};


