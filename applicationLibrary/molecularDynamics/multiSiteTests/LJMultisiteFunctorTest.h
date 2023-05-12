/**
* @file LJMultisiteFunctorTest.h
* @author S. Newcome
* @date 12/05/2022
*/

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/AutoPasDecl.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "molecularDynamicsLibrary/LJMultisiteFunctor.h"
#include "molecularDynamicsLibrary/MultisiteMoleculeLJ.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"

/**
 * Test class for LJMultisiteFunctor
 */
class LJMultisiteFunctorTest : public AutoPasTestBase {
 public:
  /**
   * Constructor
   */
  LJMultisiteFunctorTest() = default;

  /**
   * Generates a particle properties library for use with ::generateMolecules
   * @param PPL
   */
  void generatePPL(ParticlePropertiesLibrary<double, size_t> *PPL);
  /**
   * Generates a reproducible set of molecules of varying numbers of molecules.
   * @param molecules
   * @param offset position offset
   */
  void generateMolecules(std::vector<mdLib::MultisiteMoleculeLJ> *molecules, std::array<double, 3> offset);

  /**
   * Tests the correctness of the AoS functor for a given molA, molB, PPL, and cutoff.
   * @tparam newton3 Applies N3L optimization to the AoS Functor.
   * @tparam calculateGlobals The AoS Functor calculates the potential energy and virial of the system of two molecules
   * and the results are tested against expected values.
   * @tparam applyShift The AoS Functor applies a shift to the potential energy such that it becomes continuous (and the
   * expected value is adjusted accordingly).
   * @param molA
   * @param molB
   * @param PPL
   * @param cutoff
   */
  template <bool newton3, bool calculateGlobals, bool applyShift>
  void testAoSForceCalculation(mdLib::MultisiteMoleculeLJ molA, mdLib::MultisiteMoleculeLJ molB, ParticlePropertiesLibrary<double, size_t> PPL, double cutoff);

  /**
   * Runs a suite of testAoSForceCalculation functions such that the AoS Functor is tested for the given arguments with
   * all combinations of
   * * with and without Newton's 3rd law optimization.
   * * with and without the calculation of global attributes i.e. potential energy and virial
   * * if calculating global attributes, with and without applying a shift to the potential energy.
   * @param molA
   * @param molB
   * @param PPL
   * @param cutoff
   */
  void testSuiteAoSForceCalculation(mdLib::MultisiteMoleculeLJ molA, mdLib::MultisiteMoleculeLJ molB, ParticlePropertiesLibrary<double, size_t> PPL, double cutoff);

  /**
   * Compares the forces produced by the single-site functor and the multi-site functor applied to single-site molecules.
   *
   * Inputted molecules are objects of multi-site molecule class, but are expected to have molecule types which consist
   * of a single site at the center-of-mass with site type ID equal to molecule type ID.
   *
   * @tparam newton3 applies Newton's 3rd Law optimization to both functors.
   * @tparam calculateGlobals both functors calculate global attributes such as potential energy and virial.
   * @tparam applyShift both functors add a shift to the potential energy.
   * @param molA single-site molecule represented by multi-site molecule class
   * @param molB single-site molecule represented by multi-site molecule class
   * @param PPL
   * @param cutoff
   */
  template <bool newton3, bool calculateGobals, bool applyShift>
  void singleSiteSanityCheck(mdLib::MultisiteMoleculeLJ molA, mdLib::MultisiteMoleculeLJ molB, ParticlePropertiesLibrary<double, size_t> PPL, double cutoff);

  /**
   * Compares the correctness of the SoACellFunctor against that of the AoSFunctor.
   * @tparam newton3 enables N3L optimizations for both functors.
   * @tparam calculateGlobals Both functors calculate global attributes such as potential energy and virial which are
   * compared in addition to force and torque.
   * @tparam applyShift applies a shift to the potential energy calculation for both functors.
   * @param molecules vector of multi-site molecules
   * @param PPL Particle Properties Library
   * @param cutoff
   */
  template <bool newton3, bool calculateGlobals, bool applyShift>
  void testSoACellAgainstAoS(std::vector<mdLib::MultisiteMoleculeLJ> molecules, ParticlePropertiesLibrary<double, size_t> PPL, double cutoff);

  /**
   * Compares the correctness of the SoACellPairFunctor against that of the AoSFunctor.
   * @tparam newton3 enables N3L optimizations for both functors.
   * @tparam calculateGlobals Both functors calculate global attributes such as potential energy and virial which are
   * compared in addition to force and torque.
   * @tparam applyShift applies a shift to the potential energy calculation for both functors.
   * @param moleculesA vector of multi-site molecules
   * @param moleculesB vector of multi-site molecules
   * @param PPL
   * @param cutoff
   */
  template <bool newton3, bool calculateGlobals, bool applyShift>
  void testSoACellPairAgainstAoS(std::vector<mdLib::MultisiteMoleculeLJ> moleculesA, std::vector<mdLib::MultisiteMoleculeLJ> moleculesB, ParticlePropertiesLibrary<double, size_t> PPL, double cutoff);

  /**
   * Compares the correctness of the SoAVerletFunctor against that of the AoSFunctor.
   * @tparam newton3 enables N3L optimizations for both functors.
   * @tparam calculateGlobals Both functors calculate global attributes such as potential energy and virial which are
   * compared in addition to force and torque.
   * @tparam applyShift applies a shift to the potential energy calculation for both functors.
   * @param molecules vector of multi-site molecules
   * @param PPL
   * @param cutoff
   */
  template <bool newton3, bool calculateGlobals, bool applyShift>
  void testSoAVerletAgainstAoS(std::vector<mdLib::MultisiteMoleculeLJ> molecules, ParticlePropertiesLibrary<double, size_t> PPL, double cutoff);
};


