/**
 * @file LJFunctorAVX512Test.h
 * @author S. Newcome
 * @date 23/01/2024
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/AutoPasDecl.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "molecularDynamicsLibrary/LJFunctorAVX512_Mask.h"
#include "molecularDynamicsLibrary/MoleculeLJ.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"

/**
 * Test class for LJFunctorAVX512
 */
class LJFunctorAVX512Test : public AutoPasTestBase {
 public:
  /**
   * Constructor
   */
  LJFunctorAVX512Test() = default;

  /**
   * Generates a particle properties library for use with ::generateMolecules
   * @param PPL
   */
  void generatePPL(ParticlePropertiesLibrary<double, size_t> *PPL);
  /**
   * Generates a reproducible set of molecules.
   * @param molecules
   * @param offset position offset
   * @param allOwned true if all generated molecules are owned. If false, molecule ownership will alternate in the order
   * owned -> halo -> dummy -> owned -> ...
   */
  void generateMolecules(std::vector<mdLib::MoleculeLJ> *molecules, std::array<double, 3> offset,
                         bool allOwned);

  /**
   * Tests the correctness of the AoS functor for a given molA, molB, PPL, and cutoff.
   * @tparam functorType Functor tested.
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
  template <template<class, bool, bool, autopas::FunctorN3Modes, bool, bool> class functorType, bool newton3, bool calculateGlobals, bool applyShift>
  void testAoSForceCalculation(mdLib::MoleculeLJ molA, mdLib::MoleculeLJ molB,
                               ParticlePropertiesLibrary<double, size_t> &PPL, double cutoff);

  /**
   * Runs a suite of testAoSForceCalculation functions such that the AoS Functor is tested for the given arguments with
   * all combinations of
   * * with and without Newton's 3rd law optimization.
   * * with and without the calculation of global attributes i.e. potential energy and virial
   * * if calculating global attributes, with and without applying a shift to the potential energy.
   * @tparam functorType Functor tested. Only valid for functors which act on mdLib::MoleculeLJ applying LJ forces
   * @param molA
   * @param molB
   * @param PPL
   * @param cutoff
   */
  template <template<class, bool, bool, autopas::FunctorN3Modes, bool, bool> class functorType>
  void testSuiteAoSForceCalculation(mdLib::MoleculeLJ molA, mdLib::MoleculeLJ molB,
                                    ParticlePropertiesLibrary<double, size_t> &PPL, double cutoff);

  /**
   * Carries out AoS Dummy Test for some given functor F. Interacts 3 dummy molecules and expects no forces.
   * @tparam F functor to be tested.
   */
  template <class F>
  void AoSDummyTestHelper();

  /**
   * Compares the correctness of the SoACellFunctor against that of the AoSFunctor.
   * @tparam functorType Functor tested. Only valid for functors which act on mdLib::MoleculeLJ, applying
   * Lennard-Jones forces.
   * @tparam newton3 enables N3L optimizations for both functors.
   * @tparam calculateGlobals Both functors calculate global attributes such as potential energy and virial which are
   * compared in addition to force and torque.
   * @tparam applyShift applies a shift to the potential energy calculation for both functors.
   * @param molecules vector of single-site molecules
   * @param PPL Particle Properties Library
   * @param cutoff
   */
  template <template<class, bool, bool, autopas::FunctorN3Modes, bool, bool> class functorType, bool newton3, bool calculateGlobals, bool applyShift>
  void testSoACellAgainstAoS(std::vector<mdLib::MoleculeLJ> molecules,
                             ParticlePropertiesLibrary<double, size_t> &PPL, double cutoff);

  /**
   * Compares the correctness of the SoACellPairFunctor against that of the AoSFunctor.
   * @tparam functorType Functor tested. Only valid for functors which act on mdLib::MoleculeLJ, applying
   * Lennard-Jones forces.
   * @tparam newton3 enables N3L optimizations for both functors.
   * @tparam calculateGlobals Both functors calculate global attributes such as potential energy and virial which are
   * compared in addition to force and torque.
   * @tparam applyShift applies a shift to the potential energy calculation for both functors.
   * @param moleculesA vector of single-site molecules
   * @param moleculesB vector of single-site molecules
   * @param PPL
   * @param cutoff
   */
  template <template<class, bool, bool, autopas::FunctorN3Modes, bool, bool> class functorType, bool newton3, bool calculateGlobals, bool applyShift>
  void testSoACellPairAgainstAoS(std::vector<mdLib::MoleculeLJ> moleculesA,
                                 std::vector<mdLib::MoleculeLJ> moleculesB,
                                 ParticlePropertiesLibrary<double, size_t> &PPL, double cutoff);

  /**
   * Compares the correctness of the SoAVerletFunctor against that of the AoSFunctor.
   * @tparam functorType Functor tested. Only valid for functors which act on mdLib::MoleculeLJ, applying
   * Lennard-Jones forces.
   * @tparam newton3 enables N3L optimizations for both functors.
   * @tparam calculateGlobals Both functors calculate global attributes such as potential energy and virial which are
   * compared in addition to force and torque.
   * @tparam applyShift applies a shift to the potential energy calculation for both functors.
   * @param molecules vector of single-site molecules
   * @param PPL
   * @param cutoff
   */
  template <template<class, bool, bool, autopas::FunctorN3Modes, bool, bool> class functorType, bool newton3, bool calculateGlobals, bool applyShift>
  void testSoAVerletAgainstAoS(std::vector<mdLib::MoleculeLJ> molecules,
                               ParticlePropertiesLibrary<double, size_t> &PPL, double cutoff);
};
