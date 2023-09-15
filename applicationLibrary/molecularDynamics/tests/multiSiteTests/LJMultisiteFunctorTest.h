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
#include "molecularDynamicsLibrary/LJMultisiteFunctorAVX.h"
#include "molecularDynamicsLibrary/LJMultisiteFunctorAVX_STS.h"
#include "molecularDynamicsLibrary/LJMultisiteFunctorAVX512.h"
#include "molecularDynamicsLibrary/MultisiteMoleculeLJ.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"

// todo there is a lot of code duplication for trying out multiple functors. We should condense into suites.

// Number of doubles that fit into a vector register
const size_t VecLength =
#ifdef __AVX512F__
8;
#else
4;
#endif

namespace {
// Some template aliases to make the tests cleaner.

template<class Particle, bool applyShift, bool useMixing, autopas::FunctorN3Modes useNewton3, bool calculateGlobals, bool relevantForTuning>
using LJMultisiteFunctorAVX_Masks = mdLib::LJMultisiteFunctorAVX<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning, true, VecLength>;

template<class Particle, bool applyShift, bool useMixing, autopas::FunctorN3Modes useNewton3, bool calculateGlobals, bool relevantForTuning>
using LJMultisiteFunctorAVX_GatherScatter = mdLib::LJMultisiteFunctorAVX<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning, false, VecLength>;

template<class Particle, bool applyShift, bool useMixing, autopas::FunctorN3Modes useNewton3, bool calculateGlobals, bool relevantForTuning>
using LJMultisiteFunctorAVX_STS = mdlib::LJMultisiteFunctorAVX_STS<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning, true, VecLength>;

template<class Particle, bool applyShift, bool useMixing, autopas::FunctorN3Modes useNewton3, bool calculateGlobals, bool relevantForTuning>
using LJMultisiteFunctorAVX512_Masks = mdLib::LJMultisiteFunctorAVX512<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning, true, VecLength>;
}

/**
 * Test class for LJMultisiteFunctor and LJMultisiteFunctorAVX
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
   * @param allOwned true if all generated molecules are owned. If false, molecule ownership will alternate in the order
   * owned -> halo -> dummy -> owned -> ...
   */
  void generateMolecules(std::vector<mdLib::MultisiteMoleculeLJ> *molecules, std::array<double, 3> offset,
                         bool allOwned);

  /**
   * Tests the correctness of the AoS functor for a given molA, molB, PPL, and cutoff. Only suitable for the CTC functors.
   * @tparam functorType Functor tested. Only valid for functors which act on mdLib::MultisiteMoleculeLJ, applying
   * Lennard-Jones forces using a cutoff criterion based on Center-of-Mass to Center-of-Mass cutoffs.
   * @tparam newton3 Applies N3L optimization to the AoS Functor.
   * @tparam calculateGlobals The AoS Functor calculates the potential energy and virial of the system of two molecules
   * and the results are tested against expected values.
   * @tparam applyShift The AoS Functor applies a shift to the potential energy such that it becomes continuous (and the
   * expected value is adjusted accordingly).
   * @tparam useMasks if true, functorType uses 0/1-masks instead of gather/scatter approach. Returns an error if functorType
   * @param molA
   * @param molB
   * @param PPL
   * @param cutoff
   */
  template <template<class, bool, bool, autopas::FunctorN3Modes, bool, bool> class functorType, bool newton3, bool calculateGlobals, bool applyShift>
  void testAoSForceCalculation_CTC(mdLib::MultisiteMoleculeLJ molA, mdLib::MultisiteMoleculeLJ molB,
                               ParticlePropertiesLibrary<double, size_t> PPL, double cutoff);

  /**
   * Runs a suite of testAoSForceCalculation functions such that the AoS Functor is tested for the given arguments with
   * all combinations of
   * * with and without Newton's 3rd law optimization.
   * * with and without the calculation of global attributes i.e. potential energy and virial
   * * if calculating global attributes, with and without applying a shift to the potential energy.
   * Only suitable for the CTC functors.
   * @tparam functorType Functor tested. Only valid for functors which act on mdLib::MultisiteMoleculeLJ, applying
   * Lennard-Jones forces using a cutoff criterion based on Center-of-Mass to Center-of-Mass cutoffs.
   * @param molA
   * @param molB
   * @param PPL
   * @param cutoff
   */
  template <template<class, bool, bool, autopas::FunctorN3Modes, bool, bool> class functorType>
  void testSuiteAoSForceCalculation_CTC(mdLib::MultisiteMoleculeLJ molA, mdLib::MultisiteMoleculeLJ molB,
                                    ParticlePropertiesLibrary<double, size_t> PPL, double cutoff);

  /**
   * Tests the correctness of the AoS functor for a given molA, molB, PPL, and cutoff. Only suitable for the STS functors.
   * @tparam functorType Functor tested. Only valid for functors which act on mdLib::MultisiteMoleculeLJ, applying
   * Lennard-Jones forces using a cutoff criterion based on site to site cutoffs.
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
  void testAoSForceCalculation_STS(mdLib::MultisiteMoleculeLJ molA, mdLib::MultisiteMoleculeLJ molB,
                                   ParticlePropertiesLibrary<double, size_t> PPL, double cutoff);

  /**
   * Runs a suite of testAoSForceCalculation functions such that the AoS Functor is tested for the given arguments with
   * all combinations of
   * * with and without Newton's 3rd law optimization.
   * * with and without the calculation of global attributes i.e. potential energy and virial
   * * if calculating global attributes, with and without applying a shift to the potential energy.
   * Only suitable for the STS functors.
   * @tparam functorType Functor tested. Only valid for functors which act on mdLib::MultisiteMoleculeLJ, applying
   * Lennard-Jones forces using a cutoff criterion based on site to site cutoffs.
   * @param molA
   * @param molB
   * @param PPL
   * @param cutoff
   */
  template <template<class, bool, bool, autopas::FunctorN3Modes, bool, bool> class functorType>
  void testSuiteAoSForceCalculation_STS(mdLib::MultisiteMoleculeLJ molA, mdLib::MultisiteMoleculeLJ molB,
                                        ParticlePropertiesLibrary<double, size_t> PPL, double cutoff);

  /**
   * Compares the forces produced by the single-site functor (LJFunctor) and a multi-site functor applied to single-site
   * molecules.
   *
   * Inputted molecules are objects of multi-site molecule class, but are expected to have molecule types which consist
   * of a single site at the center-of-mass with site type ID equal to molecule type ID.
   *
   * @tparam functorType Functor tested. Only valid for functors which act on mdLib::MultisiteMoleculeLJ, applying
   * Lennard-Jones forces.
   * @tparam newton3 applies Newton's 3rd Law optimization to both functors.
   * @tparam calculateGlobals both functors calculate global attributes such as potential energy and virial.
   * @tparam applyShift both functors add a shift to the potential energy.
   * @param molA single-site molecule represented by multi-site molecule class
   * @param molB single-site molecule represented by multi-site molecule class
   * @param PPL
   * @param cutoff
   */
  template <template<class, bool, bool, autopas::FunctorN3Modes, bool, bool> class functorType, bool newton3, bool calculateGlobals, bool applyShift>
  void singleSiteSanityCheck(mdLib::MultisiteMoleculeLJ molA, mdLib::MultisiteMoleculeLJ molB,
                             ParticlePropertiesLibrary<double, size_t> PPL, double cutoff);

  /**
   * Compares the correctness of the SoACellFunctor against that of the AoSFunctor.
   * @tparam functorType Functor tested. Only valid for functors which act on mdLib::MultisiteMoleculeLJ, applying
   * Lennard-Jones forces.
   * @tparam newton3 enables N3L optimizations for both functors.
   * @tparam calculateGlobals Both functors calculate global attributes such as potential energy and virial which are
   * compared in addition to force and torque.
   * @tparam applyShift applies a shift to the potential energy calculation for both functors.
   * @param molecules vector of multi-site molecules
   * @param PPL Particle Properties Library
   * @param cutoff
   */
  template <template<class, bool, bool, autopas::FunctorN3Modes, bool, bool> class functorType, bool newton3, bool calculateGlobals, bool applyShift>
  void testSoACellAgainstAoS(std::vector<mdLib::MultisiteMoleculeLJ> molecules,
                             ParticlePropertiesLibrary<double, size_t> PPL, double cutoff);

  /**
   * Compares the correctness of the SoACellPairFunctor against that of the AoSFunctor.
   * @tparam functorType Functor tested. Only valid for functors which act on mdLib::MultisiteMoleculeLJ, applying
   * Lennard-Jones forces using a cutoff criterion based on Center-of-Mass to Center-of-Mass cutoffs.
   * @tparam newton3 enables N3L optimizations for both functors.
   * @tparam calculateGlobals Both functors calculate global attributes such as potential energy and virial which are
   * compared in addition to force and torque.
   * @tparam applyShift applies a shift to the potential energy calculation for both functors.
   * @param moleculesA vector of multi-site molecules
   * @param moleculesB vector of multi-site molecules
   * @param PPL
   * @param cutoff
   */
  template <template<class, bool, bool, autopas::FunctorN3Modes, bool, bool> class functorType, bool newton3, bool calculateGlobals, bool applyShift>
  void testSoACellPairAgainstAoS(std::vector<mdLib::MultisiteMoleculeLJ> moleculesA,
                                 std::vector<mdLib::MultisiteMoleculeLJ> moleculesB,
                                 ParticlePropertiesLibrary<double, size_t> PPL, double cutoff);

  /**
   * Compares the correctness of the SoAVerletFunctor against that of the AoSFunctor.
   * @tparam functorType Functor tested. Only valid for functors which act on mdLib::MultisiteMoleculeLJ, applying
   * Lennard-Jones forces using a cutoff criterion based on Center-of-Mass to Center-of-Mass cutoffs.
   * @tparam newton3 enables N3L optimizations for both functors.
   * @tparam calculateGlobals Both functors calculate global attributes such as potential energy and virial which are
   * compared in addition to force and torque.
   * @tparam applyShift applies a shift to the potential energy calculation for both functors.
   * @param molecules vector of multi-site molecules
   * @param PPL
   * @param cutoff
   */
  template <template<class, bool, bool, autopas::FunctorN3Modes, bool, bool> class functorType, bool newton3, bool calculateGlobals, bool applyShift>
  void testSoAVerletAgainstAoS(std::vector<mdLib::MultisiteMoleculeLJ> molecules,
                               ParticlePropertiesLibrary<double, size_t> PPL, double cutoff);
};
