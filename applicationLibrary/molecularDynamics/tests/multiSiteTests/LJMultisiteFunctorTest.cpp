/**
 * @file LJMultisiteFunctorTest.cpp
 * @author S. Newcome
 * @date 16/05/2022
 */

#include "LJMultisiteFunctorTest.h"

#include <gtest/gtest.h>
#include "LJMultisiteFunctorTest.h"

#define PARTICLES_PER_DIM 8
#define AOS_VS_SOA_ACCURACY 1e-8


/**
 * Tests for the correctness of the AoS functor by applying to molecules designed to test all its functionality.
 */
TEST_F(LJMultisiteFunctorTest , AoSTest) {
  using mdLib::MultisiteMoleculeLJ;

  const double cutoff = 2.5;

  ParticlePropertiesLibrary PPL(cutoff);
  PPL.addSiteType(0, 1., 1., 1.);
  PPL.addSiteType(1, 0.5, 0.5, 0.5);

  // Molecules to be used in the tests (explanation of choices documented when tests are run).
  // For ease of readability, each molecule has its own molType, even when duplicated.
  MultisiteMoleculeLJ mol0;
  mol0.setR({0., 0., 0.});
  mol0.setQuaternion({0., 0., 0., 1.});
  mol0.setF({0., 0., 0.});
  mol0.setTorque({0., 0., 0.});
  mol0.setTypeId(0);
  PPL.addMolType(0, {0}, {{0., 0., 0.}}, {1., 1., 1.});

  MultisiteMoleculeLJ mol1;
  mol1.setR({1., 0., 0.});
  mol1.setQuaternion({0., 0., 0., 1.});
  mol1.setF({0., 0., 0.});
  mol1.setTorque({0., 0., 0.});
  mol1.setTypeId(1);
  PPL.addMolType(1, {0, 0}, {{0., 0.01, 0.}, {0., -0.01, 0.}}, {1., 1., 1.});

  MultisiteMoleculeLJ mol2;
  mol2.setR({0., 1., 0.});
  mol2.setQuaternion({0., 0., 0., 1.});
  mol2.setF({0., 0., 0.});
  mol2.setTorque({0., 0., 0.});
  mol2.setTypeId(2);
  PPL.addMolType(2, {0, 0}, {{0., 0.01, 0.}, {0., -0.01, 0.}}, {1., 1., 1.});

  MultisiteMoleculeLJ mol3;
  mol3.setR({0., 0., 0.});
  mol3.setQuaternion({0., 0., 0., 1.});
  mol3.setF({0., 0., 0.});
  mol3.setTorque({0., 0., 0.});
  mol3.setTypeId(3);
  PPL.addMolType(3, {0, 0, 0}, {{-0.05, -0.05, 0.}, {0., 0.1, 0.}, {0.05, -0.05, 0.}}, {1., 1., 1.});

  MultisiteMoleculeLJ mol4;
  mol4.setR({0., 0., 1.});
  mol4.setQuaternion({0.7071067811865475, 0.7071067811865475, 0., 0.});
  mol4.setF({0., 0., 0.});
  mol4.setTorque({0., 0., 0.});
  mol4.setTypeId(4);
  PPL.addMolType(4, {0, 0, 0}, {{-0.05, -0.05, 0.}, {0., 0.1, 0.}, {0.05, -0.05, 0.}}, {1., 1., 1.});

  MultisiteMoleculeLJ mol5;
  mol5.setR({2.5, 2.5, 2.5});
  mol5.setQuaternion({0., 0., 0., 1.});
  mol5.setF({0., 0., 0.});
  mol5.setTorque({0., 0., 0.});
  mol5.setTypeId(5);
  PPL.addMolType(5, {0, 0}, {{0., 0.01, 0.}, {0., -0.01, 0.}}, {1., 1., 1.});

  MultisiteMoleculeLJ mol6;
  mol6.setR({0., 2.55, 0.});
  mol6.setQuaternion({0., 0., 0., 1.});
  mol6.setF({0., 0., 0.});
  mol6.setTorque({0., 0., 0.});
  mol6.setTypeId(6);
  PPL.addMolType(6, {0, 0}, {{0., 0.1, 0.}, {0., -0.1, 0.}}, {1., 1., 1.});

  MultisiteMoleculeLJ mol7;
  mol7.setR({0., 2.45, 0.});
  mol7.setQuaternion({0., 0., 0., 1.});
  mol7.setF({0., 0., 0.});
  mol7.setTorque({0., 0., 0.});
  mol7.setTypeId(7);
  PPL.addMolType(7, {0, 0}, {{0., 0.1, 0.}, {0., -0.1, 0.}}, {1., 1., 1.});

  MultisiteMoleculeLJ mol8;
  mol8.setR({0., 0., 0.});
  mol8.setQuaternion({0., 0., 0., 1.});
  mol8.setF({0., 0., 0.});
  mol8.setTorque({0., 0., 0.});
  mol8.setTypeId(8);
  PPL.addMolType(8, {1, 1, 0}, {{-0.05, -0.05, 0.}, {0., 0.1, 0.}, {0.05, -0.05, 0.}}, {1., 1., 1.});

  PPL.calculateMixingCoefficients();

  // tests: 1 site <-> 2 site interaction
  testSuiteAoSForceCalculation(mol0, mol1, PPL, cutoff);

  // tests: 1 site <-> 2 site interaction, where sites are aligned such that all 3 sites are along the same line
  testSuiteAoSForceCalculation(mol0, mol2, PPL, cutoff);

  // tests: 2 site <-> 3 site interaction
  testSuiteAoSForceCalculation(mol1, mol3, PPL, cutoff);

  // tests: 3 site <-> 3 site interaction, where one has a nontrivial (needs rotating) quaternion
  testSuiteAoSForceCalculation(mol3, mol4, PPL, cutoff);

  // tests: 2 site <-> 2 site, where molecules are beyond cutoff
  testSuiteAoSForceCalculation(mol1, mol5, PPL, cutoff);

  // tests: 1 site <-> 2 site, where one site is beyond cutoff, the other within; and CoM beyond cutoff
  testSuiteAoSForceCalculation(mol0, mol6, PPL, cutoff);

  // tests: 1 site <-> 2 site, where one site is beyond cutoff, the other within; and CoM within cutoff
  testSuiteAoSForceCalculation(mol0, mol7, PPL, cutoff);

  // tests: 3 site <-> 3 site, with some different site types
  testSuiteAoSForceCalculation(mol4, mol8, PPL, cutoff);
}

/**
 * Tests that the AoS functor bypasses molecules that are dummies.
 */
TEST_F(LJMultisiteFunctorTest , AoSDummyTest) {
  using mdLib::MultisiteMoleculeLJ;

  const double cutoff = 2.5;

  ParticlePropertiesLibrary PPL(cutoff);
  PPL.addSiteType(0, 1., 1., 1.);
  PPL.addMolType(0, {0}, {{0., 0., 0.}}, {1., 1., 1.});
  PPL.calculateMixingCoefficients();

  MultisiteMoleculeLJ mol0;
  mol0.setR({0., 0., 0.});
  mol0.setQuaternion({0., 0., 0., 1.});
  mol0.setF({0., 0., 0.});
  mol0.setTorque({0., 0., 0.});
  mol0.setTypeId(0);
  mol0.setOwnershipState(autopas::OwnershipState::owned);

  MultisiteMoleculeLJ mol1;
  mol1.setR({-1., 0., 0.});
  mol1.setQuaternion({0., 0., 0., 1.});
  mol1.setF({0., 0., 0.});
  mol1.setTorque({0., 0., 0.});
  mol1.setTypeId(0);
  mol1.setOwnershipState(autopas::OwnershipState::dummy);

  MultisiteMoleculeLJ mol2;
  mol2.setR({1., 0., 0.});
  mol2.setQuaternion({0., 0., 0., 1.});
  mol2.setF({0., 0., 0.});
  mol2.setTorque({0., 0., 0.});
  mol2.setTypeId(0);
  mol2.setOwnershipState(autopas::OwnershipState::dummy);

  // Interact molecules together with newton3 on and off
  // create functor
  mdLib::LJMultisiteFunctor<mdLib::MultisiteMoleculeLJ, true, true, autopas::FunctorN3Modes::Both, true, true> functor(
      cutoff, PPL);

  // newton3 on
  functor.initTraversal();
  // Test (owned, dummy)
  functor.AoSFunctor(mol0, mol1, false);
  // Test (dummy, owned)
  functor.AoSFunctor(mol1, mol0, false);
  // Test (dummy, dummy)
  functor.AoSFunctor(mol1, mol2, false);
  functor.endTraversal(false);

  // newton3 on
  functor.initTraversal();
  // Test (owned, dummy)
  functor.AoSFunctor(mol0, mol1, true);
  // Test (dummy, owned)
  functor.AoSFunctor(mol1, mol0, true);
  // Test (dummy, dummy)
  functor.AoSFunctor(mol1, mol2, true);
  functor.endTraversal(true);

  // Test all forces and torques are zero
  // mol0
  EXPECT_DOUBLE_EQ(mol0.getF()[0], 0);
  EXPECT_DOUBLE_EQ(mol0.getF()[1], 0);
  EXPECT_DOUBLE_EQ(mol0.getF()[2], 0);
  // mol1
  EXPECT_DOUBLE_EQ(mol1.getF()[0], 0);
  EXPECT_DOUBLE_EQ(mol1.getF()[1], 0);
  EXPECT_DOUBLE_EQ(mol1.getF()[2], 0);
  // mol2
  EXPECT_DOUBLE_EQ(mol2.getF()[0], 0);
  EXPECT_DOUBLE_EQ(mol2.getF()[1], 0);
  EXPECT_DOUBLE_EQ(mol2.getF()[2], 0);
}

/**
 * Tests for correctness of AoS functor by comparing the force with FunctorLJ for a single-site molecule.
 */
TEST_F(LJMultisiteFunctorTest , singleSiteSanityCheck) {
  using mdLib::MultisiteMoleculeLJ;

  const double cutoff = 3.;

  ParticlePropertiesLibrary PPL(cutoff);
  PPL.addSiteType(0, 1., 1., 1.);
  PPL.addSiteType(1, 0.5, 0.5, 0.5);

  MultisiteMoleculeLJ mol0;
  mol0.setR({0., 0., 0.});
  mol0.setQuaternion({0., 0., 0., 1.});
  mol0.setF({0., 0., 0.});
  mol0.setTorque({0., 0., 0.});
  mol0.setTypeId(0);
  PPL.addMolType(0, {0}, {{0., 0., 0.}}, {1., 1., 1.});

  MultisiteMoleculeLJ mol1;
  mol1.setR({0.5, 0., 0.});
  mol1.setQuaternion({0., 0., 0., 1.});
  mol1.setF({0., 0., 0.});
  mol1.setTorque({0., 0., 0.});
  mol1.setTypeId(1);
  PPL.addMolType(1, {1}, {{0., 0., 0.}}, {1., 1., 1.});

  PPL.calculateMixingCoefficients();

  // N3L optimization disabled, global calculation disabled.
  singleSiteSanityCheck<false, false, false>(mol0, mol1, PPL, cutoff);

  // N3L optimization enabled, global calculation disabled.
  singleSiteSanityCheck<true, false, false>(mol0, mol1, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift disabled.
  singleSiteSanityCheck<false, true, false>(mol0, mol1, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift disabled.
  singleSiteSanityCheck<true, true, false>(mol0, mol1, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift enabled.
  singleSiteSanityCheck<false, true, true>(mol0, mol1, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift enabled.
  singleSiteSanityCheck<true, true, true>(mol0, mol1, PPL, cutoff);
}

/**
 * Tests SoAFunctorSingle using AoS functor as a reference.
 */
TEST_F(LJMultisiteFunctorTest , MultisiteLJFunctorTest_AoSVsSoASingle) {
  using mdLib::MultisiteMoleculeLJ;

  const double cutoff = 3.;

  std::vector<mdLib::MultisiteMoleculeLJ> allOwnedMolecules;
  std::vector<mdLib::MultisiteMoleculeLJ> mixedOwnershipMolecules;
  ParticlePropertiesLibrary<double, size_t> PPL(cutoff);

  generatePPL(&PPL);
  generateMolecules(&allOwnedMolecules);
  generateMolecules(&mixedOwnershipMolecules, {0, 0, 0}, false);

  // tests with only owned molecules

  // N3L optimization disabled, global calculation disabled.
  testSoACellAgainstAoS<false, false, false>(allOwnedMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation disabled.
  testSoACellAgainstAoS<true, false, false>(allOwnedMolecules, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift disabled.
  testSoACellAgainstAoS<false, true, false>(allOwnedMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift disabled.
  testSoACellAgainstAoS<true, true, false>(allOwnedMolecules, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift enabled.
  testSoACellAgainstAoS<false, true, true>(allOwnedMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift enabled.
  testSoACellAgainstAoS<true, true, true>(allOwnedMolecules, PPL, cutoff);

  // tests with a mix of ownership states

  // N3L optimization disabled, global calculation disabled.
  testSoACellAgainstAoS<false, false, false>(mixedOwnershipMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation disabled.
  testSoACellAgainstAoS<true, false, false>(mixedOwnershipMolecules, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift disabled.
  testSoACellAgainstAoS<false, true, false>(mixedOwnershipMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift disabled.
  testSoACellAgainstAoS<true, true, false>(mixedOwnershipMolecules, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift enabled.
  testSoACellAgainstAoS<false, true, true>(mixedOwnershipMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift enabled.
  testSoACellAgainstAoS<true, true, true>(mixedOwnershipMolecules, PPL, cutoff);
}

/**
 * Tests SoAFunctorPair using AoS functor as a reference.
 */
TEST_F(LJMultisiteFunctorTest , MultisiteLJFunctorTest_AoSVsSoAPair) {
  using mdLib::MultisiteMoleculeLJ;

  const double cutoff = 5.;

  std::vector<mdLib::MultisiteMoleculeLJ> allOwnedMoleculesA;
  std::vector<mdLib::MultisiteMoleculeLJ> allOwnedMoleculesB;
  std::vector<mdLib::MultisiteMoleculeLJ> mixedOwnershipMoleculesA;
  std::vector<mdLib::MultisiteMoleculeLJ> mixedOwnershipMoleculesB;
  ParticlePropertiesLibrary<double, size_t> PPL(cutoff);

  generatePPL(&PPL);
  generateMolecules(&allOwnedMoleculesA, {0, 0, 0});
  generateMolecules(&allOwnedMoleculesB, {0, 0, 9});
  generateMolecules(&mixedOwnershipMoleculesA, {0, 0, 0}, false);
  generateMolecules(&mixedOwnershipMoleculesB, {0, 0, 9}, false);

  // tests with only owned molecules

  // N3L optimization disabled, global calculation disabled.
  testSoACellPairAgainstAoS<false, false, false>(allOwnedMoleculesA, allOwnedMoleculesB, PPL, cutoff);

  // N3L optimization enabled, global calculation disabled.
  testSoACellPairAgainstAoS<true, false, false>(allOwnedMoleculesA, allOwnedMoleculesB, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift disabled.
  testSoACellPairAgainstAoS<false, true, false>(allOwnedMoleculesA, allOwnedMoleculesB, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift disabled.
  testSoACellPairAgainstAoS<true, true, false>(allOwnedMoleculesA, allOwnedMoleculesB, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift enabled.
  testSoACellPairAgainstAoS<false, true, true>(allOwnedMoleculesA, allOwnedMoleculesB, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift enabled.
  testSoACellPairAgainstAoS<true, true, true>(allOwnedMoleculesA, allOwnedMoleculesB, PPL, cutoff);

  // tests with a mix of ownership states

  // N3L optimization disabled, global calculation disabled.
  testSoACellPairAgainstAoS<false, false, false>(mixedOwnershipMoleculesA, mixedOwnershipMoleculesB, PPL, cutoff);

  // N3L optimization enabled, global calculation disabled.
  testSoACellPairAgainstAoS<true, false, false>(mixedOwnershipMoleculesA, mixedOwnershipMoleculesB, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift disabled.
  testSoACellPairAgainstAoS<false, true, false>(mixedOwnershipMoleculesA, mixedOwnershipMoleculesB, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift disabled.
  testSoACellPairAgainstAoS<true, true, false>(mixedOwnershipMoleculesA, mixedOwnershipMoleculesB, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift enabled.
  testSoACellPairAgainstAoS<false, true, true>(mixedOwnershipMoleculesA, mixedOwnershipMoleculesB, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift enabled.
  testSoACellPairAgainstAoS<true, true, true>(mixedOwnershipMoleculesA, mixedOwnershipMoleculesB, PPL, cutoff);
}

/**
 * Tests SoAFunctorVerlet using AoS functor as a reference.
 */
TEST_F(LJMultisiteFunctorTest , MultisiteLJFunctorTest_AoSVsSoAVerlet) {
  using mdLib::MultisiteMoleculeLJ;

  const double cutoff = 3.1;

  std::vector<mdLib::MultisiteMoleculeLJ> allOwnedMolecules;
  std::vector<mdLib::MultisiteMoleculeLJ> mixedOwnershipMolecules;
  ParticlePropertiesLibrary<double, size_t> PPL(cutoff);

  generatePPL(&PPL);
  generateMolecules(&allOwnedMolecules);
  generateMolecules(&mixedOwnershipMolecules, {0, 0, 0}, false);

  // tests with only owned molecules

  // N3L optimization disabled, global calculation disabled.
  testSoAVerletAgainstAoS<false, false, false>(allOwnedMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation disabled.
  testSoAVerletAgainstAoS<true, false, false>(allOwnedMolecules, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift disabled.
  testSoAVerletAgainstAoS<false, true, false>(allOwnedMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift disabled.
  testSoAVerletAgainstAoS<true, true, false>(allOwnedMolecules, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift enabled.
  testSoAVerletAgainstAoS<false, true, true>(allOwnedMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift enabled.
  testSoAVerletAgainstAoS<true, true, true>(allOwnedMolecules, PPL, cutoff);

  // tests with a mix of ownership states

  // N3L optimization disabled, global calculation disabled.
  testSoAVerletAgainstAoS<false, false, false>(mixedOwnershipMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation disabled.
  testSoAVerletAgainstAoS<true, false, false>(mixedOwnershipMolecules, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift disabled.
  testSoAVerletAgainstAoS<false, true, false>(mixedOwnershipMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift disabled.
  testSoAVerletAgainstAoS<true, true, false>(mixedOwnershipMolecules, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift enabled.
  testSoAVerletAgainstAoS<false, true, true>(mixedOwnershipMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift enabled.
  testSoAVerletAgainstAoS<true, true, true>(mixedOwnershipMolecules, PPL, cutoff);
}