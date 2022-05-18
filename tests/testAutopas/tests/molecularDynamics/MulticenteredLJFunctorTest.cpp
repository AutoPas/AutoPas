/**
 * @file MulticenteredLJFunctorTest.cpp
 * @author S. Newcome
 * @date 16/05/2022
 */

#include <gtest/gtest.h>

#include "MulticenteredLJFunctorTest.h"

template<bool newton3>
void MulticenteredLJFunctorTest::testAoSForceCalculation(autopas::MulticenteredMoleculeLJ molA, autopas::MulticenteredMoleculeLJ molB, ParticlePropertiesLibrary<double, size_t> PPL, double cutoff) {
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::sub;
  using autopas::utils::ArrayMath::dot;
  using autopas::utils::ArrayMath::cross;
  using autopas::utils::quaternion::rotateVectorOfPositions;

  const auto molTypeA = molA.getTypeId();
  const auto molTypeB = molB.getTypeId();
  const auto numberOfSitesA = PPL.getNumSites(molTypeA);
  const auto numberOfSitesB = PPL.getNumSites(molTypeB);
  const auto siteTypesA = PPL.getSiteTypes(molTypeA);
  const auto siteTypesB = PPL.getSiteTypes(molTypeB);


  // determine expected forces + torques
  std::array<double,3> expectedForceA{};
  std::array<double,3> expectedTorqueA{};
  std::array<double,3> expectedForceB{};
  std::array<double,3> expectedTorqueB{};

  if constexpr (not newton3) {
    expectedForceB = {0.,0.,0.};
    expectedTorqueB = {0.,0.,0.};
  }

  // determine if within cutoff
  if (dot(sub(molA.getR(),molB.getR()),sub(molA.getR(),molB.getR())) < cutoff * cutoff) {
    // calculate exact site positions
    const auto rotatedSitePositionsA = rotateVectorOfPositions(molA.getQ(), PPL.getSitePositions(molTypeA));
    const auto rotatedSitePositionsB = rotateVectorOfPositions(molB.getQ(), PPL.getSitePositions(molTypeB));

    for (size_t siteA = 0; siteA < numberOfSitesA; ++siteA) {
      const auto exactSitePositionA = add(rotatedSitePositionsA[siteA],molA.getR());
      for (size_t siteB = 0; siteB < numberOfSitesB; ++siteB) {
        const auto exactSitePositionB = add(rotatedSitePositionsB[siteB],molB.getR());

        const auto displacement = sub(exactSitePositionA,exactSitePositionB);
        const auto distanceSquared = dot(displacement,displacement);

        const auto sigmaSquared = PPL.getMixingData(siteTypesA[siteA],siteTypesB[siteB]).sigmaSquare;
        const auto epsilon24 = PPL.getMixingData(siteTypesA[siteA],siteTypesB[siteB]).epsilon24;

        // todo shift6

        const auto invDistSquared = 1. / distanceSquared;
        const auto lj2 = sigmaSquared * invDistSquared;
        const auto lj6 = lj2 * lj2 * lj2;
        const auto lj12 = lj6 * lj6;
        const auto lj12m6 = lj12 - lj6;  // = LJ potential / (4x epsilon)
        const auto scalarMultiple = epsilon24 * (lj12 + lj12m6) * invDistSquared;
        const auto force = autopas::utils::ArrayMath::mulScalar(displacement, scalarMultiple);

        expectedForceA = add(expectedForceA,force);
        if constexpr (newton3) {
          expectedForceB = sub(expectedForceB,force);
        }

        const auto torqueOnA = cross(rotatedSitePositionsA[siteA], force);
        expectedTorqueA = add(expectedTorqueA,torqueOnA);
        if constexpr (newton3) {
          const auto torqueOnB = cross(rotatedSitePositionsB[siteB], force);
          expectedTorqueB = sub(expectedTorqueB,torqueOnB);
        }
      }
    }
  } else {
    expectedForceA = {0,0,0};
    expectedTorqueA = {0,0,0};
    expectedForceB = {0,0,0};
    expectedTorqueB = {0,0,0};
  }

  // calculate forces and torques using AoS functor

  // create functor
  // todo add options for applyShift, useMixing, calculateGlobals
  autopas::LJMulticenterFunctor<autopas::MulticenteredMoleculeLJ, false, true, autopas::FunctorN3Modes::Both, false, true> functor(cutoff, PPL);

  functor.AoSFunctor(molA,molB,newton3);

  for (size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(molA.getF()[i], expectedForceA[i], 1e-13) << "molA: Unexpected force[" << i << "] = "
                                                                  << molA.getF()[i] << " != "
                                                                  << expectedForceA[i] << " as expected.";
  }
  for (size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(molA.getTorque()[i], expectedTorqueA[i], 1e-13) << "molA: Unexpected force[" << i << "] = "
                                                                  << molA.getTorque()[i] << " != "
                                                                  << expectedTorqueA[i] << " as expected.";
  }
  for (size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(molB.getF()[i], expectedForceB[i], 1e-13) << "molB: Unexpected force[" << i << "] = "
                                                                  << molB.getF()[i] << " != "
                                                                  << expectedForceB[i] << " as expected.";
  }
  for (size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(molB.getTorque()[i], expectedTorqueB[i], 1e-13) << "molB: Unexpected force[" << i << "] = "
                                                                        << molB.getTorque()[i] << " != "
                                                                        << expectedTorqueB[i] << " as expected.";
  }
}

/**
 * Tests for the correctness of the AoS functor by applying to molecules designed to test all its functionality.
 */
TEST_F(MulticenteredLJFunctorTest, AoSTest) {
  using autopas::MulticenteredMoleculeLJ;

  ParticlePropertiesLibrary PPL(1.);
  PPL.addSiteType(0,1.,1.,1.);
  PPL.addSiteType(1,0.5,0.5,0.5);

  // Molecules to be used in the tests (explanation of choices documented when tests are run).
  // For ease of readability, each molecule has its own molType, even when duplicated.
  MulticenteredMoleculeLJ mol0;
  mol0.setR({0.,0.,0.});
  mol0.setQ({1.,1.,0.,0.});
  mol0.setF({0.,0.,0.});
  mol0.setTorque({0.,0.,0.});
  mol0.setTypeId(0);
  PPL.addMolType(0,{0},{{0.,0.,0.}});

  MulticenteredMoleculeLJ mol1;
  mol1.setR({0.1,0.,0.});
  mol1.setQ({1.,1.,0.,0.});
  mol1.setF({0.,0.,0.});
  mol1.setTorque({0.,0.,0.});
  mol1.setTypeId(1);
  PPL.addMolType(1,{0,0},{{0.,0.01,0.},{0.,-0.01,0.}});

  MulticenteredMoleculeLJ mol2;
  mol2.setR({0.,0.1,0.});
  mol2.setQ({1.,1.,0.,0.});
  mol2.setF({0.,0.,0.});
  mol2.setTorque({0.,0.,0.});
  mol2.setTypeId(2);
  PPL.addMolType(2,{0,0},{{0.,0.01,0.},{0.,-0.01,0.}});

  MulticenteredMoleculeLJ mol3;
  mol3.setR({0.,0.,0.});
  mol3.setQ({1.,1.,0.,0.});
  mol3.setF({0.,0.,0.});
  mol3.setTorque({0.,0.,0.});
  mol3.setTypeId(3);
  PPL.addMolType(3,{0,0,0},{{-0.05,-0.05,0.},{0.,0.1,0.},{0.05,-0.05,0.}});

  MulticenteredMoleculeLJ mol4;
  mol4.setR({0.,0.,0.1});
  mol4.setQ({0.5,0.25, sqrt(3)/2,0.});
  mol4.setF({0.,0.,0.});
  mol4.setTorque({0.,0.,0.});
  mol4.setTypeId(4);
  PPL.addMolType(4,{0,0,0},{{-0.05,-0.05,0.},{0.,0.1,0.},{0.05,-0.05,0.}});

  MulticenteredMoleculeLJ mol5;
  mol5.setR({2.,2.,2.});
  mol5.setQ({1.,1.,0.,0.});
  mol5.setF({0.,0.,0.});
  mol5.setTorque({0.,0.,0.});
  mol5.setTypeId(5);
  PPL.addMolType(5,{0,0},{{0.,0.01,0.},{0.,-0.01,0.}});

  MulticenteredMoleculeLJ mol6;
  mol6.setR({0.,1.05,0.});
  mol6.setQ({1.,1.,0.,0.});
  mol6.setF({0.,0.,0.});
  mol6.setTorque({0.,0.,0.});
  mol6.setTypeId(6);
  PPL.addMolType(6,{0,0},{{0.,0.1,0.},{0.,-0.1,0.}});

  MulticenteredMoleculeLJ mol7;
  mol7.setR({0.,0.95,0.});
  mol7.setQ({1.,1.,0.,0.});
  mol7.setF({0.,0.,0.});
  mol7.setTorque({0.,0.,0.});
  mol7.setTypeId(7);
  PPL.addMolType(7,{0,0},{{0.,0.1,0.},{0.,-0.1,0.}});

  MulticenteredMoleculeLJ mol8;
  mol8.setR({0.,0.,0.});
  mol8.setQ({1.,1.,0.,0.});
  mol8.setF({0.,0.,0.});
  mol8.setTorque({0.,0.,0.});
  mol8.setTypeId(8);
  PPL.addMolType(8,{1,1,0},{{-0.05,-0.05,0.},{0.,0.1,0.},{0.05,-0.05,0.}});

  PPL.calculateMixingCoefficients();


  // ----------  newton3 = true   ----------

  // tests: 1 site <-> 2 site interaction
  testAoSForceCalculation<true>(mol0, mol1,PPL,1.);

  // tests: 1 site <-> 2 site interaction, where sites are aligned such that all 3 sites are along the same line
  testAoSForceCalculation<true>(mol0, mol2,PPL,1.);

  // tests: 2 site <-> 3 site interaction
  testAoSForceCalculation<true>(mol1, mol3,PPL,1.);

  // tests: 3 site <-> 3 site interaction, where one has a nontrivial (needs rotating) quaternion
  testAoSForceCalculation<true>(mol3, mol4,PPL,1.);

  // tests: 2 site <-> 2 site, where molecules are beyond cutoff
  testAoSForceCalculation<true>(mol1, mol5,PPL,1.);

  // tests: 1 site <-> 2 site, where one site is beyond cutoff, the other within; and CoM beyond cutoff
  testAoSForceCalculation<true>(mol0, mol6,PPL,1.);

  // tests: 1 site <-> 2 site, where one site is beyond cutoff, the other within; and CoM within cutoff
  testAoSForceCalculation<true>(mol0, mol7,PPL,1.);

  // tests: 3 site <-> 3 site, with some different site types
  testAoSForceCalculation<true>(mol4, mol8,PPL,1.);


  // ----------  newton3 = false  ----------

  // tests: 1 site <-> 2 site interaction
  testAoSForceCalculation<false>(mol0, mol1,PPL,1.);
  testAoSForceCalculation<false>(mol1, mol0,PPL,1.);

  // tests: 1 site <-> 2 site interaction, where sites are aligned such that all 3 sites are along the same line
  testAoSForceCalculation<false>(mol0, mol2,PPL,1.);
  testAoSForceCalculation<false>(mol2, mol0,PPL,1.);

  // tests: 2 site <-> 3 site interaction
  testAoSForceCalculation<false>(mol1, mol3,PPL,1.);
  testAoSForceCalculation<false>(mol3, mol1,PPL,1.);

  // tests: 3 site <-> 3 site interaction, where one has a nontrivial (needs rotating) quaternion
  testAoSForceCalculation<false>(mol3, mol4,PPL,1.);
  testAoSForceCalculation<false>(mol4, mol3,PPL,1.);

  // tests: 2 site <-> 2 site, where molecules are beyond cutoff
  testAoSForceCalculation<false>(mol1, mol5,PPL,1.);
  testAoSForceCalculation<false>(mol5, mol1,PPL,1.);

  // tests: 1 site <-> 2 site, where one site is beyond cutoff, the other within; and CoM beyond cutoff
  testAoSForceCalculation<false>(mol0, mol6,PPL,1.);
  testAoSForceCalculation<false>(mol6, mol0,PPL,1.);

  // tests: 1 site <-> 2 site, where one site is beyond cutoff, the other within; and CoM within cutoff
  testAoSForceCalculation<false>(mol0, mol7,PPL,1.);
  testAoSForceCalculation<false>(mol7, mol0,PPL,1.);

  // tests: 3 site <-> 3 site, with some different site types
  testAoSForceCalculation<false>(mol4, mol8,PPL,1.);
  testAoSForceCalculation<false>(mol8, mol4,PPL,1.);

}