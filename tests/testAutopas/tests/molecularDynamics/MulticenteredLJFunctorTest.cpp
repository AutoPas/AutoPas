/**
 * @file MulticenteredLJFunctorTest.cpp
 * @author S. Newcome
 * @date 16/05/2022
 */

#include <gtest/gtest.h>

#include "MulticenteredLJFunctorTest.h"

template<bool newton3>
void MulticenteredLJFunctorTest::testAoSForceCalculation(multisiteMolecule molA, multisiteMolecule molB, double cutoff) {
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::sub;
  using autopas::utils::ArrayMath::cross;
  using autopas::utils::ArrayMath::L2Norm;
  using autopas::utils::quaternion::rotateVectorOfPositions;

  const size_t numberOfSitesA = molA.sitePositions.size();
  const size_t numberOfSitesB = molB.sitePositions.size();
  // check size of site property vectors match
  EXPECT_EQ(numberOfSitesA,molA.siteEpsilons.size()) << "molA: Number of epsilons does not match number of sites";
  EXPECT_EQ(numberOfSitesA,molA.siteSigmas.size()) << "molA: Number of sigmas does not match number of sites";
  EXPECT_EQ(numberOfSitesB,molB.siteEpsilons.size()) << "molB: Number of epsilons does not match number of sites";
  EXPECT_EQ(numberOfSitesB,molB.siteSigmas.size()) << "molB: Number of sigmas does not match number of sites";

  // determine expected forces + torques
  std::array<double,3> expectedForceA;
  std::array<double,3> expectedTorqueA;
  std::array<double,3> expectedForceB;
  std::array<double,3> expectedTorqueB;

  if constexpr (not newton3) {
    expectedForceB = {0.,0.,0.};
    expectedTorqueB = {0.,0.,0.};
  }

  // determine if within cutoff
  if (L2Norm(sub(molA.CoMPosition,molB.CoMPosition)) < cutoff * cutoff) {
    // calculate exact site positions
    const auto rotatedSitePositionsA = rotateVectorOfPositions(molA.quaternion, molA.sitePositions);
    const auto rotatedSitePositionsB = rotateVectorOfPositions(molB.quaternion, molB.sitePositions);

    for (size_t siteA = 0; siteA < numberOfSitesA; ++siteA) {
      const auto exactSitePositionA = add(rotatedSitePositionsA[siteA],molA.CoMPosition);
      for (size_t siteB = 0; siteB < numberOfSitesB; ++siteB) {
        const auto exactSitePositionB = add(rotatedSitePositionsB[siteB],molB.CoMPosition);

        const auto displacement = sub(exactSitePositionA,exactSitePositionB);
        const auto distanceSquared = L2Norm(displacement);

        const auto sigma = (molA.siteSigmas[siteA] + molB.siteSigmas[siteB]) / 2.0;
        const auto sigmaSquared = sigma * sigma;
        const auto epsilon24 = 24 * sqrt(molA.siteEpsilons[siteA] * molB.siteEpsilons[siteB]);

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

  // create + fill PPL
  ParticlePropertiesLibrary PPL(cutoff);
  std::vector<size_t> molASiteIds;
  for (int i = 0; i < molA.sitePositions.size(); ++i) {
    PPL.addSiteType(i,molA.siteEpsilons[i],molA.siteSigmas[i],1);
    molASiteIds[i] = i;
  }
  PPL.addMolType(0, molASiteIds, molA.sitePositions);

  std::vector<size_t> molBSiteIds;
  for (int i = 0; i < molB.sitePositions.size(); ++i) {
    PPL.addSiteType(i+molA.sitePositions.size(),molB.siteEpsilons[i],molB.siteSigmas[i],1);
    molBSiteIds[i] = i+molA.sitePositions.size();
  }
  PPL.addMolType(1, molBSiteIds, molB.sitePositions);

  // create molecules
  autopas::MulticenteredMoleculeLJ molAParticle, molBParticle;

  molAParticle.setR(molA.CoMPosition);
  molAParticle.setQ(molA.quaternion);
  molAParticle.setF({0,0,0});
  molAParticle.setTorque({0,0,0});
  molAParticle.setTypeId(0);
  molAParticle.setID(0);

  molBParticle.setR(molB.CoMPosition);
  molBParticle.setQ(molB.quaternion);
  molBParticle.setF({0,0,0});
  molBParticle.setTorque({0,0,0});
  molBParticle.setTypeId(0);
  molBParticle.setID(1);

  // create autopass container, functor, and add particles
  // todo add options for applyShift, useMixing, calculateGlobals
  autopas::LJMulticenterFunctor<autopas::MulticenteredMoleculeLJ, false, true, autopas::FunctorN3Modes::Both, false, true> functor(cutoff, PPL);

  functor.AoSFunctor(molAParticle,molBParticle,newton3);

  for (size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(molAParticle.getF()[i], expectedForceA[i], 1e-13) << "molA: Unexpected force[" << i << "] = "
                                                                  << molAParticle.getF()[i] << " != "
                                                                  << expectedForceA[i] << " as expected.";
  }
  for (size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(molAParticle.getTorque()[i], expectedTorqueA[i], 1e-13) << "molA: Unexpected force[" << i << "] = "
                                                                  << molAParticle.getTorque()[i] << " != "
                                                                  << expectedTorqueA[i] << " as expected.";
  }
  for (size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(molBParticle.getF()[i], expectedForceB[i], 1e-13) << "molB: Unexpected force[" << i << "] = "
                                                                  << molBParticle.getF()[i] << " != "
                                                                  << expectedForceB[i] << " as expected.";
  }
  for (size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(molBParticle.getTorque()[i], expectedTorqueB[i], 1e-13) << "molB: Unexpected force[" << i << "] = "
                                                                        << molBParticle.getTorque()[i] << " != "
                                                                        << expectedTorqueB[i] << " as expected.";
  }
}

TEST_F(MulticenteredLJFunctorTest, AoSTest) {
  std::vector<multisiteMolecule> mols;

  mols[0].CoMPosition = {0.,0.,0.};
  mols[0].quaternion = {0.,1.,0.,0.};
  mols[0].sitePositions = {{0.,0.,0.}};
  mols[0].siteEpsilons = {1.};
  mols[0].siteSigmas = {1.};

  mols[1].CoMPosition = {0.1,0.,0.};
  mols[1].quaternion = {0.,1.,0.,0.};
  mols[1].sitePositions = {{0.,0.01,0.},{0.,-0.01,0.}};
  mols[1].siteEpsilons = {1.,1.};
  mols[1].siteSigmas = {1.,1.};

  mols[2].CoMPosition = {0.,0.1,0.};
  mols[2].quaternion = {0.,1.,0.,0.};
  mols[2].sitePositions = {{0.,0.01,0.},{0.,-0.01,0.}};
  mols[2].siteEpsilons = {1.,1.};
  mols[2].siteSigmas = {1.,1.};

  mols[3].CoMPosition = {0.,0.,0.};
  mols[3].quaternion = {0.,1.,0.,0.};
  mols[3].sitePositions = {{-0.05,-0.05,0.},{0.,0.1,0.},{0.05,-0.05,0.}};
  mols[3].siteEpsilons = {1.,1.,1.};
  mols[3].siteSigmas = {1.,1.,1.};

  mols[4].CoMPosition = {0.,0.,0.1};
  mols[4].quaternion = {0.5,0.25, sqrt(3)/2,0.};
  mols[4].sitePositions = {{-0.05,-0.05,0.},{0.,0.1,0.},{0.05,-0.05,0.}};
  mols[4].siteEpsilons = {1.,1.,1.};
  mols[4].siteSigmas = {1.,1.,1.};

  // tests: 1 site <-> 2 site interaction
  testAoSForceCalculation<true>(mols[0],mols[1],1.);

  // tests: 1 site <-> 2 site interaction, where sites are aligned such that all 3 sites are along the same line
  testAoSForceCalculation<true>(mols[0],mols[2],1.);

  // tests: 2 site <-> 3 site interaction
  testAoSForceCalculation<true>(mols[1],mols[3],1.);

  // tests: 3 site <-> 3 site interaction, where one has a nontrivial (needs rotating) quaternion
  testAoSForceCalculation<true>(mols[3],mols[4],1.);

  // tests: 2 site <-> 2 site, where molecules are beyond cutoff
  testAoSForceCalculation<true>(mols[1],mols[2],0.001);

  // tests: 1 site <-> 2 site, where one site is beyond cutoff, the other within; and CoM beyond cutoff
  testAoSForceCalculation<true>(mols[0],mols[2],0.999);

  // tests: 1 site <-> 2 site, where one site is beyond cutoff, the other within; and CoM within cutoff
  testAoSForceCalculation<true>(mols[0],mols[2],1.001);

}