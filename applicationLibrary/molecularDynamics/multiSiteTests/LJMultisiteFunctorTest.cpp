/**
 * @file LJMultisiteFunctorTest.cpp
 * @author S. Newcome
 * @date 16/05/2022
 */

#include "LJMultisiteFunctorTest.h"

#include <gtest/gtest.h>

#define PARTICLES_PER_DIM 8
#define AOS_VS_SOA_ACCURACY 1e-8

void LJMultisiteFunctorTest::generatePPL(ParticlePropertiesLibrary<double, size_t> *PPL) {
  PPL->addSiteType(0, 1, 1, 1);
  PPL->addSiteType(1, 0.5, 0.5, 0.7);
  PPL->addMolType(0, {0}, {{0, 0, 0}}, {1, 1, 1});
  PPL->addMolType(1, {1, 0}, {{-0.05, 0, 0}, {0.05, 0, 0}}, {1, 1, 1});
  PPL->addMolType(2, {1, 0, 1, 0}, {{-0.025, 0, -0.025}, {-0.025, 0, 0.025}, {0.025, 0, -0.025}, {0.025, 0, 0.025}},
                  {1, 1, 1});

  PPL->calculateMixingCoefficients();
}

void LJMultisiteFunctorTest::generateMolecules(std::vector<mdLib::MultisiteMoleculeLJ> *molecules,
                                               std::array<double, 3> offset = {0, 0, 0}, bool allOwned = true) {
  molecules->resize(PARTICLES_PER_DIM * PARTICLES_PER_DIM * PARTICLES_PER_DIM);

  for (unsigned int i = 0; i < PARTICLES_PER_DIM; ++i) {
    for (unsigned int j = 0; j < PARTICLES_PER_DIM; ++j) {
      for (unsigned int k = 0; k < PARTICLES_PER_DIM; ++k) {
        const auto index = i * PARTICLES_PER_DIM * PARTICLES_PER_DIM + j * PARTICLES_PER_DIM + k;
        molecules->at(index).setID(index);
        molecules->at(index).setR({(double)i + offset[0], (double)j + offset[1], (double)k + offset[2]});
        // Generate quaternion deterministically but arbitrarily with fair variation
        const std::array<double, 4> qNonNormalized{1., (double)i + offset[0], (double)j + offset[1],
                                                   (double)k + offset[2]};
        molecules->at(index).setQ(autopas::utils::ArrayMath::normalize(qNonNormalized));
        molecules->at(index).setF({0, 0, 0});
        molecules->at(index).setTorque({0, 0, 0});
        molecules->at(index).setV({0, 0, 0});
        molecules->at(index).setAngularVel({0, 0, 0});
        molecules->at(index).setTypeId(index % 3);
        if (allOwned) {
          molecules->at(index).setOwnershipState(autopas::OwnershipState::owned);
        } else {
          if (index % 3 == 0) {
            molecules->at(index).setOwnershipState(autopas::OwnershipState::owned);
          } else if (index % 3 == 1) {
            molecules->at(index).setOwnershipState(autopas::OwnershipState::halo);
          } else {
            molecules->at(index).setOwnershipState(autopas::OwnershipState::dummy);
          }
        }
      }
    }
  }
}

template <bool newton3, bool calculateGlobals, bool applyShift>
void LJMultisiteFunctorTest::testAoSForceCalculation(mdLib::MultisiteMoleculeLJ molA, mdLib::MultisiteMoleculeLJ molB,
                                                     ParticlePropertiesLibrary<double, size_t> PPL, double cutoff) {
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::cross;
  using autopas::utils::ArrayMath::dot;
  using autopas::utils::ArrayMath::sub;
  using autopas::utils::quaternion::rotateVectorOfPositions;

  const auto molTypeA = molA.getTypeId();
  const auto molTypeB = molB.getTypeId();
  const auto numberOfSitesA = PPL.getNumSites(molTypeA);
  const auto numberOfSitesB = PPL.getNumSites(molTypeB);
  const auto siteTypesA = PPL.getSiteTypes(molTypeA);
  const auto siteTypesB = PPL.getSiteTypes(molTypeB);

  // determine expected forces + torques (+ globals)
  std::array<double, 3> expectedForceA{};
  std::array<double, 3> expectedTorqueA{};
  std::array<double, 3> expectedForceB{};
  std::array<double, 3> expectedTorqueB{};

  if constexpr (not newton3) {
    expectedForceB = {0., 0., 0.};
    expectedTorqueB = {0., 0., 0.};
  }

  double expectedPotentialEnergySum{0.};
  double expectedVirialSum{0.};

  // determine if within cutoff
  if (dot(sub(molA.getR(), molB.getR()), sub(molA.getR(), molB.getR())) < cutoff * cutoff) {
    // calculate exact site positions
    const auto rotatedSitePositionsA = rotateVectorOfPositions(molA.getQ(), PPL.getSitePositions(molTypeA));
    const auto rotatedSitePositionsB = rotateVectorOfPositions(molB.getQ(), PPL.getSitePositions(molTypeB));

    for (size_t siteA = 0; siteA < numberOfSitesA; ++siteA) {
      const auto exactSitePositionA = add(rotatedSitePositionsA[siteA], molA.getR());
      for (size_t siteB = 0; siteB < numberOfSitesB; ++siteB) {
        const auto exactSitePositionB = add(rotatedSitePositionsB[siteB], molB.getR());

        const auto displacement = sub(exactSitePositionA, exactSitePositionB);
        const auto distanceSquared = dot(displacement, displacement);

        const auto sigmaSquared = PPL.getMixingData(siteTypesA[siteA], siteTypesB[siteB]).sigmaSquared;
        const auto epsilon24 = PPL.getMixingData(siteTypesA[siteA], siteTypesB[siteB]).epsilon24;

        const auto invDistSquared = 1. / distanceSquared;
        const auto lj2 = sigmaSquared * invDistSquared;
        const auto lj6 = lj2 * lj2 * lj2;
        const auto lj12 = lj6 * lj6;
        const auto lj12m6 = lj12 - lj6;  // = LJ potential / (4x epsilon)
        const auto scalarMultiple = epsilon24 * (lj12 + lj12m6) * invDistSquared;
        const auto force = autopas::utils::ArrayMath::mulScalar(displacement, scalarMultiple);

        expectedForceA = add(expectedForceA, force);
        if constexpr (newton3) {
          expectedForceB = sub(expectedForceB, force);
        }

        const auto torqueOnA = cross(rotatedSitePositionsA[siteA], force);
        expectedTorqueA = add(expectedTorqueA, torqueOnA);
        if constexpr (newton3) {
          const auto torqueOnB = cross(rotatedSitePositionsB[siteB], force);
          expectedTorqueB = sub(expectedTorqueB, torqueOnB);
        }

        if constexpr (calculateGlobals) {
          const auto shift6 = applyShift ? PPL.getMixingData(siteTypesA[siteA], siteTypesB[siteB]).shift6 : 0;
          const auto shift = shift6 / 6.;
          const auto epsilon4 = epsilon24 / 6.;

          // only add half the potential energy if no newton3 is used (because we are missing half the interaction
          // between the two molecules).
          const auto potentialEnergy = newton3 ? epsilon4 * lj12m6 + shift : 0.5 * (epsilon4 * lj12m6 + shift);
          const auto virialDimensionwiseContributions =
              newton3 ? autopas::utils::ArrayMath::mul(displacement, force)
                      : autopas::utils::ArrayMath::mulScalar(autopas::utils::ArrayMath::mul(displacement, force), 0.5);

          expectedPotentialEnergySum += potentialEnergy;
          expectedVirialSum += virialDimensionwiseContributions[0] + virialDimensionwiseContributions[1] +
                               virialDimensionwiseContributions[2];
        }
      }
    }
  } else {
    expectedForceA = {0, 0, 0};
    expectedTorqueA = {0, 0, 0};
    expectedForceB = {0, 0, 0};
    expectedTorqueB = {0, 0, 0};
  }

  // calculate forces and torques using AoS functor

  // create functor
  mdLib::LJMultisiteFunctor<mdLib::MultisiteMoleculeLJ, applyShift, true, autopas::FunctorN3Modes::Both,
                            calculateGlobals, true>
      functor(cutoff, PPL);

  functor.initTraversal();
  functor.AoSFunctor(molA, molB, newton3);
  functor.endTraversal(newton3);

  for (size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(molA.getF()[i], expectedForceA[i], 1e-13)
        << "molA: Unexpected force[" << i << "] = " << molA.getF()[i] << " != " << expectedForceA[i]
        << " as expected with newton3 = " << newton3 << ", calculateGlobals = " << calculateGlobals
        << ", and applyShift = " << applyShift << ".";
  }
  for (size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(molA.getTorque()[i], expectedTorqueA[i], 1e-13)
        << "molA: Unexpected force[" << i << "] = " << molA.getTorque()[i] << " != " << expectedTorqueA[i]
        << " as expected with newton3 = " << newton3 << ", calculateGlobals = " << calculateGlobals
        << ", and applyShift = " << applyShift << ".";
  }
  for (size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(molB.getF()[i], expectedForceB[i], 1e-13)
        << "molB: Unexpected force[" << i << "] = " << molB.getF()[i] << " != " << expectedForceB[i]
        << " as expected with newton3 = " << newton3 << ", calculateGlobals = " << calculateGlobals
        << ", and applyShift = " << applyShift << ".";
  }
  for (size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(molB.getTorque()[i], expectedTorqueB[i], 1e-13)
        << "molB: Unexpected force[" << i << "] = " << molB.getTorque()[i] << " != " << expectedTorqueB[i]
        << " as expected with newton3 = " << newton3 << ", calculateGlobals = " << calculateGlobals
        << ", and applyShift = " << applyShift << ".";
  }
  if constexpr (calculateGlobals) {
    EXPECT_NEAR(functor.getPotentialEnergy(), expectedPotentialEnergySum, 1e-13)
        << "Unexpected potential energy = " << functor.getPotentialEnergy() << " != " << expectedPotentialEnergySum
        << " as expected with newton3 = " << newton3 << ", calculateGlobals = " << calculateGlobals
        << ", and applyShift = " << applyShift << ".";

    EXPECT_NEAR(functor.getVirial(), expectedVirialSum, 1e-13)
        << "Unexpected virial = " << functor.getVirial() << " != " << expectedVirialSum
        << " as expected with newton3 = " << newton3 << ", calculateGlobals = " << calculateGlobals
        << ", and applyShift = " << applyShift << ".";
  }
}

void LJMultisiteFunctorTest::testSuiteAoSForceCalculation(mdLib::MultisiteMoleculeLJ molA,
                                                          mdLib::MultisiteMoleculeLJ molB,
                                                          ParticlePropertiesLibrary<double, size_t> PPL,
                                                          double cutoff) {
  // N3L Disabled, No Calculating Globals
  testAoSForceCalculation<false, false, false>(molA, molB, PPL, cutoff);

  // N3L Enabled, No Calculating Globals
  testAoSForceCalculation<true, false, false>(molA, molB, PPL, cutoff);

  // N3L Disabled, Calculating Globals, no shift applied
  testAoSForceCalculation<false, true, false>(molA, molB, PPL, cutoff);

  // N3L Disabled, Calculating Globals, shift applied
  testAoSForceCalculation<false, true, true>(molA, molB, PPL, cutoff);

  // N3L Enabled, Calculating Globals, no shift applied
  testAoSForceCalculation<true, true, false>(molA, molB, PPL, cutoff);

  // N3L Enabled, Calculating Globals, shift applied
  testAoSForceCalculation<true, true, true>(molA, molB, PPL, cutoff);
}

template <bool newton3, bool calculateGlobals, bool applyShift>
void LJMultisiteFunctorTest::singleSiteSanityCheck(mdLib::MultisiteMoleculeLJ molA, mdLib::MultisiteMoleculeLJ molB,
                                                   ParticlePropertiesLibrary<double, size_t> PPL, double cutoff) {
  using mdLib::MoleculeLJ;

  // create functors
  mdLib::LJMultisiteFunctor<mdLib::MultisiteMoleculeLJ, applyShift, true, autopas::FunctorN3Modes::Both,
                            calculateGlobals, true>
      multiSiteFunctor(cutoff, PPL);
  mdLib::LJFunctor<mdLib::MoleculeLJ, applyShift, true, autopas::FunctorN3Modes::Both, calculateGlobals, true>
      singleSiteFunctor(cutoff, PPL);

  // create single site versions of the molecules
  mdLib::MoleculeLJ molASimple;
  molASimple.setTypeId(PPL.getSiteTypes(molA.getTypeId())[0]);
  molASimple.setR(molA.getR());
  molASimple.setV(molA.getV());
  molASimple.setF(molA.getF());
  molASimple.setOldF(molA.getOldF());
  molASimple.setID(molA.getID());
  molASimple.setOwnershipState(molA.getOwnershipState());

  mdLib::MoleculeLJ molBSimple;
  molBSimple.setTypeId(PPL.getSiteTypes(molB.getTypeId())[0]);
  molBSimple.setR(molB.getR());
  molBSimple.setV(molB.getV());
  molBSimple.setF(molB.getF());
  molBSimple.setOldF(molB.getOldF());
  molBSimple.setID(molB.getID());
  molBSimple.setOwnershipState(molB.getOwnershipState());

  // initialise traversals
  singleSiteFunctor.initTraversal();
  multiSiteFunctor.initTraversal();

  // apply functors
  singleSiteFunctor.AoSFunctor(molASimple, molBSimple, newton3);
  multiSiteFunctor.AoSFunctor(molA, molB, newton3);

  // end traversals
  singleSiteFunctor.endTraversal(newton3);
  multiSiteFunctor.endTraversal(newton3);

  for (size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(molA.getF()[i], molASimple.getF()[i], 1e-13);
  }
  for (size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(molB.getF()[i], molBSimple.getF()[i], 1e-13);
  }

  if constexpr (calculateGlobals) {
    EXPECT_NEAR(singleSiteFunctor.getPotentialEnergy(), multiSiteFunctor.getPotentialEnergy(), 1e-13);
    EXPECT_NEAR(singleSiteFunctor.getVirial(), multiSiteFunctor.getVirial(), 1e-13);
  }
}

template <bool newton3, bool calculateGlobals, bool applyShift>
void LJMultisiteFunctorTest::testSoACellAgainstAoS(std::vector<mdLib::MultisiteMoleculeLJ> molecules,
                                                   ParticlePropertiesLibrary<double, size_t> PPL, double cutoff) {
  using mdLib::MultisiteMoleculeLJ;

  mdLib::LJMultisiteFunctor<MultisiteMoleculeLJ, applyShift, true, autopas::FunctorN3Modes::Both, calculateGlobals,
                            true>
      functorAoS(cutoff, PPL);
  mdLib::LJMultisiteFunctor<MultisiteMoleculeLJ, applyShift, true, autopas::FunctorN3Modes::Both, calculateGlobals,
                            true>
      functorSoA(cutoff, PPL);

  auto moleculesAoS = molecules;
  auto moleculesSoA = molecules;
  const auto numberMolecules = molecules.size();

  // init traversal for functor
  functorAoS.initTraversal();

  // Apply AoSFunctor to molecules
  for (size_t i = 0; i < numberMolecules; ++i) {
    for (size_t j = newton3 ? i + 1 : 0; j < numberMolecules; ++j) {
      if (i == j) continue;
      functorAoS.AoSFunctor(moleculesAoS[i], moleculesAoS[j], newton3);
    }
  }

  // end traversal for functor
  functorAoS.endTraversal(newton3);

  // generate SoA Cell
  autopas::FullParticleCell<MultisiteMoleculeLJ> cellSoA;
  for (auto &&mol : moleculesSoA) {
    cellSoA.addParticle(mol);
  }

  // init traversal for functor
  functorSoA.initTraversal();

  functorSoA.SoALoader(cellSoA, cellSoA._particleSoABuffer, 0);

  // apply functor
  functorSoA.SoAFunctorSingle(cellSoA._particleSoABuffer, newton3);

  // copy back to original particle array
  moleculesSoA.clear();

  functorSoA.SoAExtractor(cellSoA, cellSoA._particleSoABuffer, 0);

  // end traversal for functor
  functorSoA.endTraversal(newton3);

  // compare for consistency
  EXPECT_EQ(moleculesAoS.size(), cellSoA.numParticles());

  for (size_t i = 0; i < numberMolecules; ++i) {
    EXPECT_NEAR(moleculesAoS[i].getF()[0], cellSoA._particles[i].getF()[0], AOS_VS_SOA_ACCURACY)
        << "Incorrect x-force for molecule " << i << " with newton3 = " << newton3;
    EXPECT_NEAR(moleculesAoS[i].getF()[1], cellSoA._particles[i].getF()[1], AOS_VS_SOA_ACCURACY)
        << "Incorrect y-force for molecule " << i << " with newton3 = " << newton3;
    EXPECT_NEAR(moleculesAoS[i].getF()[2], cellSoA._particles[i].getF()[2], AOS_VS_SOA_ACCURACY)
        << "Incorrect z-force for molecule " << i << " with newton3 = " << newton3;
  }

  for (size_t i = 0; i < numberMolecules; ++i) {
    EXPECT_NEAR(moleculesAoS[i].getTorque()[0], cellSoA._particles[i].getTorque()[0], AOS_VS_SOA_ACCURACY)
        << "Incorrect x-torque for molecule " << i << " with newton3 = " << newton3;
    EXPECT_NEAR(moleculesAoS[i].getTorque()[1], cellSoA._particles[i].getTorque()[1], AOS_VS_SOA_ACCURACY)
        << "Incorrect y-torque for molecule " << i << " with newton3 = " << newton3;
    EXPECT_NEAR(moleculesAoS[i].getTorque()[2], cellSoA._particles[i].getTorque()[2], AOS_VS_SOA_ACCURACY)
        << "Incorrect z-torque for molecule " << i << " with newton3 = " << newton3;
  }

  if constexpr (calculateGlobals) {
    const auto potentialEnergyAoS = functorAoS.getPotentialEnergy();
    const auto virialAoS = functorAoS.getVirial();
    const auto potentialEnergySoA = functorSoA.getPotentialEnergy();
    const auto virialSoA = functorSoA.getVirial();

    EXPECT_NEAR(potentialEnergyAoS, potentialEnergySoA, AOS_VS_SOA_ACCURACY)
        << "Incorrect potential energy with newton3 = " << newton3 << " and applyShift = " << applyShift;
    EXPECT_NEAR(virialAoS, virialSoA, AOS_VS_SOA_ACCURACY)
        << "Incorrect potential energy with newton3 = " << newton3 << " and applyShift = " << applyShift;
  }
}

template <bool newton3, bool calculateGlobals, bool applyShift>
void LJMultisiteFunctorTest::testSoACellPairAgainstAoS(std::vector<mdLib::MultisiteMoleculeLJ> moleculesA,
                                                       std::vector<mdLib::MultisiteMoleculeLJ> moleculesB,
                                                       ParticlePropertiesLibrary<double, size_t> PPL, double cutoff) {
  using mdLib::MultisiteMoleculeLJ;

  mdLib::LJMultisiteFunctor<MultisiteMoleculeLJ, applyShift, true, autopas::FunctorN3Modes::Both, calculateGlobals,
                            true>
      functorAoS(cutoff, PPL);
  mdLib::LJMultisiteFunctor<MultisiteMoleculeLJ, applyShift, true, autopas::FunctorN3Modes::Both, calculateGlobals,
                            true>
      functorSoA(cutoff, PPL);

  auto moleculesAoSA = moleculesA;
  auto moleculesSoAA = moleculesA;
  const auto numberMoleculesA = moleculesA.size();

  auto moleculesAoSB = moleculesB;
  auto moleculesSoAB = moleculesB;
  const auto numberMoleculesB = moleculesB.size();

  // init traversal for functor
  functorAoS.initTraversal();

  // Apply AoSFunctor to molecules
  for (size_t i = 0; i < numberMoleculesA; ++i) {
    for (size_t j = 0; j < numberMoleculesB; ++j) {
      functorAoS.AoSFunctor(moleculesAoSA[i], moleculesAoSB[j], newton3);
    }
  }

  // end traversal for functor
  functorAoS.endTraversal(newton3);

  // generate SoA Cells
  autopas::FullParticleCell<MultisiteMoleculeLJ> cellSoAA;
  for (auto &&mol : moleculesSoAA) {
    cellSoAA.addParticle(mol);
  }
  autopas::FullParticleCell<MultisiteMoleculeLJ> cellSoAB;
  for (auto &&mol : moleculesSoAB) {
    cellSoAB.addParticle(mol);
  }

  // init traversal for functor
  functorSoA.initTraversal();

  functorSoA.SoALoader(cellSoAA, cellSoAA._particleSoABuffer, 0);
  functorSoA.SoALoader(cellSoAB, cellSoAB._particleSoABuffer, 0);

  // apply functor
  functorSoA.SoAFunctorPair(cellSoAA._particleSoABuffer, cellSoAB._particleSoABuffer, newton3);

  // copy back to original particle array
  moleculesSoAA.clear();
  moleculesSoAB.clear();

  functorSoA.SoAExtractor(cellSoAA, cellSoAA._particleSoABuffer, 0);
  functorSoA.SoAExtractor(cellSoAB, cellSoAB._particleSoABuffer, 0);

  // end traversal for functor
  functorSoA.endTraversal(newton3);

  // compare for consistency
  EXPECT_EQ(moleculesAoSA.size(), cellSoAA.numParticles());
  EXPECT_EQ(moleculesAoSB.size(), cellSoAB.numParticles());

  for (size_t i = 0; i < numberMoleculesA; ++i) {
    EXPECT_NEAR(moleculesAoSA[i].getF()[0], cellSoAA._particles[i].getF()[0], AOS_VS_SOA_ACCURACY)
        << "Incorrect x-force for molecule " << i << " with newton3 = " << newton3;
    EXPECT_NEAR(moleculesAoSA[i].getF()[1], cellSoAA._particles[i].getF()[1], AOS_VS_SOA_ACCURACY)
        << "Incorrect y-force for molecule " << i << " with newton3 = " << newton3;
    EXPECT_NEAR(moleculesAoSA[i].getF()[2], cellSoAA._particles[i].getF()[2], AOS_VS_SOA_ACCURACY)
        << "Incorrect z-force for molecule " << i << " with newton3 = " << newton3;
  }
  for (size_t i = 0; i < numberMoleculesA; ++i) {
    EXPECT_NEAR(moleculesAoSA[i].getTorque()[0], cellSoAA._particles[i].getTorque()[0], AOS_VS_SOA_ACCURACY)
        << "Incorrect x-torque for molecule " << i << " with newton3 = " << newton3;
    EXPECT_NEAR(moleculesAoSA[i].getTorque()[1], cellSoAA._particles[i].getTorque()[1], AOS_VS_SOA_ACCURACY)
        << "Incorrect y-torque for molecule " << i << " with newton3 = " << newton3;
    EXPECT_NEAR(moleculesAoSA[i].getTorque()[2], cellSoAA._particles[i].getTorque()[2], AOS_VS_SOA_ACCURACY)
        << "Incorrect z-torque for molecule " << i << " with newton3 = " << newton3;
  }

  for (size_t i = 0; i < numberMoleculesB; ++i) {
    EXPECT_NEAR(moleculesAoSB[i].getF()[0], cellSoAB._particles[i].getF()[0], AOS_VS_SOA_ACCURACY)
        << "Incorrect x-force for molecule " << i << " with newton3 = " << newton3;
    EXPECT_NEAR(moleculesAoSB[i].getF()[1], cellSoAB._particles[i].getF()[1], AOS_VS_SOA_ACCURACY)
        << "Incorrect y-force for molecule " << i << " with newton3 = " << newton3;
    EXPECT_NEAR(moleculesAoSB[i].getF()[2], cellSoAB._particles[i].getF()[2], AOS_VS_SOA_ACCURACY)
        << "Incorrect z-force for molecule " << i << " with newton3 = " << newton3;
  }
  for (size_t i = 0; i < numberMoleculesB; ++i) {
    EXPECT_NEAR(moleculesAoSB[i].getTorque()[0], cellSoAB._particles[i].getTorque()[0], AOS_VS_SOA_ACCURACY)
        << "Incorrect x-torque for molecule " << i << " with newton3 = " << newton3;
    EXPECT_NEAR(moleculesAoSB[i].getTorque()[1], cellSoAB._particles[i].getTorque()[1], AOS_VS_SOA_ACCURACY)
        << "Incorrect y-torque for molecule " << i << " with newton3 = " << newton3;
    EXPECT_NEAR(moleculesAoSB[i].getTorque()[2], cellSoAB._particles[i].getTorque()[2], AOS_VS_SOA_ACCURACY)
        << "Incorrect z-torque for molecule " << i << " with newton3 = " << newton3;
  }

  if constexpr (calculateGlobals) {
    const auto potentialEnergyAoS = functorAoS.getPotentialEnergy();
    const auto virialAoS = functorAoS.getVirial();
    const auto potentialEnergySoA = functorSoA.getPotentialEnergy();
    const auto virialSoA = functorSoA.getVirial();

    EXPECT_NEAR(potentialEnergyAoS, potentialEnergySoA, AOS_VS_SOA_ACCURACY)
        << "Incorrect potential energy with newton3 = " << newton3 << " and applyShift = " << applyShift;
    EXPECT_NEAR(virialAoS, virialSoA, AOS_VS_SOA_ACCURACY)
        << "Incorrect potential energy with newton3 = " << newton3 << " and applyShift = " << applyShift;
  }
}

template <bool newton3, bool calculateGlobals, bool applyShift>
void LJMultisiteFunctorTest::testSoAVerletAgainstAoS(std::vector<mdLib::MultisiteMoleculeLJ> molecules,
                                                     ParticlePropertiesLibrary<double, size_t> PPL, double cutoff) {
  using mdLib::MultisiteMoleculeLJ;

  mdLib::LJMultisiteFunctor<MultisiteMoleculeLJ, applyShift, true, autopas::FunctorN3Modes::Both, calculateGlobals,
                            true>
      functorAoS(cutoff, PPL);
  mdLib::LJMultisiteFunctor<MultisiteMoleculeLJ, applyShift, true, autopas::FunctorN3Modes::Both, calculateGlobals,
                            true>
      functorSoA(cutoff, PPL);

  auto moleculesAoS = molecules;
  auto moleculesSoA = molecules;
  const auto numberMolecules = molecules.size();

  // init traversal for functor
  functorAoS.initTraversal();

  // Apply AoSFunctor to molecules
  for (size_t i = 0; i < numberMolecules; ++i) {
    for (size_t j = newton3 ? i + 1 : 0; j < numberMolecules; ++j) {
      if (i != j) {
        functorAoS.AoSFunctor(moleculesAoS[i], moleculesAoS[j], newton3);
      }
    }
  }

  // end traversal for functor
  functorAoS.endTraversal(newton3);

  // generate neighbor lists
  std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> neighborLists(numberMolecules);
  for (size_t i = 0; i < numberMolecules; ++i) {
    for (size_t j = i + 1; j < numberMolecules; ++j) {
      if (i == j) {
        continue;
      }
      auto displacement = autopas::utils::ArrayMath::sub(moleculesSoA[i].getR(), moleculesSoA[j].getR());
      double distanceSquared = autopas::utils::ArrayMath::dot(displacement, displacement);
      if (distanceSquared < cutoff * cutoff) {
        neighborLists[i].push_back(j);
        if constexpr (not newton3) {
          neighborLists[j].push_back(i);
        }
      }
    }
  }

  // generate SoA Cell
  autopas::FullParticleCell<MultisiteMoleculeLJ> cellSoA;
  for (auto &&mol : moleculesSoA) {
    cellSoA.addParticle(mol);
  }

  // init traversal for functor
  functorSoA.initTraversal();

  functorSoA.SoALoader(cellSoA, cellSoA._particleSoABuffer, 0);

  // apply functor
  for (size_t i = 0; i < numberMolecules; ++i) {
    functorSoA.SoAFunctorVerlet(cellSoA._particleSoABuffer, i, neighborLists[i], newton3);
  }

  // copy back to original particle array
  moleculesSoA.clear();

  functorSoA.SoAExtractor(cellSoA, cellSoA._particleSoABuffer, 0);

  // end traversal for functor
  functorSoA.endTraversal(newton3);

  // compare for consistency
  EXPECT_EQ(moleculesAoS.size(), cellSoA.numParticles());

  for (size_t i = 0; i < numberMolecules; ++i) {
    EXPECT_NEAR(moleculesAoS[i].getF()[0], cellSoA._particles[i].getF()[0], AOS_VS_SOA_ACCURACY)
        << "Incorrect x-force for molecule " << i << " with newton3 = " << newton3;
    EXPECT_NEAR(moleculesAoS[i].getF()[1], cellSoA._particles[i].getF()[1], AOS_VS_SOA_ACCURACY)
        << "Incorrect y-force for molecule " << i << " with newton3 = " << newton3;
    EXPECT_NEAR(moleculesAoS[i].getF()[2], cellSoA._particles[i].getF()[2], AOS_VS_SOA_ACCURACY)
        << "Incorrect z-force for molecule " << i << " with newton3 = " << newton3;
  }

  for (size_t i = 0; i < numberMolecules; ++i) {
    EXPECT_NEAR(moleculesAoS[i].getTorque()[0], cellSoA._particles[i].getTorque()[0], AOS_VS_SOA_ACCURACY)
        << "Incorrect x-torque for molecule " << i << " with newton3 = " << newton3;
    EXPECT_NEAR(moleculesAoS[i].getTorque()[1], cellSoA._particles[i].getTorque()[1], AOS_VS_SOA_ACCURACY)
        << "Incorrect y-torque for molecule " << i << " with newton3 = " << newton3;
    EXPECT_NEAR(moleculesAoS[i].getTorque()[2], cellSoA._particles[i].getTorque()[2], AOS_VS_SOA_ACCURACY)
        << "Incorrect z-torque for molecule " << i << " with newton3 = " << newton3;
  }

  if constexpr (calculateGlobals) {
    const auto potentialEnergyAoS = functorAoS.getPotentialEnergy();
    const auto virialAoS = functorAoS.getVirial();
    const auto potentialEnergySoA = functorSoA.getPotentialEnergy();
    const auto virialSoA = functorSoA.getVirial();

    EXPECT_NEAR(potentialEnergyAoS, potentialEnergySoA, AOS_VS_SOA_ACCURACY)
        << "Incorrect potential energy with newton3 = " << newton3 << " and applyShift = " << applyShift;
    EXPECT_NEAR(virialAoS, virialSoA, AOS_VS_SOA_ACCURACY)
        << "Incorrect potential energy with newton3 = " << newton3 << " and applyShift = " << applyShift;
  }
}

/**
 * Tests for the correctness of the AoS functor by applying to molecules designed to test all its functionality.
 */
TEST_F(LJMultisiteFunctorTest, AoSTest) {
  using mdLib::MultisiteMoleculeLJ;

  const double cutoff = 2.5;

  ParticlePropertiesLibrary PPL(cutoff);
  PPL.addSiteType(0, 1., 1., 1.);
  PPL.addSiteType(1, 0.5, 0.5, 0.5);

  // Molecules to be used in the tests (explanation of choices documented when tests are run).
  // For ease of readability, each molecule has its own molType, even when duplicated.
  MultisiteMoleculeLJ mol0;
  mol0.setR({0., 0., 0.});
  mol0.setQ({0., 0., 0., 1.});
  mol0.setF({0., 0., 0.});
  mol0.setTorque({0., 0., 0.});
  mol0.setTypeId(0);
  PPL.addMolType(0, {0}, {{0., 0., 0.}}, {1., 1., 1.});

  MultisiteMoleculeLJ mol1;
  mol1.setR({1., 0., 0.});
  mol1.setQ({0., 0., 0., 1.});
  mol1.setF({0., 0., 0.});
  mol1.setTorque({0., 0., 0.});
  mol1.setTypeId(1);
  PPL.addMolType(1, {0, 0}, {{0., 0.01, 0.}, {0., -0.01, 0.}}, {1., 1., 1.});

  MultisiteMoleculeLJ mol2;
  mol2.setR({0., 1., 0.});
  mol2.setQ({0., 0., 0., 1.});
  mol2.setF({0., 0., 0.});
  mol2.setTorque({0., 0., 0.});
  mol2.setTypeId(2);
  PPL.addMolType(2, {0, 0}, {{0., 0.01, 0.}, {0., -0.01, 0.}}, {1., 1., 1.});

  MultisiteMoleculeLJ mol3;
  mol3.setR({0., 0., 0.});
  mol3.setQ({0., 0., 0., 1.});
  mol3.setF({0., 0., 0.});
  mol3.setTorque({0., 0., 0.});
  mol3.setTypeId(3);
  PPL.addMolType(3, {0, 0, 0}, {{-0.05, -0.05, 0.}, {0., 0.1, 0.}, {0.05, -0.05, 0.}}, {1., 1., 1.});

  MultisiteMoleculeLJ mol4;
  mol4.setR({0., 0., 1.});
  mol4.setQ({0.7071067811865475, 0.7071067811865475, 0., 0.});
  mol4.setF({0., 0., 0.});
  mol4.setTorque({0., 0., 0.});
  mol4.setTypeId(4);
  PPL.addMolType(4, {0, 0, 0}, {{-0.05, -0.05, 0.}, {0., 0.1, 0.}, {0.05, -0.05, 0.}}, {1., 1., 1.});

  MultisiteMoleculeLJ mol5;
  mol5.setR({2.5, 2.5, 2.5});
  mol5.setQ({0., 0., 0., 1.});
  mol5.setF({0., 0., 0.});
  mol5.setTorque({0., 0., 0.});
  mol5.setTypeId(5);
  PPL.addMolType(5, {0, 0}, {{0., 0.01, 0.}, {0., -0.01, 0.}}, {1., 1., 1.});

  MultisiteMoleculeLJ mol6;
  mol6.setR({0., 2.55, 0.});
  mol6.setQ({0., 0., 0., 1.});
  mol6.setF({0., 0., 0.});
  mol6.setTorque({0., 0., 0.});
  mol6.setTypeId(6);
  PPL.addMolType(6, {0, 0}, {{0., 0.1, 0.}, {0., -0.1, 0.}}, {1., 1., 1.});

  MultisiteMoleculeLJ mol7;
  mol7.setR({0., 2.45, 0.});
  mol7.setQ({0., 0., 0., 1.});
  mol7.setF({0., 0., 0.});
  mol7.setTorque({0., 0., 0.});
  mol7.setTypeId(7);
  PPL.addMolType(7, {0, 0}, {{0., 0.1, 0.}, {0., -0.1, 0.}}, {1., 1., 1.});

  MultisiteMoleculeLJ mol8;
  mol8.setR({0., 0., 0.});
  mol8.setQ({0., 0., 0., 1.});
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
TEST_F(LJMultisiteFunctorTest, AoSDummyTest) {
  using mdLib::MultisiteMoleculeLJ;

  const double cutoff = 2.5;

  ParticlePropertiesLibrary PPL(cutoff);
  PPL.addSiteType(0, 1., 1., 1.);
  PPL.addMolType(0, {0}, {{0., 0., 0.}}, {1., 1., 1.});
  PPL.calculateMixingCoefficients();

  MultisiteMoleculeLJ mol0;
  mol0.setR({0., 0., 0.});
  mol0.setQ({0., 0., 0., 1.});
  mol0.setF({0., 0., 0.});
  mol0.setTorque({0., 0., 0.});
  mol0.setTypeId(0);
  mol0.setOwnershipState(autopas::OwnershipState::owned);

  MultisiteMoleculeLJ mol1;
  mol1.setR({-1., 0., 0.});
  mol1.setQ({0., 0., 0., 1.});
  mol1.setF({0., 0., 0.});
  mol1.setTorque({0., 0., 0.});
  mol1.setTypeId(0);
  mol1.setOwnershipState(autopas::OwnershipState::dummy);

  MultisiteMoleculeLJ mol2;
  mol2.setR({1., 0., 0.});
  mol2.setQ({0., 0., 0., 1.});
  mol2.setF({0., 0., 0.});
  mol2.setTorque({0., 0., 0.});
  mol2.setTypeId(0);
  mol2.setOwnershipState(autopas::OwnershipState::dummy);

  // Interact molecules together with newton3 on and off
  // create functor
  mdLib::LJMultisiteFunctor<mdLib::MultisiteMoleculeLJ, true, true, autopas::FunctorN3Modes::Both, true, true>
      functor(cutoff, PPL);

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
TEST_F(LJMultisiteFunctorTest, singleSiteSanityCheck) {
  using mdLib::MultisiteMoleculeLJ;

  const double cutoff = 3.;

  ParticlePropertiesLibrary PPL(cutoff);
  PPL.addSiteType(0, 1., 1., 1.);
  PPL.addSiteType(1, 0.5, 0.5, 0.5);

  MultisiteMoleculeLJ mol0;
  mol0.setR({0., 0., 0.});
  mol0.setQ({0., 0., 0., 1.});
  mol0.setF({0., 0., 0.});
  mol0.setTorque({0., 0., 0.});
  mol0.setTypeId(0);
  PPL.addMolType(0, {0}, {{0., 0., 0.}}, {1., 1., 1.});

  MultisiteMoleculeLJ mol1;
  mol1.setR({0.5, 0., 0.});
  mol1.setQ({0., 0., 0., 1.});
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
TEST_F(LJMultisiteFunctorTest, MultisiteLJFunctorTest_AoSVsSoASingle) {
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
TEST_F(LJMultisiteFunctorTest, MultisiteLJFunctorTest_AoSVsSoAPair) {
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
TEST_F(LJMultisiteFunctorTest, MultisiteLJFunctorTest_AoSVsSoAVerlet) {
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