//
// Created by johnny on 03.11.23.
//

#include "LJParentMultiSiteFunctorTest.h"

#include <gtest/gtest.h>

#define PARTICLES_PER_DIM 8
#define AOS_VS_SOA_ACCURACY 1e-8

template <class MoleculeType>
void LJParentMultiSiteFunctorTest<MoleculeType>::generatePPL(ParticlePropertiesLibrary<double, size_t> *PPL) {
  PPL->addSiteType(0, 1, 1, 1);
  PPL->addSiteType(1, 0.5, 0.5, 0.7);
  PPL->addMolType(0, {0}, {{0, 0, 0}}, {1, 1, 1});
  PPL->addMolType(1, {1, 0}, {{-0.05, 0, 0}, {0.05, 0, 0}}, {1, 1, 1});
  PPL->addMolType(2, {1, 0, 1, 0}, {{-0.025, 0, -0.025}, {-0.025, 0, 0.025}, {0.025, 0, -0.025}, {0.025, 0, 0.025}},
                  {1, 1, 1});

  PPL->calculateMixingCoefficients();
}

template <class MoleculeType>
void LJParentMultiSiteFunctorTest<MoleculeType>::generateMolecules(std::vector<MoleculeType> *molecules,
                                                             std::array<double, 3> offset, bool allOwned, const ParticlePropertiesLibrary<double, size_t>* ppl) {
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
        molecules->at(index).setQuaternion(autopas::utils::ArrayMath::normalize(qNonNormalized));
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

template <class MoleculeType>
template <bool newton3, bool calculateGlobals, bool applyShift>
void LJParentMultiSiteFunctorTest<MoleculeType>::testAoSForceCalculation(MoleculeType molA, MoleculeType molB,
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
    const auto rotatedSitePositionsA = rotateVectorOfPositions(molA.getQuaternion(), PPL.getSitePositions(molTypeA));
    const auto rotatedSitePositionsB = rotateVectorOfPositions(molB.getQuaternion(), PPL.getSitePositions(molTypeB));

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

template <class MoleculeType>
template <bool newton3, bool calculateGlobals, bool applyShift>
void LJParentMultiSiteFunctorTest<MoleculeType>::singleSiteSanityCheck(MoleculeType molA, MoleculeType molB,
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

template<class MoleculeType>
template <bool newton3, bool calculateGlobals, bool applyShift>
void LJParentMultiSiteFunctorTest<MoleculeType>::testSoACellAgainstAoS(std::vector<MoleculeType> molecules,
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
  EXPECT_EQ(moleculesAoS.size(), cellSoA.size());

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

template <class MoleculeType>
template <bool newton3, bool calculateGlobals, bool applyShift>
void LJParentMultiSiteFunctorTest<MoleculeType>::testSoACellPairAgainstAoS(std::vector<MoleculeType> moleculesA,
                                                                     std::vector<MoleculeType> moleculesB,
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
  EXPECT_EQ(moleculesAoSA.size(), cellSoAA.size());
  EXPECT_EQ(moleculesAoSB.size(), cellSoAB.size());

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

template <class MoleculeType>
template <bool newton3, bool calculateGlobals, bool applyShift>
void LJParentMultiSiteFunctorTest<MoleculeType>::testSoAVerletAgainstAoS(std::vector<MoleculeType> molecules,
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
  EXPECT_EQ(moleculesAoS.size(), cellSoA.size());

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

template <class MoleculeType>
void LJParentMultiSiteFunctorTest<MoleculeType>::testSuiteAoSForceCalculation(MoleculeType molA,
                                                                        MoleculeType molB,
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