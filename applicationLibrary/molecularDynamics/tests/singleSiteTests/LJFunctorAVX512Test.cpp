/**
* @file LJFunctorAVX512Test.cpp
* @author S. Newcome
* @date 23/01/2024
*/

#include "LJFunctorAVX512Test.h"

#include <gtest/gtest.h>

#define PARTICLES_PER_DIM 7 // This should not be a multiple of 2
#define AOS_VS_SOA_ACCURACY 1e-8

void LJFunctorAVX512Test::generatePPL(ParticlePropertiesLibrary<double, size_t> *PPL) {
  PPL->addSiteType(0, 1, 1, 1);
  PPL->addSiteType(1, 0.5, 0.5, 0.7);
  PPL->addSiteType(2, 0.2, 1.4, 1.1);
  PPL->calculateMixingCoefficients();
}

void LJFunctorAVX512Test::generateMolecules(std::vector<mdLib::MoleculeLJ> *molecules,
                                               std::array<double, 3> offset = {0, 0, 0}, bool allOwned = true) {
  molecules->resize(PARTICLES_PER_DIM * PARTICLES_PER_DIM * PARTICLES_PER_DIM);

  for (unsigned int i = 0; i < PARTICLES_PER_DIM; ++i) {
    for (unsigned int j = 0; j < PARTICLES_PER_DIM; ++j) {
      for (unsigned int k = 0; k < PARTICLES_PER_DIM; ++k) {
        const auto index = i * PARTICLES_PER_DIM * PARTICLES_PER_DIM + j * PARTICLES_PER_DIM + k;
        molecules->at(index).setID(index);
        molecules->at(index).setR({(double)i + offset[0], (double)j + offset[1], (double)k + offset[2]});
        molecules->at(index).setF({0, 0, 0});
        molecules->at(index).setV({0, 0, 0});
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

template <template<class, bool, bool, autopas::FunctorN3Modes, bool, bool> class functorType, bool newton3, bool calculateGlobals, bool applyShift>
void LJFunctorAVX512Test::testAoSForceCalculation(mdLib::MoleculeLJ molA, mdLib::MoleculeLJ molB,
                                                     ParticlePropertiesLibrary<double, size_t> &PPL, double cutoff) {
  using namespace autopas::utils::ArrayMath::literals;
  using autopas::utils::ArrayMath::dot;

  const auto molTypeA = molA.getTypeId();
  const auto molTypeB = molB.getTypeId();

  // determine expected forces (+ globals)
  std::array<double, 3> expectedForceA{};
  std::array<double, 3> expectedForceB{};

  if constexpr (not newton3) {
    expectedForceB = {0., 0., 0.};
  }

  double expectedPotentialEnergySum{0.};
  double expectedVirialSum{0.};

  // determine if within cutoff

  const auto posA = molA.getR();
  const auto posB = molB.getR();

  const auto displacement = posA - posB;
  const auto distanceSquared = dot(displacement, displacement);

  if (distanceSquared < cutoff * cutoff) {

    const auto sigmaSquared = PPL.getMixingData(molA.getTypeId(), molB.getID()).sigmaSquared;
    const auto epsilon24 = PPL.getMixingData(molA.getTypeId(), molB.getID()).epsilon24;

    const auto invDistSquared = 1. / distanceSquared;
    const auto lj2 = sigmaSquared * invDistSquared;
    const auto lj6 = lj2 * lj2 * lj2;
    const auto lj12 = lj6 * lj6;
    const auto lj12m6 = lj12 - lj6;  // = LJ potential / (4x epsilon)
    const auto scalarMultiple = epsilon24 * (lj12 + lj12m6) * invDistSquared;
    const auto force = autopas::utils::ArrayMath::mulScalar(displacement, scalarMultiple);

    expectedForceA = force;
    if constexpr (newton3) {
      expectedForceB = force * -1.;
    }

    if constexpr (calculateGlobals) {
      const auto shift6 = applyShift ? PPL.getMixingData(molA.getTypeId(), molB.getTypeId()).shift6 : 0;
      const auto shift = shift6 / 6.;
      const auto epsilon4 = epsilon24 / 6.;

      // only add half the potential energy if no newton3 is used (because we are missing half the interaction
      // between the two molecules).
      const auto potentialEnergy = newton3 ? epsilon4 * lj12m6 + shift : 0.5 * (epsilon4 * lj12m6 + shift);
      const auto virialDimensionwiseContributions =
          newton3 ? autopas::utils::ArrayMath::mul(displacement, force)
                  : autopas::utils::ArrayMath::mulScalar(autopas::utils::ArrayMath::mul(displacement, force), 0.5);

      expectedPotentialEnergySum = potentialEnergy;
      expectedVirialSum = virialDimensionwiseContributions[0] + virialDimensionwiseContributions[1] +
                           virialDimensionwiseContributions[2];
    }
  } else {
    expectedForceA = {0, 0, 0};
    expectedForceB = {0, 0, 0};
  }

  // calculate forces and torques using AoS functor

  // create functor
  functorType<mdLib::MoleculeLJ, applyShift, true, autopas::FunctorN3Modes::Both,
                            calculateGlobals, true> functor(cutoff, PPL);

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
    EXPECT_NEAR(molB.getF()[i], expectedForceB[i], 1e-13)
        << "molB: Unexpected force[" << i << "] = " << molB.getF()[i] << " != " << expectedForceB[i]
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

template <template<class, bool, bool, autopas::FunctorN3Modes, bool, bool> class functorType>
void LJFunctorAVX512Test::testSuiteAoSForceCalculation(mdLib::MoleculeLJ molA,
                                                          mdLib::MoleculeLJ molB,
                                                          ParticlePropertiesLibrary<double, size_t> &PPL,
                                                          double cutoff) {
  // N3L Disabled, No Calculating Globals
  testAoSForceCalculation<functorType, false, false, false>(molA, molB, PPL, cutoff);

  // N3L Enabled, No Calculating Globals
  testAoSForceCalculation<functorType, true, false, false>(molA, molB, PPL, cutoff);

  // N3L Disabled, Calculating Globals, no shift applied
  testAoSForceCalculation<functorType, false, true, false>(molA, molB, PPL, cutoff);

  // N3L Disabled, Calculating Globals, shift applied
  testAoSForceCalculation<functorType, false, true, true>(molA, molB, PPL, cutoff);

  // N3L Enabled, Calculating Globals, no shift applied
  testAoSForceCalculation<functorType, true, true, false>(molA, molB, PPL, cutoff);

  // N3L Enabled, Calculating Globals, shift applied
  testAoSForceCalculation<functorType, true, true, true>(molA, molB, PPL, cutoff);
}

template <template<class, bool, bool, autopas::FunctorN3Modes, bool, bool> class functorType, bool newton3, bool calculateGlobals, bool applyShift>
void LJFunctorAVX512Test::testSoACellAgainstAoS(std::vector<mdLib::MoleculeLJ> molecules,
                                                   ParticlePropertiesLibrary<double, size_t> &PPL, double cutoff) {
  using mdLib::MoleculeLJ;

  functorType<MoleculeLJ, applyShift, true, autopas::FunctorN3Modes::Both, calculateGlobals, true>
      functorAoS(cutoff, PPL);
  functorType<MoleculeLJ, applyShift, true, autopas::FunctorN3Modes::Both, calculateGlobals, true>
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
  autopas::FullParticleCell<MoleculeLJ> cellSoA;
  for (auto &&mol : moleculesSoA) {
    cellSoA.addParticle(mol);
  }

  // init traversal for functor
  functorSoA.initTraversal();

  functorSoA.SoALoader(cellSoA, cellSoA._particleSoABuffer, 0, /*skipSoAResize*/ false);

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

template <template<class, bool, bool, autopas::FunctorN3Modes, bool, bool> class functorType, bool newton3, bool calculateGlobals, bool applyShift>
void LJFunctorAVX512Test::testSoACellPairAgainstAoS(std::vector<mdLib::MoleculeLJ> moleculesA,
                                                       std::vector<mdLib::MoleculeLJ> moleculesB,
                                                       ParticlePropertiesLibrary<double, size_t> &PPL, double cutoff) {
  using mdLib::MoleculeLJ;

  functorType<MoleculeLJ, applyShift, true, autopas::FunctorN3Modes::Both, calculateGlobals, true>
      functorAoS(cutoff, PPL);
  functorType<MoleculeLJ, applyShift, true, autopas::FunctorN3Modes::Both, calculateGlobals, true>
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
  autopas::FullParticleCell<MoleculeLJ> cellSoAA;
  for (auto &&mol : moleculesSoAA) {
    cellSoAA.addParticle(mol);
  }
  autopas::FullParticleCell<MoleculeLJ> cellSoAB;
  for (auto &&mol : moleculesSoAB) {
    cellSoAB.addParticle(mol);
  }

  // init traversal for functor
  functorSoA.initTraversal();

  functorSoA.SoALoader(cellSoAA, cellSoAA._particleSoABuffer, 0, /*skipSoAResize*/ false);
  functorSoA.SoALoader(cellSoAB, cellSoAB._particleSoABuffer, 0, /*skipSoAResize*/ false);

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


  for (size_t i = 0; i < numberMoleculesB; ++i) {
    EXPECT_NEAR(moleculesAoSB[i].getF()[0], cellSoAB._particles[i].getF()[0], AOS_VS_SOA_ACCURACY)
        << "Incorrect x-force for molecule " << i << " with newton3 = " << newton3;
    EXPECT_NEAR(moleculesAoSB[i].getF()[1], cellSoAB._particles[i].getF()[1], AOS_VS_SOA_ACCURACY)
        << "Incorrect y-force for molecule " << i << " with newton3 = " << newton3;
    EXPECT_NEAR(moleculesAoSB[i].getF()[2], cellSoAB._particles[i].getF()[2], AOS_VS_SOA_ACCURACY)
        << "Incorrect z-force for molecule " << i << " with newton3 = " << newton3;
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

template <template<class, bool, bool, autopas::FunctorN3Modes, bool, bool> class functorType, bool newton3, bool calculateGlobals, bool applyShift>
void LJFunctorAVX512Test::testSoAVerletAgainstAoS(std::vector<mdLib::MoleculeLJ> molecules,
                                                     ParticlePropertiesLibrary<double, size_t> &PPL, double cutoff) {
  using mdLib::MoleculeLJ;

  functorType<MoleculeLJ, applyShift, true, autopas::FunctorN3Modes::Both, calculateGlobals, true>
      functorAoS(cutoff, PPL);
  functorType<MoleculeLJ, applyShift, true, autopas::FunctorN3Modes::Both, calculateGlobals, true>
      functorSoA(cutoff, PPL);

  auto moleculesAoS = molecules;
  auto moleculesSoA = molecules;
  const auto numberMolecules = molecules.size();

  // generate neighbor lists
  std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> neighborLists(numberMolecules);
  for (size_t i = 0; i < numberMolecules; ++i) {
    for (size_t j = i + 1; j < numberMolecules; ++j) {
      if (i == j) {
        continue;
      }
      auto displacement = autopas::utils::ArrayMath::sub(moleculesAoS[i].getR(), moleculesAoS[j].getR());
      double distanceSquared = autopas::utils::ArrayMath::dot(displacement, displacement);
      if (distanceSquared < cutoff * cutoff) {
        neighborLists[i].push_back(j);
        if constexpr (not newton3) {
          neighborLists[j].push_back(i);
        }
      }
    }
  }

  // init traversal for functor
  functorAoS.initTraversal();

  // Apply AoSFunctor to molecules
  for (size_t i = 0; i < numberMolecules; ++i) {
    const auto neighborsOfI = neighborLists[i];
    for (unsigned long neighborIndex : neighborsOfI) {
      functorAoS.AoSFunctor(moleculesAoS[i], moleculesAoS[neighborIndex], newton3);
    }
  }

  // end traversal for functor
  functorAoS.endTraversal(newton3);

  // generate SoA Cell
  autopas::FullParticleCell<MoleculeLJ> cellSoA;
  for (auto &&mol : moleculesSoA) {
    cellSoA.addParticle(mol);
  }

  // init traversal for functor
  functorSoA.initTraversal();

  functorSoA.SoALoader(cellSoA, cellSoA._particleSoABuffer, 0, /*skipSoAResize*/ false);

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
TEST_F(LJFunctorAVX512Test, AoSTest) {
  using mdLib::MoleculeLJ;

  const double cutoff = 2.5;

  ParticlePropertiesLibrary PPL(cutoff);

  // Molecules to be used in the tests (explanation of choices documented when tests are run).
  // For ease of readability, each molecule has its own molType, even when duplicated.
  MoleculeLJ mol0;
  mol0.setR({0., 0., 0.});
  mol0.setF({0., 0., 0.});
  mol0.setTypeId(0);
  PPL.addSiteType(0, 1., 1., 1.);
  mol0.setOwnershipState(autopas::OwnershipState::owned);

  MoleculeLJ mol1;
  mol1.setR({1., 1., 1.});
  mol1.setF({0., 0., 0.});
  mol1.setTypeId(1);
  PPL.addSiteType(1, 1., 1., 1.);
  mol1.setOwnershipState(autopas::OwnershipState::owned);

  MoleculeLJ mol2;
  mol2.setR({1., 0.5, 0.});
  mol2.setF({0., 0., 0.});
  mol2.setTypeId(2);
  PPL.addSiteType(2, 0.5, 0.5, 1.);
  mol2.setOwnershipState(autopas::OwnershipState::owned);

  MoleculeLJ mol3;
  mol3.setR({0., 0., 1.});
  mol3.setF({0., 0., 0.});
  mol3.setTypeId(3);
  PPL.addSiteType(3, 1., 1., 1.);
  mol3.setOwnershipState(autopas::OwnershipState::halo);

  MoleculeLJ mol4;
  mol4.setR({0., 0.8, 0.});
  mol4.setF({0., 0., 0.});
  mol4.setTypeId(4);
  PPL.addSiteType(4, 1., 1., 1.);
  mol4.setOwnershipState(autopas::OwnershipState::dummy);

  MoleculeLJ mol5;
  mol5.setR({2.5, 2.5, 2.5});
  mol5.setF({0., 0., 0.});
  mol5.setTypeId(5);
  PPL.addSiteType(5, 1., 1., 1.);
  mol5.setOwnershipState(autopas::OwnershipState::dummy);

  PPL.calculateMixingCoefficients();

  testSuiteAoSForceCalculation<mdLib::LJFunctorAVX512_Mask>(mol0, mol1, PPL, 2.5);


}

/**
 * Tests that the AoS functor bypasses molecules that are dummies. Tests AutoVec, AVX512_Mask, AVX512_GS
 */
TEST_F(LJFunctorAVX512Test, AoSDummyTest) {
  using mdLib::MoleculeLJ;

  const double cutoff = 2.5;

  ParticlePropertiesLibrary PPL(cutoff);
  PPL.addSiteType(0, 1., 1., 1.);
  PPL.calculateMixingCoefficients();

  MoleculeLJ mol0;
  mol0.setR({0., 0., 0.});
  mol0.setF({0., 0., 0.});
  mol0.setTypeId(0);
  mol0.setOwnershipState(autopas::OwnershipState::owned);

  MoleculeLJ mol1;
  mol1.setR({-1., 0., 0.});
  mol1.setF({0., 0., 0.});
  mol1.setTypeId(0);
  mol1.setOwnershipState(autopas::OwnershipState::dummy);

 MoleculeLJ mol2;
  mol2.setR({1., 0., 0.});
  mol2.setF({0., 0., 0.});
  mol2.setTypeId(0);
  mol2.setOwnershipState(autopas::OwnershipState::dummy);

  // AutoVec

  // create functor
  mdLib::LJFunctorAVX512_Mask<mdLib::MoleculeLJ, true, true, autopas::FunctorN3Modes::Both, true, true> functor(cutoff, PPL);

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

  // Test all forces are zero
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
 * Tests SoAFunctorSingle using AoS functor as a reference.
 */
TEST_F(LJFunctorAVX512Test, MultisiteLJFunctorTest_AoSVsSoASingle) {
  using mdLib::MoleculeLJ;

  const double cutoff = 3.1;

  std::vector<mdLib::MoleculeLJ> allOwnedMolecules;
  std::vector<mdLib::MoleculeLJ> mixedOwnershipMolecules;
  ParticlePropertiesLibrary<double, size_t> PPL(cutoff);

  generatePPL(&PPL);
  generateMolecules(&allOwnedMolecules);
  generateMolecules(&mixedOwnershipMolecules, {0, 0, 0}, false);


  // tests with only owned molecules

  // N3L optimization disabled, global calculation disabled.
  testSoACellAgainstAoS<mdLib::LJFunctorAVX512_Mask, false, false, false>(allOwnedMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation disabled.
  testSoACellAgainstAoS<mdLib::LJFunctorAVX512_Mask, true, false, false>(allOwnedMolecules, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift disabled.
  testSoACellAgainstAoS<mdLib::LJFunctorAVX512_Mask, false, true, false>(allOwnedMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift disabled.
  testSoACellAgainstAoS<mdLib::LJFunctorAVX512_Mask, true, true, false>(allOwnedMolecules, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift enabled.
  testSoACellAgainstAoS<mdLib::LJFunctorAVX512_Mask, false, true, true>(allOwnedMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift enabled.
  testSoACellAgainstAoS<mdLib::LJFunctorAVX512_Mask, true, true, true>(allOwnedMolecules, PPL, cutoff);

  // tests with a mix of ownership states

  // N3L optimization disabled, global calculation disabled.
  testSoACellAgainstAoS<mdLib::LJFunctorAVX512_Mask, false, false, false>(mixedOwnershipMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation disabled.
  testSoACellAgainstAoS<mdLib::LJFunctorAVX512_Mask, true, false, false>(mixedOwnershipMolecules, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift disabled.
  testSoACellAgainstAoS<mdLib::LJFunctorAVX512_Mask, false, true, false>(mixedOwnershipMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift disabled.
  testSoACellAgainstAoS<mdLib::LJFunctorAVX512_Mask, true, true, false>(mixedOwnershipMolecules, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift enabled.
  testSoACellAgainstAoS<mdLib::LJFunctorAVX512_Mask, false, true, true>(mixedOwnershipMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift enabled.
  testSoACellAgainstAoS<mdLib::LJFunctorAVX512_Mask, true, true, true>(mixedOwnershipMolecules, PPL, cutoff);



}

/**
 * Tests SoAFunctorPair using AoS functor as a reference.
 */
TEST_F(LJFunctorAVX512Test, MultisiteLJFunctorTest_AoSVsSoAPair) {
  using mdLib::MoleculeLJ;

  const double cutoff = 5.1;

  std::vector<mdLib::MoleculeLJ> allOwnedMoleculesA;
  std::vector<mdLib::MoleculeLJ> allOwnedMoleculesB;
  std::vector<mdLib::MoleculeLJ> mixedOwnershipMoleculesA;
  std::vector<mdLib::MoleculeLJ> mixedOwnershipMoleculesB;
  ParticlePropertiesLibrary<double, size_t> PPL(cutoff);

  generatePPL(&PPL);
  generateMolecules(&allOwnedMoleculesA, {0, 0, 0});
  generateMolecules(&allOwnedMoleculesB, {0, 0, 9});
  generateMolecules(&mixedOwnershipMoleculesA, {0, 0, 0}, false);
  generateMolecules(&mixedOwnershipMoleculesB, {0, 0, 9}, false);


  // tests with only owned molecules

  // N3L optimization disabled, global calculation disabled.
  testSoACellPairAgainstAoS<mdLib::LJFunctorAVX512_Mask, false, false, false>(allOwnedMoleculesA, allOwnedMoleculesB, PPL, cutoff);

  // N3L optimization enabled, global calculation disabled.
  testSoACellPairAgainstAoS<mdLib::LJFunctorAVX512_Mask, true, false, false>(allOwnedMoleculesA, allOwnedMoleculesB, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift disabled.
  testSoACellPairAgainstAoS<mdLib::LJFunctorAVX512_Mask, false, true, false>(allOwnedMoleculesA, allOwnedMoleculesB, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift disabled.
  testSoACellPairAgainstAoS<mdLib::LJFunctorAVX512_Mask, true, true, false>(allOwnedMoleculesA, allOwnedMoleculesB, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift enabled.
  testSoACellPairAgainstAoS<mdLib::LJFunctorAVX512_Mask, false, true, true>(allOwnedMoleculesA, allOwnedMoleculesB, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift enabled.
  testSoACellPairAgainstAoS<mdLib::LJFunctorAVX512_Mask, true, true, true>(allOwnedMoleculesA, allOwnedMoleculesB, PPL, cutoff);

  // tests with a mix of ownership states

  // N3L optimization disabled, global calculation disabled.
  testSoACellPairAgainstAoS<mdLib::LJFunctorAVX512_Mask, false, false, false>(mixedOwnershipMoleculesA, mixedOwnershipMoleculesB, PPL, cutoff);

  // N3L optimization enabled, global calculation disabled.
  testSoACellPairAgainstAoS<mdLib::LJFunctorAVX512_Mask, true, false, false>(mixedOwnershipMoleculesA, mixedOwnershipMoleculesB, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift disabled.
  testSoACellPairAgainstAoS<mdLib::LJFunctorAVX512_Mask, false, true, false>(mixedOwnershipMoleculesA, mixedOwnershipMoleculesB, PPL, cutoff);

}

/**
 * Tests SoAFunctorVerlet using AoS functor as a reference.
 */
TEST_F(LJFunctorAVX512Test, MultisiteLJFunctorTest_AoSVsSoAVerlet) {
  using mdLib::MoleculeLJ;

  const double cutoff = 3.1;

  std::vector<mdLib::MoleculeLJ> allOwnedMolecules;
  std::vector<mdLib::MoleculeLJ> mixedOwnershipMolecules;
  ParticlePropertiesLibrary<double, size_t> PPL(cutoff);

  generatePPL(&PPL);
  generateMolecules(&allOwnedMolecules);
  generateMolecules(&mixedOwnershipMolecules, {0, 0, 0}, false);

  // AutoVec Tests

  // tests with only owned molecules

  // N3L optimization disabled, global calculation disabled.
  testSoAVerletAgainstAoS<mdLib::LJFunctorAVX512_Mask, false, false, false>(allOwnedMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation disabled.
  testSoAVerletAgainstAoS<mdLib::LJFunctorAVX512_Mask, true, false, false>(allOwnedMolecules, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift disabled.
  testSoAVerletAgainstAoS<mdLib::LJFunctorAVX512_Mask, false, true, false>(allOwnedMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift disabled.
  testSoAVerletAgainstAoS<mdLib::LJFunctorAVX512_Mask, true, true, false>(allOwnedMolecules, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift enabled.
  testSoAVerletAgainstAoS<mdLib::LJFunctorAVX512_Mask, false, true, true>(allOwnedMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift enabled.
  testSoAVerletAgainstAoS<mdLib::LJFunctorAVX512_Mask, true, true, true>(allOwnedMolecules, PPL, cutoff);

  // tests with a mix of ownership states

  // N3L optimization disabled, global calculation disabled.
  testSoAVerletAgainstAoS<mdLib::LJFunctorAVX512_Mask, false, false, false>(mixedOwnershipMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation disabled.
  testSoAVerletAgainstAoS<mdLib::LJFunctorAVX512_Mask, true, false, false>(mixedOwnershipMolecules, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift disabled.
  testSoAVerletAgainstAoS<mdLib::LJFunctorAVX512_Mask, false, true, false>(mixedOwnershipMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift disabled.
  testSoAVerletAgainstAoS<mdLib::LJFunctorAVX512_Mask, true, true, false>(mixedOwnershipMolecules, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift enabled.
  testSoAVerletAgainstAoS<mdLib::LJFunctorAVX512_Mask, false, true, true>(mixedOwnershipMolecules, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift enabled.
  testSoAVerletAgainstAoS<mdLib::LJFunctorAVX512_Mask, true, true, true>(mixedOwnershipMolecules, PPL, cutoff);



}