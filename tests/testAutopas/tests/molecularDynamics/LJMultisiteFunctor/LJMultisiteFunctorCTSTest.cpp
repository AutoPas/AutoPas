/**
 * @file LJMultisiteFunctorCTSTest.cpp
 * @author Q. Behrami
 * @date 06/05/2023
 */
#include "LJMultisiteFunctorCTSTest.h"

#include <gtest/gtest.h>

#define AOS_VS_SOA_ACCURACY 1e-8
#define PARTICLES_PER_DIM 8

#include <chrono>

void LJMultisiteFunctorCTSTest::generatePPL(ParticlePropertiesLibrary<double, size_t> *PPL) {
  PPL->addSiteType(0, 1, 1, 1);
  PPL->addSiteType(1, 0.5, 0.5, 0.7);
  PPL->addMolType(0, {0}, {{0, 0, 0}}, {1, 1, 1});
  PPL->addMolType(1, {1, 0}, {{-0.05, 0, 0}, {0.05, 0, 0}}, {1, 1, 1});
  PPL->addMolType(2, {1, 0, 1, 0}, {{-0.025, 0, -0.025}, {-0.025, 0, 0.025}, {0.025, 0, -0.025}, {0.025, 0, 0.025}},
                  {1, 1, 1});

  PPL->calculateMixingCoefficients();
}

void LJMultisiteFunctorCTSTest::generateMolecules(std::vector<autopas::MultisiteMoleculeLJ> *molecules,
                                                  std::array<double, 3> offset = {0, 0, 0}) {
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
      }
    }
  }
}

template <bool newton3, bool calculateGlobals, bool applyShift>
void LJMultisiteFunctorCTSTest::testSoACellAgainstAoS(std::vector<autopas::MultisiteMoleculeLJ> molecules,
                                                      ParticlePropertiesLibrary<double, size_t> PPL, double cutoff) {
  using autopas::MultisiteMoleculeLJ;

  autopas::LJMultisiteFunctorCTS<MultisiteMoleculeLJ, applyShift, true, autopas::FunctorN3Modes::Both, calculateGlobals,
                                 true>
      functor(cutoff, PPL);

  auto moleculesAoS = molecules;
  auto moleculesSoA = molecules;
  const auto numberMolecules = molecules.size();

  // init traversal for functor
  functor.initTraversal();

  // Apply AoSFunctor to molecules
  for (size_t i = 0; i < numberMolecules; ++i) {
    for (size_t j = i + 1; j < numberMolecules; ++j) {
      functor.AoSFunctor(moleculesAoS[i], moleculesAoS[j], newton3);
    }
  }

  // end traversal for functor and get globals
  functor.endTraversal(newton3);
  const auto potentialEnergyAoS = calculateGlobals ? functor.getPotentialEnergy() : 0;
  const auto virialAoS = calculateGlobals ? functor.getVirial() : 0;

  // generate SoA Cell
  autopas::FullParticleCell<MultisiteMoleculeLJ> cellSoA;
  for (auto &&mol : moleculesSoA) {
    cellSoA.addParticle(mol);
  }

  // init traversal for functor
  functor.initTraversal();

  functor.SoALoader(cellSoA, cellSoA._particleSoABuffer, 0);
  // apply functor
  using namespace std::chrono;
  auto start = high_resolution_clock::now();
  functor.SoAFunctorSingle(cellSoA._particleSoABuffer, newton3);
  auto end = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(end - start);
  std::cout << "SoA functor took " << duration.count() << " microseconds." << std::endl;

  // copy back to original particle array
  moleculesSoA.clear();

  functor.SoAExtractor(cellSoA, cellSoA._particleSoABuffer, 0);

  // end traversal for functor and get globals
  functor.endTraversal(newton3);
  const auto potentialEnergySoA = calculateGlobals ? functor.getPotentialEnergy() : 0;
  const auto virialSoA = calculateGlobals ? functor.getVirial() : 0;

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

//  for (size_t i = 0; i < numberMolecules; ++i) {
//    EXPECT_NEAR(moleculesAoS[i].getTorque()[0], cellSoA._particles[i].getTorque()[0], AOS_VS_SOA_ACCURACY)
//        << "Incorrect x-torque for molecule " << i << " with newton3 = " << newton3;
//    EXPECT_NEAR(moleculesAoS[i].getTorque()[1], cellSoA._particles[i].getTorque()[1], AOS_VS_SOA_ACCURACY)
//        << "Incorrect y-torque for molecule " << i << " with newton3 = " << newton3;
//    EXPECT_NEAR(moleculesAoS[i].getTorque()[2], cellSoA._particles[i].getTorque()[2], AOS_VS_SOA_ACCURACY)
//        << "Incorrect z-torque for molecule " << i << " with newton3 = " << newton3;
//  }
//
//  if constexpr (calculateGlobals) {
//    EXPECT_NEAR(potentialEnergyAoS, potentialEnergySoA, AOS_VS_SOA_ACCURACY)
//        << "Incorrect potential energy with newton3 = " << newton3 << " and applyShift = " << applyShift;
//    EXPECT_NEAR(virialAoS, virialSoA, AOS_VS_SOA_ACCURACY)
//        << "Incorrect potential energy with newton3 = " << newton3 << " and applyShift = " << applyShift;
//  }
}

template <bool newton3, bool calculateGlobals, bool applyShift>
void LJMultisiteFunctorCTSTest::testSoACellPairAgainstAoS(std::vector<autopas::MultisiteMoleculeLJ> moleculesA,
                                                          std::vector<autopas::MultisiteMoleculeLJ> moleculesB,
                                                          ParticlePropertiesLibrary<double, size_t> PPL,
                                                          double cutoff) {
  using autopas::MultisiteMoleculeLJ;

  autopas::LJMultisiteFunctorCTS<MultisiteMoleculeLJ, applyShift, true, autopas::FunctorN3Modes::Both, calculateGlobals,
                                 true>
      functor(cutoff, PPL);

  auto moleculesAoSA = moleculesA;
  auto moleculesSoAA = moleculesA;
  const auto numberMoleculesA = moleculesA.size();

  auto moleculesAoSB = moleculesB;
  auto moleculesSoAB = moleculesB;
  const auto numberMoleculesB = moleculesB.size();

  // init traversal for functor
  functor.initTraversal();

  // Apply AoSFunctor to molecules
  for (size_t i = 0; i < numberMoleculesA; ++i) {
    for (size_t j = 0; j < numberMoleculesB; ++j) {
      functor.AoSFunctor(moleculesAoSA[i], moleculesAoSB[j], newton3);
    }
  }

  // end traversal for functor and get globals
  functor.endTraversal(newton3);
  const auto potentialEnergyAoS = calculateGlobals ? functor.getPotentialEnergy() : 0;
  const auto virialAoS = calculateGlobals ? functor.getVirial() : 0;

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
  functor.initTraversal();

  functor.SoALoader(cellSoAA, cellSoAA._particleSoABuffer, 0);
  functor.SoALoader(cellSoAB, cellSoAB._particleSoABuffer, 0);

  // apply functor
  using namespace std::chrono;
  auto start = high_resolution_clock::now();
  functor.SoAFunctorPair(cellSoAA._particleSoABuffer,cellSoAB._particleSoABuffer, newton3);
  auto end = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(end - start);
  std::cout << "SoA functor took " << duration.count() << " microseconds." << std::endl;

  // copy back to original particle array
  moleculesSoAA.clear();
  moleculesSoAB.clear();

  functor.SoAExtractor(cellSoAA, cellSoAA._particleSoABuffer, 0);
  functor.SoAExtractor(cellSoAB, cellSoAB._particleSoABuffer, 0);

  // end traversal for functor and get globals
  functor.endTraversal(newton3);
  const auto potentialEnergySoA = calculateGlobals ? functor.getPotentialEnergy() : 0;
  const auto virialSoA = calculateGlobals ? functor.getVirial() : 0;

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
    EXPECT_NEAR(potentialEnergyAoS, potentialEnergySoA, AOS_VS_SOA_ACCURACY)
        << "Incorrect potential energy with newton3 = " << newton3 << " and applyShift = " << applyShift;
    EXPECT_NEAR(virialAoS, virialSoA, AOS_VS_SOA_ACCURACY)
        << "Incorrect potential energy with newton3 = " << newton3 << " and applyShift = " << applyShift;
  }
}

template <bool newton3, bool calculateGlobals, bool applyShift>
void LJMultisiteFunctorCTSTest::testSoAVerletAgainstAoS(std::vector<autopas::MultisiteMoleculeLJ> molecules,
                                                        ParticlePropertiesLibrary<double, size_t> PPL, double cutoff) {
  using autopas::MultisiteMoleculeLJ;

  autopas::LJMultisiteFunctorCTS<MultisiteMoleculeLJ, applyShift, true, autopas::FunctorN3Modes::Both, calculateGlobals,
                                 true>
      functor(cutoff, PPL);

  auto moleculesAoS = molecules;
  auto moleculesSoA = molecules;
  const auto numberMolecules = molecules.size();

  // init traversal for functor
  functor.initTraversal();

  // Apply AoSFunctor to molecules
  for (size_t i = 0; i < numberMolecules; ++i) {
    for (size_t j = newton3 ? i + 1 : 0; j < numberMolecules; ++j) {
      if (i != j) {
        functor.AoSFunctor(moleculesAoS[i], moleculesAoS[j], newton3);
      }
    }
  }

  // end traversal for functor and get globals
  functor.endTraversal(newton3);
  const auto potentialEnergyAoS = calculateGlobals ? functor.getPotentialEnergy() : 0;
  const auto virialAoS = calculateGlobals ? functor.getVirial() : 0;

  // generate neighbor lists
  std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> neighborLists(numberMolecules);
  for (size_t i = 0; i < numberMolecules; ++i) {
    for (size_t j = i + 1; j < numberMolecules; ++j) {
      auto displacement = autopas::utils::ArrayMath::sub(moleculesSoA[i].getR(), moleculesSoA[j].getR());
      double distanceSquared = autopas::utils::ArrayMath::dot(displacement, displacement);
      if (distanceSquared <= cutoff * cutoff) {
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
  functor.initTraversal();

  functor.SoALoader(cellSoA, cellSoA._particleSoABuffer, 0);

  // apply functor
  using namespace std::chrono;
  auto start = high_resolution_clock::now();
  for (size_t i = 0; i < numberMolecules; ++i) {
    functor.SoAFunctorVerlet(cellSoA._particleSoABuffer, i, neighborLists[i], newton3);
  }
  auto end = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(end - start);
  std::cout << "SoA functor took " << duration.count() << " microseconds." << std::endl;

  // copy back to original particle array
  moleculesSoA.clear();

  functor.SoAExtractor(cellSoA, cellSoA._particleSoABuffer, 0);

  // end traversal for functor and get globals
  functor.endTraversal(newton3);
  const auto potentialEnergySoA = calculateGlobals ? functor.getPotentialEnergy() : 0;
  const auto virialSoA = calculateGlobals ? functor.getVirial() : 0;

  //   compare for consistency
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
      EXPECT_NEAR(potentialEnergyAoS, potentialEnergySoA, AOS_VS_SOA_ACCURACY)
          << "Incorrect potential energy with newton3 = " << newton3 << " and applyShift = " << applyShift;
      EXPECT_NEAR(virialAoS, virialSoA, AOS_VS_SOA_ACCURACY)
          << "Incorrect potential energy with newton3 = " << newton3 << " and applyShift = " << applyShift;
    }
}

/*
   Test functions
*/

/*
 * @note No newton3 disabled as SoACell always uses newton3 optimisation
 */
TEST_F(LJMultisiteFunctorCTSTest, AoSVsSoACell) {
  using autopas::MultisiteMoleculeLJ;

  const double cutoff = 1.5;

  std::vector<autopas::MultisiteMoleculeLJ> molecules;
  ParticlePropertiesLibrary<double, size_t> PPL(cutoff);

  generatePPL(&PPL);
  generateMolecules(&molecules);

  // N3L optimization enabled, global calculation disabled.
  testSoACellAgainstAoS<true, false, false>(molecules, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift disabled.
  testSoACellAgainstAoS<true, true, false>(molecules, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift enabled.
  testSoACellAgainstAoS<true, true, true>(molecules, PPL, cutoff);
}

TEST_F(LJMultisiteFunctorCTSTest, AoSVsSoACellPair) {
  using autopas::MultisiteMoleculeLJ;

  const double cutoff = 1.5;

  std::vector<autopas::MultisiteMoleculeLJ> moleculesA;
  std::vector<autopas::MultisiteMoleculeLJ> moleculesB;
  ParticlePropertiesLibrary<double, size_t> PPL(cutoff);

  generatePPL(&PPL);
  generateMolecules(&moleculesA, {0, 0, 0});
  generateMolecules(&moleculesB, {0, 0, 9});

  // N3L optimization disabled, global calculation disabled.
  testSoACellPairAgainstAoS<false, false, false>(moleculesA, moleculesB, PPL, cutoff);

  // N3L optimization enabled, global calculation disabled.
  testSoACellPairAgainstAoS<true, false, false>(moleculesA, moleculesB, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift disabled.
  testSoACellPairAgainstAoS<false, true, false>(moleculesA, moleculesB, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift disabled.
  testSoACellPairAgainstAoS<true, true, false>(moleculesA, moleculesB, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift enabled.
  testSoACellPairAgainstAoS<false, true, true>(moleculesA, moleculesB, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift enabled.
  testSoACellPairAgainstAoS<true, true, true>(moleculesA, moleculesB, PPL, cutoff);
}

TEST_F(LJMultisiteFunctorCTSTest, AoSVsSoAVerlet) {
  using autopas::MultisiteMoleculeLJ;

  const double cutoff = 1.5;

  std::vector<autopas::MultisiteMoleculeLJ> molecules;
  ParticlePropertiesLibrary<double, size_t> PPL(cutoff);

  generatePPL(&PPL);
  generateMolecules(&molecules);

  // N3L optimization disabled, global calculation disabled.
  testSoAVerletAgainstAoS<false, false, false>(molecules, PPL, cutoff);

  // N3L optimization enabled, global calculation disabled.
  testSoAVerletAgainstAoS<true, false, false>(molecules, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift disabled.
  testSoAVerletAgainstAoS<false, true, false>(molecules, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift disabled.
  testSoAVerletAgainstAoS<true, true, false>(molecules, PPL, cutoff);

  // N3L optimization disabled, global calculation enabled, apply shift enabled.
  testSoAVerletAgainstAoS<false, true, true>(molecules, PPL, cutoff);

  // N3L optimization enabled, global calculation enabled, apply shift enabled.
  testSoAVerletAgainstAoS<true, true, true>(molecules, PPL, cutoff);
}