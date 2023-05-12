/**
 * @file LJMultisiteFunctorCTSTest.cpp
 * @author Q. Behrami
 * @date 06/05/2023
 */
#include "LJMultisiteFunctorCTSTest.h"

#include <gtest/gtest.h>

#define AOS_VS_SOA_ACCURACY 1e-8
#define PARTICLES_PER_DIM 8

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
  // const auto potentialEnergyAoS = calculateGlobals ? functor.getPotentialEnergy() : 0;
  // const auto virialAoS = calculateGlobals ? functor.getVirial() : 0;

  // generate SoA Cell
  autopas::FullParticleCell<MultisiteMoleculeLJ> cellSoA;
  for (auto &&mol : moleculesSoA) {
    cellSoA.addParticle(mol);
  }

  // init traversal for functor
  functor.initTraversal();

  functor.SoALoader(cellSoA, cellSoA._particleSoABuffer, 0);
  // apply functor
  functor.SoAFunctorSingle(cellSoA._particleSoABuffer, newton3);

  // copy back to original particle array
  moleculesSoA.clear();

  functor.SoAExtractor(cellSoA, cellSoA._particleSoABuffer, 0);

  // end traversal for functor and get globals
  functor.endTraversal(newton3);
  // const auto potentialEnergySoA = calculateGlobals ? functor.getPotentialEnergy() : 0;
  // const auto virialSoA = calculateGlobals ? functor.getVirial() : 0;

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
//    std::cout << "Torque " << i << ": " << moleculesAoS[i].getTorque()[0] << " " << moleculesAoS[i].getTorque()[1]
//              << " " << moleculesAoS[i].getTorque()[2] << "| " << cellSoA._particles[i].getTorque()[0] << " "
//              << cellSoA._particles[i].getTorque()[1] << " " << cellSoA._particles[i].getTorque()[2] << std::endl;
//  }

   for (size_t i = 0; i < numberMolecules; ++i) {
     EXPECT_NEAR(moleculesAoS[i].getTorque()[0], cellSoA._particles[i].getTorque()[0], AOS_VS_SOA_ACCURACY)
         << "Incorrect x-torque for molecule " << i << " with newton3 = " << newton3;
     EXPECT_NEAR(moleculesAoS[i].getTorque()[1], cellSoA._particles[i].getTorque()[1], AOS_VS_SOA_ACCURACY)
         << "Incorrect y-torque for molecule " << i << " with newton3 = " << newton3;
     EXPECT_NEAR(moleculesAoS[i].getTorque()[2], cellSoA._particles[i].getTorque()[2], AOS_VS_SOA_ACCURACY)
         << "Incorrect z-torque for molecule " << i << " with newton3 = " << newton3;
   }

  // This is currently only working if shifting is applied
  // if constexpr (calculateGlobals && applyShift) {
  //   EXPECT_NEAR(potentialEnergyAoS, potentialEnergySoA, AOS_VS_SOA_ACCURACY)
  //       << "Incorrect potential energy with newton3 = " << newton3 << " and applyShift = " << applyShift;
  //   EXPECT_NEAR(virialAoS, virialSoA, AOS_VS_SOA_ACCURACY)
  //       << "Incorrect potential energy with newton3 = " << newton3 << " and applyShift = " << applyShift;
  // }
}

/*
   Test functions
*/

/*
 * @note No newton3 disabled as SoACell always uses newton3 optimisation
 */
TEST_F(LJMultisiteFunctorCTSTest, MulticenteredLJFunctorTest_AoSVsSoACell) {
  using autopas::MultisiteMoleculeLJ;

  const double cutoff = 1.5;

  std::vector<autopas::MultisiteMoleculeLJ> molecules;
  ParticlePropertiesLibrary<double, size_t> PPL(cutoff);

  generatePPL(&PPL);
  generateMolecules(&molecules);

  // N3L optimization enabled, global calculation disabled.
  testSoACellAgainstAoS<true, false, false>(molecules, PPL, cutoff);
  //
  // // N3L optimization enabled, global calculation enabled, apply shift disabled.
  // testSoACellAgainstAoS<true, true, false>(molecules, PPL, cutoff);
  //
  // // N3L optimization enabled, global calculation enabled, apply shift enabled.
  // testSoACellAgainstAoS<true, true, true>(molecules, PPL, cutoff);
}