/**
 * @file ParticlePropertiesLibraryTest.cpp
 * @author N. Fottner
 * @date 02/08/19
 */
#include "ParticlePropertiesLibraryTest.h"

#include "autopas/AutoPas.h"
#include "autopas/molecularDynamics/LJFunctor.h"
#include "src/Simulation.h"
#include "src/parsing/YamlParser.h"
#include "testingHelpers/commonTypedefs.h"

double ParticlePropertiesLibraryTest::mixingE(double e1, double e2) { return std::sqrt(e1 * e2); }
double ParticlePropertiesLibraryTest::mixingS(double s1, double s2) { return ((s1 + s2) / 2); }

TEST_F(ParticlePropertiesLibraryTest, massTest) {
  ASSERT_EQ(PPL.getTypes().size(), masses.size());
  EXPECT_EQ(PPL.getMass(0), masses[0]);
  EXPECT_EQ(PPL.getMass(1), masses[1]);
}

TEST_F(ParticlePropertiesLibraryTest, epsilonTest) {
  ASSERT_EQ(PPL.getTypes().size(), epsilons.size());

  EXPECT_EQ(PPL.get24Epsilon(0), epsilons[0] * 24);
  EXPECT_EQ(PPL.get24Epsilon(1), epsilons[1] * 24);

  for (unsigned int i = 0; i < epsilons.size(); ++i) {
    for (unsigned int j = 0; j < epsilons.size(); ++j) {
      auto expectedVal = mixingE(epsilons[i], epsilons[j]) * 24;
      EXPECT_EQ(PPL.mixing24Epsilon(i, j), expectedVal) << "For i=" << i << " j=" << j;
    }
  }
}

TEST_F(ParticlePropertiesLibraryTest, sigmaTest) {
  ASSERT_EQ(PPL.getTypes().size(), sigmas.size());

  EXPECT_EQ(PPL.getSigmaSquare(0), sigmas[0] * sigmas[0]);
  EXPECT_EQ(PPL.getSigmaSquare(1), sigmas[1] * sigmas[1]);

  for (unsigned int i = 0; i < sigmas.size(); ++i) {
    for (unsigned int j = 0; j < sigmas.size(); ++j) {
      auto expectedVal = mixingS(sigmas[i], sigmas[j]);
      expectedVal *= expectedVal;
      EXPECT_EQ(PPL.mixingSigmaSquare(i, j), expectedVal) << "For i=" << i << " j=" << j;
    }
  }
}

TEST_F(ParticlePropertiesLibraryTest, shiftTest) {
  ASSERT_EQ(PPL.getTypes().size(), shifts.size());

  EXPECT_EQ(PPL.mixingShift6(0, 0), shifts[0]);
  EXPECT_EQ(PPL.mixingShift6(1, 1), shifts[1]);
}

/**
 * Idea: Two particles with distance of (almost) cutoff should produce (almost) zero shifted potential.
 */
TEST_F(ParticlePropertiesLibraryTest, mixedShiftTestUpot) {
  Molecule m1({0, 0, 0}, {0, 0, 0}, 0, 0);
  Molecule m2({cutoff - 1e-14, 0, 0}, {0, 0, 0}, 1, 1);

  autopas::LJFunctor<Molecule, autopas::FullParticleCell<Molecule>, /* shifting */ true, /*mixing*/ true,
                     autopas::FunctorN3Modes::Both,
                     /*globals*/ true>
      functor(cutoff, PPL);

  functor.initTraversal();
  functor.AoSFunctor(m1, m2, true);
  functor.endTraversal(true);
  EXPECT_NE(functor.getUpot(), 0);
  EXPECT_NEAR(functor.getUpot(), 0, 1e-10);
}

TEST_F(ParticlePropertiesLibraryTest, ParticlePropertiesInitialization) {
  Simulation<Molecule, FMCell> simulation;
  // this test need to be adapted if the input file changes
  MDFlexConfig config;
  autopas::AutoPas<Molecule, FMCell> autopas;
  config.yamlFilename = std::string(YAMLDIRECTORY) + "multipleObjectsWithMultipleTypesTest.yaml";
  YamlParser::parseYamlFile(config);
  config.calcSimulationBox();
  simulation.initialize(config, autopas);
  simulation.initializeParticlePropertiesLibrary();
  EXPECT_EQ(simulation.getPpl()->getMass(0), 1.0);
  EXPECT_EQ(simulation.getPpl()->get24Epsilon(0), 24.0);
  EXPECT_EQ(simulation.getPpl()->getSigmaSquare(0), 1.0);
  EXPECT_EQ(simulation.getPpl()->getMass(1), 2.);
  EXPECT_EQ(simulation.getPpl()->get24Epsilon(1), 48.0);
  EXPECT_EQ(simulation.getPpl()->getSigmaSquare(1), 4.0);
  EXPECT_EQ(simulation.getPpl()->getMass(2), 3.0);
  EXPECT_EQ(simulation.getPpl()->get24Epsilon(2), 72.0);
  EXPECT_EQ(simulation.getPpl()->getSigmaSquare(2), 9.0);
  EXPECT_EQ(simulation.getPpl()->getMass(3), 4.0);
  EXPECT_EQ(simulation.getPpl()->get24Epsilon(3), 96.0);
  EXPECT_EQ(simulation.getPpl()->getSigmaSquare(3), 16.0);
  ASSERT_ANY_THROW(simulation.getPpl()->get24Epsilon(10));
  ASSERT_ANY_THROW(simulation.getPpl()->getSigmaSquare(10));
  ASSERT_ANY_THROW(simulation.getPpl()->getMass(10));
}

TEST_F(ParticlePropertiesLibraryTest, ParticlePropertiesInitializationDefault) {
  // tests ParticleProperties initialization with only one Type
  Simulation<Molecule, FMCell> simulation;
  MDFlexConfig config;
  config.calcSimulationBox();
  autopas::AutoPas<Molecule, FMCell> autopas;
  simulation.initialize(config, autopas);
  simulation.initializeParticlePropertiesLibrary();
  // default values: epsilon=sigma=mass=1.0
  EXPECT_EQ(simulation.getPpl()->getMass(0), 1.0);
  EXPECT_EQ(simulation.getPpl()->get24Epsilon(0), 24.0);
  EXPECT_EQ(simulation.getPpl()->getSigmaSquare(0), 1.0);
}
