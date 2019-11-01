/**
 * @file ParticlePropertiesLibraryTest.cpp
 * @author N. Fottner
 * @date 02/08/19
 */
#include "ParticlePropertiesLibraryTest.h"

double ParticlePropertiesLibraryTest::mixingE(double e1, double e2) { return std::sqrt(e1 * e2); }
double ParticlePropertiesLibraryTest::mixingS(double s1, double s2) { return ((s1 + s2) / 2); }

TEST_F(ParticlePropertiesLibraryTest, ClassFunctions) {
  PPL.addType(0, epsilon, sigma, mass);
  PPL.addType(1, epsilon2, sigma2, mass);
  PrintableMolecule p1({0., 0., 0.}, {0., 0., 0.}, 0);
  PrintableMolecule p2({0., 0., 0.}, {0., 0., 0.}, 1);
  p2.setTypeId(1);
  ASSERT_EQ(mass, PPL.getMass(p1.getTypeId()));
  ASSERT_EQ(PPL.mixing24Epsilon(p1.getTypeId(), p2.getTypeId()), 24 * mixingE(epsilon, epsilon2));
  ASSERT_EQ(PPL.mixingSigmaSquare(p1.getTypeId(), p2.getTypeId()), mixingS(sigma, sigma2) * mixingS(sigma, sigma2));
}

TEST_F(ParticlePropertiesLibraryTest, ParticlePropertiesInitialization) {
  Simulation<Molecule, FMCell> simulation;
  // this test need to be adapted if the input file changes
  MDFlexConfig config;
  config.yamlFilename = std::string(YAMLDIRECTORY) + "multipleObjectsWithMultipleTypesTest.yaml";
  YamlParser::parseYamlFile(config);
  config.calcSimulationBox();
  simulation.initialize(config);
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
  simulation.initialize(config);
  simulation.initializeParticlePropertiesLibrary();
  // default values: epsilon=sigma=mass=1.0
  EXPECT_EQ(simulation.getPpl()->getMass(0), 1.0);
  EXPECT_EQ(simulation.getPpl()->get24Epsilon(0), 24.0);
  EXPECT_EQ(simulation.getPpl()->getSigmaSquare(0), 1.0);
}