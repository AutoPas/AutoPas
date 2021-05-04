/**
 * @file GeneratorsTest.cpp
 * @author N. Fottner
 * @date 3/8/19
 */

#include "GeneratorsTest.h"

#include "autopas/AutoPas.h"
#include "autopasTools/generators/GridGenerator.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "src/parsing/YamlParser.h"
#include "testingHelpers/commonTypedefs.h"

TEST_F(GeneratorsTest, GridFillwithBoxMin) {
  auto autoPas = autopas::AutoPas<Molecule>(std::cout);
  std::array<double, 3> boxmin = {5., 5., 5.};
  std::array<double, 3> boxmax = {10., 10., 10.};
  autoPas.setBoxMax(boxmax);
  autoPas.setBoxMin(boxmin);
  Molecule dummy;

  autoPas.init();
  autopasTools::generators::GridGenerator::fillWithParticles(autoPas, {5, 5, 5}, dummy, {1, 1, 1}, boxmin);
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    EXPECT_TRUE(autopas::utils::inBox(iter->getR(), boxmin, boxmax));
  }
}

/**
 * This test expects multipleObjectsWithMultipleTypesTest.yaml to be placed in md-flexible/tests/yamlTestFiles
 */
TEST_F(GeneratorsTest, MultipleObjectGeneration) {
  auto autoPas = autopas::AutoPas<ParticleType>(std::cout);
  MDFlexConfig config;
  config.yamlFilename.value = std::string(YAMLDIRECTORY) + "multipleObjectsWithMultipleTypesTest.yaml";
  MDFlexParser::YamlParser::parseYamlFile(config);
  config.calcSimulationBox();
  autoPas.setBoxMax(config.boxMax.value);
  autoPas.setBoxMin(config.boxMin.value);
  autoPas.init();
  std::array<ParticleType::ParticleSoAFloatPrecision, 3> velocity = {0., 0., 0.};
  // parses the multiple Objects input of "multipleObjectsWithMultipleTypesTest" and generates a VTK File from the Input
  // yaml input file with all 4 different Object types
  auto cubeGrid(config.cubeGridObjects);
  auto cubeGauss(config.cubeGaussObjects);
  auto cubeUniform(config.cubeUniformObjects);
  auto sphere(config.sphereObjects);

  EXPECT_THAT(cubeGrid, ::testing::SizeIs(1));
  EXPECT_THAT(cubeGauss, ::testing::SizeIs(1));
  EXPECT_THAT(cubeUniform, ::testing::SizeIs(1));
  EXPECT_THAT(sphere, ::testing::SizeIs(1));

  size_t idcounter = 0;  // to avoid multiple particles with the same ids

  cubeGrid[0].generate(autoPas);
  idcounter += cubeGrid.at(0).getParticlesTotal();
  EXPECT_EQ(autoPas.getNumberOfParticles(), idcounter) << "CubeGrid generator added a wrong number of particles!";

  cubeGauss[0].generate(autoPas);
  idcounter += cubeGauss.at(0).getParticlesTotal();
  EXPECT_EQ(autoPas.getNumberOfParticles(), idcounter) << "CubeGauss generator added a wrong number of particles!";

  cubeUniform[0].generate(autoPas);
  idcounter += cubeUniform.at(0).getParticlesTotal();
  EXPECT_EQ(autoPas.getNumberOfParticles(), idcounter) << "CubeRandom generator added a wrong number of particles!";

  sphere[0].generate(autoPas);
  idcounter += sphere.at(0).getParticlesTotal();
  ASSERT_EQ(autoPas.getNumberOfParticles(), idcounter) << "Sphere generator added a wrong number of particles!";

  // counters to checks if all particles types are well initialized for different Objects:
  int gridCounter = 0;
  int gaussCounter = 0;
  int uniformCounter = 0;
  int sphereCounter = 0;

  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    EXPECT_EQ(velocity, iter->getV());  // velocity set to {0.,0.,0.} in parsingFile
    switch (iter->getTypeId()) {
      case 0: {
        gridCounter++;
        break;
      }
      case 1: {
        gaussCounter++;
        break;
      }
      case 2: {
        uniformCounter++;
        break;
      }
      case 3: {
        sphereCounter++;
        break;
      }
      default: {
        throw std::runtime_error("something went wrong with the Types");
      }
    }
  }
  EXPECT_EQ(gridCounter, cubeGrid.at(0).getParticlesTotal());
  EXPECT_EQ(gaussCounter, cubeGauss.at(0).getParticlesTotal());
  EXPECT_EQ(uniformCounter, cubeUniform.at(0).getParticlesTotal());
  EXPECT_EQ(sphereCounter, sphere.at(0).getParticlesTotal());
  // check if during initialization, not 2 Particles were initialized with same id
  std::set<size_t> ids;
  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    ASSERT_EQ(ids.count(iter->getID()), 0) << "Two particles have the same ID: " << iter->toString();
    ids.insert(iter->getID());
  }
  EXPECT_EQ(ids.size(), autoPas.getNumberOfParticles());
}
