/**
 * @file GeneratorsTest.cpp
 * @author N. Fottner
 * @date 3/8/19
 */

#include "GeneratorsTest.h"

#include "autopasTools/generators/GridGenerator.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "src/configuration/YamlParser.h"
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
  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename",
                                        std::string(YAMLDIRECTORY) + "multipleObjectsWithMultipleTypesTest.yaml"};

  char *argv[3] = {arguments[0].data(), arguments[1].data(), arguments[2].data()};

  MDFlexConfig configuration(3, argv);

  EXPECT_THAT(configuration.cubeGridObjects, ::testing::SizeIs(1));
  EXPECT_THAT(configuration.cubeGaussObjects, ::testing::SizeIs(1));
  EXPECT_THAT(configuration.cubeUniformObjects, ::testing::SizeIs(1));
  EXPECT_THAT(configuration.sphereObjects, ::testing::SizeIs(1));

  // counters to checks if all particles types are well initialized for different Objects:
  int gridCounter = 0;
  int gaussCounter = 0;
  int uniformCounter = 0;
  int sphereCounter = 0;

  std::array<double, 3> velocity = {0., 0., 0.};
  for (auto &particle : configuration.getParticles()) {
    EXPECT_EQ(velocity, particle.getV());  // velocity set to {0.,0.,0.} in parsingFile
    switch (particle.getTypeId()) {
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

  EXPECT_EQ(gridCounter, configuration.cubeGridObjects.at(0).getParticlesTotal());
  EXPECT_EQ(gaussCounter, configuration.cubeGaussObjects.at(0).getParticlesTotal());
  EXPECT_EQ(uniformCounter, configuration.cubeUniformObjects.at(0).getParticlesTotal());
  EXPECT_EQ(sphereCounter, configuration.sphereObjects.at(0).getParticlesTotal());
  // check if during initialization, not 2 Particles were initialized with same id
  std::set<size_t> ids;
  for (auto &particle : configuration.getParticles()) {
    int particleId = particle.getID();
    ASSERT_EQ(ids.count(particleId), 0) << "Two particles have the same ID " << particleId;
    ids.insert(particleId);
  }
  EXPECT_EQ(ids.size(), configuration.getParticles().size());
}
