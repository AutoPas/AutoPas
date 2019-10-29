/**
 * @file GeneratorsTest.h
 * @author N. Fottner
 * @date 3/8/19
 */

#include "GeneratorsTest.h"
#include "autopas/utils/inBox.h"
// the following test only work if testParsing.yaml is well set
// all ParticleVelocities = {0.,0.,0.}

TEST_F(GeneratorsTest, GridFillwithBoxMin) {
  auto autoPas = autopas::AutoPas<autopas::MoleculeLJ<>, FMCell>(std::cout);
  std::array<double, 3> boxmin = {5., 5., 5.};
  std::array<double, 3> boxmax = {10., 10., 10.};
  autoPas.setBoxMax(boxmax);
  autoPas.setBoxMin(boxmin);
  Molecule dummy;

  autoPas.init();
  GridGenerator::fillWithParticles(autoPas, {5, 5, 5}, 0, 0, dummy, {1, 1, 1}, boxmin);
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    EXPECT_TRUE(autopas::utils::inBox(iter->getR(), boxmin, boxmax));
  }
}

TEST_F(GeneratorsTest, MultipleObjectGeneration) {
  auto autoPas = autopas::AutoPas<autopas::MoleculeLJ<>, FMCell>(std::cout);
  MDFlexConfig config;
  config.yamlFilename = std::string(YAMLDIRECTORY) + "multipleObjectsWithMultipleTypesTest.yaml";
  YamlParser::parseYamlFile(config);
  config.calcSimulationBox();
  autoPas.setBoxMax(config.boxMax);
  autoPas.setBoxMin(config.boxMin);
  autoPas.init();
  std::array<double, 3> velocity = {0., 0., 0.};
  // parses the multiple Objects input of "multipleObjectsWithMultipleTypesTest" and generates a VTK File from the Input
  // yaml input file with all 4 different Object types
  auto cubeGrid(config.cubeGridObjects);
  auto cubeGauss(config.cubeGaussObjects);
  auto cubeUniform(config.cubeUniformObjects);
  auto sphere(config.sphereObjects);
  size_t idcounter = 0;  // to avoid multiple particles with the same ids

  Generator::CubeGrid(autoPas, cubeGrid.at(0).getTypeId(), idcounter, cubeGrid.at(0).getBoxMin(),
                      cubeGrid.at(0).getParticlesPerDim(), cubeGrid.at(0).getParticleSpacing(),
                      cubeGrid.at(0).getVelocity());
  idcounter += cubeGrid.at(0).getParticlesTotal();

  Generator::CubeGauss(autoPas, cubeGauss.at(0).getTypeId(), idcounter, cubeGauss.at(0).getBoxMin(),
                       cubeGauss.at(0).getBoxMax(), cubeGauss.at(0).getParticlesTotal(),
                       cubeGauss.at(0).getDistributionMean(), cubeGauss.at(0).getDistributionStdDev(),
                       cubeGauss.at(0).getVelocity());
  idcounter += cubeGauss.at(0).getParticlesTotal();

  Generator::CubeRandom(autoPas, cubeUniform.at(0).getTypeId(), idcounter, cubeUniform.at(0).getBoxMin(),
                        cubeUniform.at(0).getBoxMax(), cubeUniform.at(0).getParticlesTotal(),
                        cubeUniform.at(0).getVelocity());
  idcounter += cubeUniform.at(0).getParticlesTotal();

  Generator::Sphere(autoPas, sphere.at(0).getCenter(), sphere.at(0).getRadius(), sphere.at(0).getParticleSpacing(),
                    idcounter, sphere.at(0).getTypeId(), sphere.at(0).getVelocity());
  idcounter += sphere.at(0).getParticlesTotal();

  //  EXPECT_EQ(config.particlesTotal(), autoPas.getNumberOfParticles());
  //  EXPECT_EQ(idcounter, config.particlesTotal());
  // counters to checks if all particles types are well initialized for different Objects:
  int gridCounter = 0;
  int gaussCounter = 0;
  int uniformCounter = 0;
  int sphereCounter = 0;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
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
      default: { throw std::runtime_error("something went wrong with the Types"); }
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
