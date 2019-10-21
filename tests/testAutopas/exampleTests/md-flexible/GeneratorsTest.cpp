//
// Created by nicola on 216.06.19.
//

#include "GeneratorsTest.h"
#include "autopas/utils/inBox.h"
// the following test only work if testParsing.yaml is well set
// all ParticleVelocities = {0.,0.,0.}

TEST_F(GeneratorsTest, Gauss) {
  auto autoPas = autopas::AutoPas<autopas::MoleculeLJ<>, FMCell>(std::cout);
  Molecule dummyParticle;
  autoPas.setBoxMax({20., 20., 20.});
  autoPas.setBoxMin({0., 0., 0.});
  std::array<double, 3> velocity = {0., 0., 0.};
  autoPas.init();
  Generator::CubeGauss(autoPas, 0, 0, {0., 0., 0.}, {20., 20., 20.}, 100, 5, 2, {0., 0., 0.});
  ASSERT_EQ(autoPas.getNumberOfParticles(), 100);
//    std::string CubeGeneration = "GaussTestGeneration.vtu";
//    writeVTKFile<decltype(autoPas)>(CubeGeneration, autoPas.getNumberOfParticles(), autoPas);
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    EXPECT_EQ(velocity, iter->getV());
  }
}

TEST_F(GeneratorsTest, CubeGenerator) {
  //@tood check: Fabio wieso zeigt paraview ein extra Particle an?
  auto autoPas = autopas::AutoPas<autopas::MoleculeLJ<>, FMCell>(std::cout);
  unsigned long particlesPerDim = 5;
  autoPas.setBoxMax(boxmax);
  autoPas.setBoxMin(boxmin);
  autoPas.init();
  std::array<double, 3> velocity = {0., 0., 0.};
  std::array<size_t, 3> cube = {particlesPerDim, particlesPerDim, particlesPerDim};
  Generator::CubeGrid(autoPas, 0, 0, {0., 0., 0.}, cube, .5, {0., 0., 0.});
  //  std::string CubeGeneration = "CubeGeneration.vtu";
  //  writeVTKFile<decltype(autoPas)>(CubeGeneration, autoPas.getNumberOfParticles(), autoPas);
  EXPECT_EQ(autoPas.getNumberOfParticles(), (5 * 5 * 5));
  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    EXPECT_EQ(velocity, iter->getV());
  }
}

TEST_F(GeneratorsTest, GridFillwithBoxMin) {
  auto autoPas = autopas::AutoPas<autopas::MoleculeLJ<>, FMCell>(std::cout);
  std::array<double, 3> boxmin = {5., 5., 5.};
  std::array<double, 3> boxmax = {10., 10., 10.};
  autoPas.setBoxMax(boxmax);
  autoPas.setBoxMin(boxmin);
  Molecule dummy;

  autoPas.init();
  GridGenerator::fillWithParticles(autoPas, {5, 5, 5}, 0, 0, dummy, {1, 1, 1}, boxmin);
//    std::string CubeGeneration = "FillGrid-boxMin.vtu";
//      writeVTKFile<decltype(autoPas)>(CubeGeneration, autoPas.getNumberOfParticles(), autoPas);
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    EXPECT_TRUE(autopas::utils::inBox(iter->getR(), boxmin, boxmax));
  }
}

TEST_F(GeneratorsTest, Sphere) {
  auto autoPas = autopas::AutoPas<autopas::MoleculeLJ<>, FMCell>(std::cout);
  autoPas.setBoxMax({3.1, 3.1, 3.1});
  autoPas.setBoxMin({-3., -3., -3.});
  autoPas.init();
  std::array<double, 3> velocity = {0., 0., 0.};
  Generator::Sphere(autoPas, {0., 0., 0.}, 10, .3, 0, 0);
  //  std::string SphereGeneration = "SphereGeneration.vtu";
  //  writeVTKFile<decltype(autoPas)>(SphereGeneration, autoPas.getNumberOfParticles(), autoPas);
  EXPECT_EQ(autoPas.getNumberOfParticles(), 5569);
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    EXPECT_EQ(velocity, iter->getV());
  }
}

TEST_F(GeneratorsTest, MultipleObjectGeneration) {
  auto autoPas = autopas::AutoPas<autopas::MoleculeLJ<>, FMCell>(std::cout);
  MDFlexConfig config;
  config.yamlFilename = "multipleObjectsWithMultipleTypesTest.yaml";
  parser.parseYamlFile(config);
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
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    if (ids.count(iter->getID()) != 0) {
      FAIL();  // 2 Particles with same id
    }
    ids.insert(iter->getID());
  }
  EXPECT_EQ(ids.size(), autoPas.getNumberOfParticles());
}
