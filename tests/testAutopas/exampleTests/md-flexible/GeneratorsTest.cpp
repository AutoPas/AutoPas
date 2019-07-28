//
// Created by nicola on 216.06.19.
//

#include "GeneratorsTest.h"
#include "autopas/utils/inBox.h"
TEST_F(GeneratorsTest, Gauss) {
  auto autoPas = autopas::AutoPas<Particle, FPCell>(std::cout);
  Particle dummyParticle;
  autoPas.setBoxMax({20., 20., 20.});
  autoPas.setBoxMin({0., 0., 0.});
  std::array<double, 3> velocity = {0., 0., 0.};
  autoPas.init();
  Generator::CubeGauss(autoPas, {0., 0., 0.}, {20., 20., 20.}, 100, 5, 2, {0., 0., 0.});
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
  auto autoPas = autopas::AutoPas<Particle, FPCell>(std::cout);
  unsigned long particlesPerDim = 5;
  autoPas.setBoxMax(boxmax);
  autoPas.setBoxMin(boxmin);
  autoPas.init();
  std::array<double, 3> velocity = {0., 0., 0.};
  std::array<size_t, 3> cube = {particlesPerDim, particlesPerDim, particlesPerDim};
  Generator::CubeGrid(autoPas, {0., 0., 0.}, cube, .5, {0., 0., 0.});
  //  std::string CubeGeneration = "CubeGeneration.vtu";
  //  writeVTKFile<decltype(autoPas)>(CubeGeneration, autoPas.getNumberOfParticles(), autoPas);
  EXPECT_EQ(autoPas.getNumberOfParticles(), (5 * 5 * 5));
  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    EXPECT_EQ(velocity, iter->getV());
  }
}

TEST_F(GeneratorsTest, GridFillwithBoxMin) {
  auto autoPas = autopas::AutoPas<Particle, FPCell>(std::cout);
  std::array<double, 3> boxmin = {5., 5., 5.};
  std::array<double, 3> boxmax = {10., 10., 10.};
  autoPas.setBoxMax(boxmax);
  autoPas.setBoxMin(boxmin);
  Particle dummy;

  autoPas.init();
  GridGenerator::fillWithParticles(autoPas, {5, 5, 5}, dummy, {1, 1, 1}, boxmin, {0., 0., 0.});
//    std::string CubeGeneration = "FillGrid-BoxMin.vtu";
//      writeVTKFile<decltype(autoPas)>(CubeGeneration, autoPas.getNumberOfParticles(), autoPas);
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    EXPECT_TRUE(autopas::utils::inBox(iter->getR(), boxmin, boxmax));
  }
}

TEST_F(GeneratorsTest, Sphere) {
  auto autoPas = autopas::AutoPas<Particle, FPCell>(std::cout);
  autoPas.setBoxMax({3., 3., 3.});
  autoPas.setBoxMin({-3., -3., -3.});
  autoPas.init();
  std::array<double, 3> velocity = {0., 0., 0.};
  Generator::Sphere(autoPas, {0., 0., 0.}, 10, .3, 0);
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
  auto autoPas = autopas::AutoPas<Particle, FPCell>(std::cout);
  autoPas.setBoxMax(parser.getBoxMax());
  autoPas.setBoxMin(parser.getBoxMin());
  autoPas.init();
  std::array<double, 3> velocity = {0., 0., 0.};
  // parses the multiple Objects input of "testParsing.yaml" and generates a VTK File from the Input
  auto CubeGrid(parser.getCubeGrid());
  auto CubeGauss(parser.getCubeGauss());
  auto CubeUniform(parser.getCubeUniform());
  auto Sphere(parser.getSphere());
  for (auto C : CubeGrid) {
    Generator::CubeGrid(autoPas, C.getBoxMin(), C.getParticlesPerDim(), C.getParticleSpacing(), C.getVelocity());
  }
  for (auto C : CubeGauss) {
    Generator::CubeGauss(autoPas, C.getBoxMin(), C.getBoxMax(), C.getNumParticles(), C.getDistributionMean(),
                         C.getDistributionStdDev(), C.getVelocity());
  }
  for (auto C : CubeUniform) {
    Generator::CubeRandom(autoPas, C.getBoxMin(), C.getBoxMax(), C.getNumParticles(), C.getVelocity());
  }
  for (auto S : Sphere) {
    Generator::Sphere(autoPas, S.getCenter(), S.getRadius(), S.getParticleSpacing(), S.getId(), S.getVelocity());
  }
  // to see output:
  //  std::string SphereGeneration = "MultipleGeneration.vtu";
  //  writeVTKFile<decltype(autoPas)>(SphereGeneration, autoPas.getNumberOfParticles(), autoPas);
  // checked VTK File
  EXPECT_EQ(parser.particlesTotal(), autoPas.getNumberOfParticles());
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    EXPECT_EQ(velocity, iter->getV());
  }
}
