//
// Created by nicola on 216.06.19.
//

#include "GeneratorsTest.h"

// fehler Erstellt Umgebung von Blatt 2 Task 3 des MolSim Praktikums
void GeneratorsTest::MolSimTaskGeneration(autopas::AutoPas<Particle, FPCell> &autopas) {
  std::array<double, 3> boxMax({50., 30., 50.});

  std::array<double, 3> smallGridBoxMin({15., 15., 0});
  //  std::array<double, 3> bigGridBoxMin({0., 0., 0.});

  std::array<double, 3> initialVelocitySmall({0., -10., 0.});

  autopas.setBoxMin(boxmin);
  autopas.setBoxMax(boxMax);
  autopas.setCutoff(cutoff);
  autopas.init();
  Particle dummyParticle;
  // small Grid
  GridGenerator::fillWithParticlesOnR(autopas, smallGridBoxMin, {8, 8, 1}, dummyParticle, {.5, .5, .5}, {0.3, 0.3, 0.3},
                                      initialVelocitySmall);
  // Big grid
  // initialize
}

TEST_F(GeneratorsTest, Gauss) {
  auto autoPas = autopas::AutoPas<Particle, FPCell>(std::cout);
  PrintableMolecule::setEpsilon(epsilon);
  PrintableMolecule::setSigma(sigma);
  PrintableMolecule::setMass(1.0);
  Particle dummyParticle;
  Generator::CubeGauss(autoPas, {5., 5., 5.}, 100, 5, 2);
  ASSERT_EQ(autoPas.getNumberOfParticles(), 100);
}

TEST_F(GeneratorsTest, CubeGenerator) {
  //@tood check: Fabio wieso zeigt paraview ein extra Particle an?
  auto autoPas = autopas::AutoPas<Particle, FPCell>(std::cout);
  PrintableMolecule::setEpsilon(epsilon);
  PrintableMolecule::setSigma(sigma);
  PrintableMolecule::setMass(1.0);
  unsigned long particlesPerDim = 5;
  std::array<size_t, 3> cube = {particlesPerDim, particlesPerDim, particlesPerDim};
  Generator::CubeGrid(autoPas, cube, .5);
  string CubeGeneration = "CubeGeneration.vtu";
  writeVTKFile<decltype(autoPas)>(CubeGeneration, autoPas.getNumberOfParticles(), autoPas);
  EXPECT_EQ(autoPas.getNumberOfParticles(), (5 * 5 * 5));
}

TEST_F(GeneratorsTest, MolSimTask) {
  auto autoPas = autopas::AutoPas<Particle, FPCell>(std::cout);
  PrintableMolecule::setEpsilon(epsilon);
  PrintableMolecule::setSigma(sigma);
  PrintableMolecule::setMass(1.0);
  MolSimTaskGeneration(autoPas);
  std::cout << "Number of particles generated " << autoPas.getNumberOfParticles() << std::endl;  // 64
  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    // std::cout << iter->toString() << std::endl;
  }
  // iterationen beginnend
  //@todo schauen was man hier machen kann zum testen: VTK File ausgabe wäre eine idee
  ASSERT_TRUE(true);
}

TEST_F(GeneratorsTest, Sphere) {
  auto autoPas = autopas::AutoPas<Particle, FPCell>(std::cout);
  PrintableMolecule::setEpsilon(epsilon);
  PrintableMolecule::setSigma(sigma);
  PrintableMolecule::setMass(1.0);
  Generator::Sphere(autoPas, {0., 0., 0.}, 10, .3, 0);
  cout << "Number of particles generated " << autoPas.getNumberOfParticles() << endl;
  string SphereGeneration = "SphereGeneration.vtu";
  writeVTKFile<decltype(autoPas)>(SphereGeneration, autoPas.getNumberOfParticles(), autoPas);
  cout << "Vtk output in Autopas/tests/testAutopas" << endl;
  // check if no Particle is outside of radius range:
}