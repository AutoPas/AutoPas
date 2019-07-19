//
// Created by nicola on 216.06.19.
//

#include "GeneratorsTest.h"

// fehler Erstellt Umgebung von Blatt 2 Task 3 des MolSim Praktikums
void GeneratorsTest::MolSimTaskGeneration(autopas::AutoPas<Particle, FPCell> &autopas) {
  std::array<double, 3> boxMax({50., 30., 50.});

  std::array<double, 3> smallGridBoxMin({15., 15., 0});
  std::array<double, 3> bigGridBoxMin({0., 0., 0.});

  std::array<double, 3> initialVelocitySmall({0., -10., 0.});

  autopas.setBoxMin(boxmin);
  autopas.setBoxMax(boxMax);
  autopas.setCutoff(cutoff);
  autopas.init();
  Particle dummyParticle;
  // small Grid
  GridGenerator::fillWithParticlesOnRInitialV(autopas, smallGridBoxMin, initialVelocitySmall, {8, 8, 1}, dummyParticle,
                                              {.5, .5, .5}, {0.3, 0.3, 0.3});
  // Big grid
  GridGenerator::fillWithParticlesOnR(autopas, bigGridBoxMin, {40, 10, 1}, dummyParticle, {.5, .5, .5},
                                      {0.3, 0.3, 0.3});
}

TEST_F(GeneratorsTest, fillWithParticlesOnPosition) {
  auto autoPas = autopas::AutoPas<Particle, FPCell>(std::cout);
  PrintableMolecule::setEpsilon(epsilon);
  PrintableMolecule::setSigma(sigma);
  PrintableMolecule::setMass(1.0);
  Particle dummyParticle;
  std::array<double, 3> BoxPos = {2., 2., 1.};
  autoPas.setBoxMin(boxmin);
  autoPas.setBoxMax(boxmax);
  autoPas.setCutoff(cutoff);
  autoPas.init();
  GridGenerator::fillWithParticlesOnR(autoPas, BoxPos, {3, 3, 1}, dummyParticle, {1., 1., 1.}, {0.3, 0.3, 0.3});
  // GridGenerator::fillWithParticles(*autoPas,{3,3,3},dummyParticle,{0.5,0.5,.5},{0.3,0.3,0.3});
  std::cout << "Number of particles generated " << autoPas.getNumberOfParticles() << std::endl;
  ASSERT_TRUE(true);
}

TEST_F(GeneratorsTest, Behavior) {
  auto autoPas = autopas::AutoPas<Particle, FPCell>(std::cout);
  PrintableMolecule::setEpsilon(epsilon);
  PrintableMolecule::setSigma(sigma);
  PrintableMolecule::setMass(1.0);
  MolSimTaskGeneration(autoPas);
  // for GRID generator:
  //  int particlePerDim = 5;
  //  double particleSpacing = 0.5;
  // initContainerGrid(*autoPas,particlePerDim,particleSpacing);
  // Uniform generator
  // initContainerUniform(*autoPas,5.,125);
  // Gauß generator
  // initContainerGauss(*autoPas,5.,125,5,2);
  //  std::cout << "Number of particles generated " << autoPas.getNumberOfParticles() << std::endl;
  //    for (auto iter = autoPas->getContainer()->begin() ; iter.isValid(); ++iter) {
  //        cout << iter->toString() << endl;
  //    }

  // print State -- ugly code , I know
  size_t numParticles = autoPas.getNumberOfParticles();
  string filename = "VtkTestOutput.vtu";
  std::ofstream vtkFile;
  vtkFile.open(filename);
  vtkFile << "# vtk DataFile Version 2.0" << std::endl;
  vtkFile << "Timestep" << std::endl;
  vtkFile << "ASCII" << std::endl;
  vtkFile << "DATASET STRUCTURED_GRID" << std::endl;
  vtkFile << "DIMENSIONS 1 1 1" << std::endl;
  vtkFile << "POINTS " << numParticles << " double" << std::endl;
  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    auto pos = iter->getR();
    vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
  }
  vtkFile.close();
  ASSERT_TRUE(true);
}
TEST_F(GeneratorsTest, MolSimTask) {
  auto autoPas = autopas::AutoPas<Particle, FPCell>(std::cout);
  PrintableMolecule::setEpsilon(epsilon);
  PrintableMolecule::setSigma(sigma);
  PrintableMolecule::setMass(1.0);
  MolSimTaskGeneration(autoPas);
  std::cout << "Number of particles generated " << autoPas.getNumberOfParticles() << std::endl;
  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    // std::cout << iter->toString() << std::endl;
  }
  double particleD = 0.01;
  int iterations = 0;
  // iterationen beginnend
  TimeDiscretization<decltype(autoPas)> td(particleD);
  // domain vorbeireiten: -Force initialisieren
  autoPas.iteratePairwise(functor);
  while (iterations < 10) {
    td.VSCalculateX(autoPas);
    autoPas.iteratePairwise(functor);
    td.VSCalculateV(autoPas);
    iterations++;
  }
  ASSERT_TRUE(true);
}
