//
// Created by nicola on 16.06.19.
//

#include <gtest/gtest.h>
#include <math.h>
#include <vector>
#include "../../../../examples/md-flexible/MDFlexParser.h"
#include "../../../../examples/md-flexible/PrintableMolecule.h"
#include "../../../../examples/md-flexible/TimeDiscretization.h"
#include "../../../../src/autopas/utils/ArrayMath.h"
#include "../../testingHelpers/GaussianGenerator.h"
#include "../../testingHelpers/GridGenerator.h"
#include "../../testingHelpers/RandomGenerator.h"
#include "autopas/AutoPas.h"
#include "testingHelpers/commonTypedefs.h"

using namespace std;
using namespace autopas;

template <class AutoPasTemplate>
void writeVTKFile(int iteration, size_t numParticles, AutoPasTemplate &autopas) {
  string filename = "VtkOutput";
  stringstream strstr;
  strstr << filename << "_" << setfill('0') << setw(4) << iteration << ".vtu";
  // string path = "./vtk";
  std::ofstream vtkFile;
  vtkFile.open(strstr.str());

  vtkFile << "# vtk DataFile Version 2.0" << endl;
  vtkFile << "Timestep" << endl;
  vtkFile << "ASCII" << endl;
  vtkFile << "DATASET STRUCTURED_GRID" << endl;
  vtkFile << "DIMENSIONS 1 1 1" << endl;
  vtkFile << "POINTS " << numParticles << " double" << endl;

  for (auto iter = autopas->begin(); iter.isValid(); ++iter) {
    auto pos = iter->getR();
    vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << endl;
  }

  vtkFile.close();
}

void initContainerGrid(autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                       size_t particlesPerDim, double particelSpacing) {
  std::array<double, 3> boxMin({0., 0., 0.});
  std::array<double, 3> boxMax(
      {(particlesPerDim)*particelSpacing, (particlesPerDim)*particelSpacing, (particlesPerDim)*particelSpacing});

  autopas.setBoxMin(boxMin);
  autopas.setBoxMax(boxMax);

  autopas.init();

  PrintableMolecule dummyParticle;
  GridGenerator::fillWithParticles(autopas, {particlesPerDim, particlesPerDim, particlesPerDim}, dummyParticle,
                                   {particelSpacing, particelSpacing, particelSpacing},
                                   {particelSpacing / 2, particelSpacing / 2, particelSpacing / 2});
}

void initContainerGauss(autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                        double boxLength, size_t numParticles, double distributionMean, double distributionStdDev) {
  std::array<double, 3> boxMin({0., 0., 0.});
  std::array<double, 3> boxMax({boxLength, boxLength, boxLength});

  autopas.setBoxMin(boxMin);
  autopas.setBoxMax(boxMax);

  autopas.init();

  PrintableMolecule dummyParticle;
  GaussianGenerator::fillWithParticles(autopas, numParticles, dummyParticle, distributionMean, distributionStdDev);
}

void initContainerUniform(autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                          double boxLength, size_t numParticles) {
  std::array<double, 3> boxMin({0., 0., 0.});
  std::array<double, 3> boxMax({boxLength, boxLength, boxLength});

  autopas.setBoxMin(boxMin);
  autopas.setBoxMax(boxMax);

  autopas.init();

  PrintableMolecule dummyParticle;
  RandomGenerator::fillWithParticles(autopas, dummyParticle, numParticles);
}

// FRAGE-FABIO
// wieso kann ich nicht Particle und ParticleCell im konstruktor behalten-> wenn ichs behalte-> no matching constructor
// fehler Erstellt Umgebung von Blatt 2 Task 3 des MolSim Praktikums
void MolSimTaskGeneration(autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas) {
  std::array<double, 3> boxMin({0., 0., 0.});
  std::array<double, 3> boxMax({50., 30., 50.});

  std::array<double, 3> smallGridBoxMin({15., 15., 0});
  std::array<double, 3> bigGridBoxMin({0., 0., 0.});

  std::array<double, 3> smallGridBoxMax({23., 23., 0});
  std::array<double, 3> bigGridBoxMax({40., 10., 0});

  std::array<double, 3> initialVelocitySmall({0., -10., 0.});

  autopas.setBoxMin(boxMin);
  autopas.setBoxMax(boxMax);
  autopas.setCutoff(1.5);
  autopas.init();
  PrintableMolecule dummyParticle;
  // small Grid
  GridGenerator::fillWithParticlesOnRInitialV(autopas, smallGridBoxMin, initialVelocitySmall, {8, 8, 1}, dummyParticle,
                                              {.5, .5, .5}, {0.3, 0.3, 0.3});
  // Big grid
  GridGenerator::fillWithParticlesOnR(autopas, bigGridBoxMax, {40, 10, 1}, dummyParticle, {.5, .5, .5},
                                      {0.3, 0.3, 0.3});
}

TEST(Generater, fillWithParticlesOnPosition) {
  auto *autoPas = new autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>(std::cout);
  double epsilon = 5.0;
  double sigma = 1.0;
  double cutoff = 2;
  array<double, 3> boxmin = {0., 0., 0.};
  array<double, 3> boxmax = {5., 5., 5.};
  PrintableMolecule::setEpsilon(epsilon);
  PrintableMolecule::setSigma(sigma);
  PrintableMolecule::setMass(1.0);
  PrintableMolecule dummyParticle;
  std::array<double, 3> BoxPos = {2., 2., 1.};
  autoPas->setBoxMin(boxmin);
  autoPas->setBoxMax(boxmax);
  autoPas->setCutoff(cutoff);
  autoPas->init();
  GridGenerator::fillWithParticlesOnR(*autoPas, BoxPos, {3, 3, 1}, dummyParticle, {1., 1., 1.}, {0.3, 0.3, 0.3});
  // GridGenerator::fillWithParticles(*autoPas,{3,3,3},dummyParticle,{0.5,0.5,.5},{0.3,0.3,0.3});
  cout << "Number of particles generated " << autoPas->getNumberOfParticles() << endl;

  delete autoPas;
  ASSERT_TRUE(true);
}

TEST(Generater, Behavior) {
  auto *autoPas = new autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>(std::cout);
  double epsilon = 5.0;
  double sigma = 1.0;
  double cutoff = 2;
  array<double, 3> boxmin = {0., 0., 0.};
  array<double, 3> boxmax = {5., 5., 5.};
  PrintableMolecule::setEpsilon(epsilon);
  PrintableMolecule::setSigma(sigma);
  PrintableMolecule::setMass(1.0);
  MolSimTaskGeneration(*autoPas);
  //    autoPas->init();
  // for GRID generator:
  int particlePerDim = 5;
  double particleSpacing = 0.5;
  // initContainerGrid(*autoPas,particlePerDim,particleSpacing);
  // Uniform generator
  // initContainerUniform(*autoPas,5.,125);
  // Gauß generator
  // initContainerGauss(*autoPas,5.,125,5,2);

  MolSimTaskGeneration(*autoPas);
  cout << "Number of particles generated " << autoPas->getNumberOfParticles() << endl;
  //    for (auto iter = autoPas->getContainer()->begin() ; iter.isValid(); ++iter) {
  //        cout << iter->toString() << endl;
  //    }

  // print State -- ugly code , I know
  size_t numParticles = autoPas->getNumberOfParticles();
  string filename = "VtkTestOutput.vtu";
  std::ofstream vtkFile;
  vtkFile.open(filename);
  vtkFile << "# vtk DataFile Version 2.0" << endl;
  vtkFile << "Timestep" << endl;
  vtkFile << "ASCII" << endl;
  vtkFile << "DATASET STRUCTURED_GRID" << endl;
  vtkFile << "DIMENSIONS 1 1 1" << endl;
  vtkFile << "POINTS " << numParticles << " double" << endl;
  for (auto iter = autoPas->begin(); iter.isValid(); ++iter) {
    auto pos = iter->getR();
    vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << endl;
  }
  vtkFile.close();

  delete autoPas;
  ASSERT_TRUE(true);
}
TEST(Generator, MolSimTask) {
  auto *autoPas = new autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>>(std::cout);
  double epsilon = 5.0;
  double sigma = 1.0;
  double cutoff = 0.5;
  std::array<double, 3> boxmin({0., 0., 0.});
  std::array<double, 3> boxmax({50., 30., 50.});
  PrintableMolecule::setEpsilon(epsilon);
  PrintableMolecule::setSigma(sigma);
  PrintableMolecule::setMass(1.0);
  MolSimTaskGeneration(*autoPas);
  // initContainerGrid(*autoPas,20,.5);
  cout << "Number of particles generated " << autoPas->getNumberOfParticles() << endl;
  for (auto iter = autoPas->getContainer()->begin(); iter.isValid(); ++iter) {
    cout << iter->toString() << endl;

    double particleD = 0.01;
    int iterations = 0;
    // iterationen beginnend
    TimeDiscretization<decltype(autoPas)> td(particleD);
    auto *functor =
        new autopas::LJFunctor<PrintableMolecule, autopas::ParticleCell<PrintableMolecule>,
                               autopas::FunctorN3Modes::Both, true>(cutoff, epsilon, sigma, 0.0, boxmin, boxmax, true);
    // domain vorbeireiten: -Force initialisieren
    autoPas->iteratePairwise(functor);
    for (auto iter = autoPas->getContainer()->begin(); iter.isValid(); ++iter) {
      cout << iter->toString() << endl;
    }

    writeVTKFile<decltype(autoPas)>(iterations, autoPas->getNumberOfParticles(), autoPas);
    while (iterations < 10) {
      td.VSCalculateX(autoPas);
      autoPas->iteratePairwise(functor);
      td.VSCalculateV(autoPas);
      iterations++;
      writeVTKFile<decltype(autoPas)>(iterations, autoPas->getNumberOfParticles(), autoPas);
    }

    delete autoPas;
    delete functor;
    ASSERT_TRUE(true);
  }
}
