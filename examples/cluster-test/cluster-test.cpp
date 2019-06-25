/**
 * This is just a quick and dirty implementation to visualize Verlet cluster. A VerletClusterLists container is filled
 * and the generated clusters are printed into a vtu files.
 * @file cluster-test.cpp
 * @date 29.11.18
 * @author nguyen
 */

#include <array>
#include <fstream>
#include <iostream>
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"
#include "autopas/autopasIncludes.h"
#include "autopas/utils/Timer.h"

#define CLUSTER_SIZE 4

void addParticles(autopas::VerletClusterLists<autopas::MoleculeLJ> &lj_system, int numParticles) {
  srand(10032);  // fixed seedpoint

  std::array<double, 3> boxMin(lj_system.getBoxMin()), boxMax(lj_system.getBoxMax());

  for (int i = 0; i < numParticles; ++i) {
    auto id = static_cast<unsigned long>(i);
    autopas::MoleculeLJ particle(RandomGenerator::randomPosition(boxMin, boxMax), {0., 0., 0.}, id);
    lj_system.addParticle(particle);
  }
}

void writeVTKFile(std::string filename, size_t numParticles, autopas::VerletClusterLists<autopas::MoleculeLJ> &cont) {
  std::ofstream vtkFile;

  size_t numClusters = numParticles / CLUSTER_SIZE;

  vtkFile.open(filename);

  vtkFile << "# vtk DataFile Version 2.0" << std::endl;
  vtkFile << "Timestep" << std::endl;
  vtkFile << "ASCII" << std::endl;
  vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
  vtkFile << "POINTS " << numParticles << " double" << std::endl;

  for (auto iter = cont.begin(); iter.isValid(); ++iter) {
    auto pos = iter->getR();
    vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
  }

  vtkFile << std::endl << "CELLS " << numClusters << " " << numClusters * (CLUSTER_SIZE + 1) << std::endl;

  // set every four points as cell
  for (size_t i = 0; i < numClusters; i++) {
    vtkFile << "4";
    for (size_t j = 0; j < CLUSTER_SIZE; j++) {
      vtkFile << " " << i * 4 + j;
    }
    vtkFile << std::endl;
  }

  // define cells as tetraeder
  vtkFile << std::endl << "CELL_TYPES " << numClusters << std::endl;
  for (size_t i = 0; i < numClusters; i++) {
    vtkFile << "10" << std::endl;
  }

  vtkFile.close();
}

int main(int argc, char *argv[]) {
  autopas::Logger::create();

  autopas::MoleculeLJ::setEpsilon(1.0);
  autopas::MoleculeLJ::setSigma(1.0);

  std::array<double, 3> boxMin({0., 0., 0.}), boxMax{};
  boxMax[0] = 0.15;
  boxMax[1] = boxMax[2] = boxMax[0] / 1.0;
  double cutoff = .03;

  int numParticles = 16;
  double skin = 0.;
  int rebuildFrequency = 1;
  if (argc == 4) {
    numParticles = atoi(argv[1]);
    boxMax[0] = boxMax[1] = boxMax[2] = atof(argv[2]);
    skin = atof(argv[3]);
  } else {
    std::cerr << "ERROR: wrong number of arguments given. " << std::endl
              << "cluster-test requires the following arguments:" << std::endl
              << "numParticles boxSize skin" << std::endl;
    exit(1);
  }

  autopas::VerletClusterLists<autopas::MoleculeLJ> cont(boxMin, boxMax, cutoff, skin * cutoff, rebuildFrequency,
                                                        CLUSTER_SIZE);

  autopas::LJFunctor<autopas::MoleculeLJ, autopas::FullParticleCell<autopas::MoleculeLJ>> func(
      cutoff, autopas::MoleculeLJ::getEpsilon(), autopas::MoleculeLJ::getSigma(), 0.0);

  addParticles(cont, numParticles);

  autopas::C01Traversal<autopas::FullParticleCell<autopas::MoleculeLJ>,
                        autopas::LJFunctor<autopas::MoleculeLJ, autopas::FullParticleCell<autopas::MoleculeLJ>>,
                        autopas::DataLayoutOption::aos, false>
      dummyTraversal({0, 0, 0}, &func);

  // iterate to rebuild
  cont.iteratePairwise(&func, &dummyTraversal);

  int newNumParticles = 0;
  for (auto iter = cont.begin(); iter.isValid(); ++iter) {
    newNumParticles++;
  }

  writeVTKFile("cluster.vtu", newNumParticles, cont);
}
