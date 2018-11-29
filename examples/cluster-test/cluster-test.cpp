/**
 * @file lj-traversals.cpp
 * @date 04.11.18
 * @author nguyen
 */

#include <array>
#include <iostream>
#include <fstream>
#include "../md/mdutils.h"
#include "autopas/autopasIncludes.h"
#include "autopas/utils/Timer.h"

#define CLUSTER_SIZE 4

void addParticles(
    autopas::VerletClusterLists<PrintableMolecule> &lj_system,
    int numParticles) {

  srand(10032);  // fixed seedpoint

  std::array<double, 3> boxMin(lj_system.getBoxMin()), boxMax(lj_system.getBoxMax());

  for (int i = 0; i < numParticles; ++i) {
    auto id = static_cast<unsigned long>(i);
    PrintableMolecule particle(randomPosition(boxMin, boxMax), {0., 0., 0.}, id);
    lj_system.addParticle(particle);
  }
}

void writeVTKFile(string filename, size_t numParticles,
                  autopas::VerletClusterLists<PrintableMolecule> &cont) {
  std::ofstream vtkFile;

  size_t numClusters = numParticles / CLUSTER_SIZE;

  vtkFile.open(filename);

  vtkFile << "# vtk DataFile Version 2.0" << endl;
  vtkFile << "Timestep" << endl;
  vtkFile << "ASCII" << endl;
  vtkFile << "DATASET UNSTRUCTURED_GRID" << endl;
  vtkFile << "POINTS " << numParticles << " double" << endl;

  for (auto iter = cont.begin(); iter.isValid(); ++iter) {
    auto pos = iter->getR();
    vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << endl;
  }

  vtkFile << endl << "CELLS " << numClusters << " " << numClusters * (CLUSTER_SIZE + 1) <<  endl;

  // set every four points as cell
  for (size_t i = 0; i < numClusters; i++) {
    vtkFile << "4";
    for (size_t j = 0; j < CLUSTER_SIZE; j++) {
	vtkFile << " " << i * 4 + j;
    }
    vtkFile << endl;
  }

  // define cells as tetraeder
  vtkFile << endl << "CELL_TYPES " << numClusters  <<  endl;
  for (size_t i = 0; i < numClusters; i++) {
    vtkFile << "10" << endl;
  }

  vtkFile.close();
}

int main(int argc, char *argv[]) {
  autopas::Logger::create();

  PrintableMolecule::setEpsilon(1.0);
  PrintableMolecule::setSigma(1.0);

  std::array<double, 3> boxMin({0., 0., 0.}), boxMax{};
  boxMax[0] = 0.15;
  boxMax[1] = boxMax[2] = boxMax[0] / 1.0;
  double cutoff = .03;


  autopas::LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>> func;

  int numParticles = 16;
  int numIterations = 1;
  bool useNewton3 = false;
  double skin = 0.;
  int rebuildFrequency = 1;
  if (argc == 4) {
    numParticles = atoi(argv[1]);
    boxMax[0] = boxMax[1] = boxMax[2] = atof(argv[2]);
    skin = atof(argv[3]);
  } else {
    std::cerr << "ERROR: wrong number of arguments given. " << std::endl
              << "cluster-test requires the following arguments:" << std::endl
              << "numParticles boxSize skin"
              << std::endl;
    exit(1);
  }

  autopas::VerletClusterLists<PrintableMolecule> cont(boxMin, boxMax, cutoff, skin * cutoff,
                                                                           rebuildFrequency, CLUSTER_SIZE);

  addParticles(cont, numParticles);

  C01Traversal<FullParticleCell<PrintableMolecule>, autopas::LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>>, false>
              dummyTraversal({0, 0, 0}, &func);

  // iterate to rebuild
  cont.iteratePairwiseAoS(&func, &dummyTraversal, useNewton3);

  int newNumParticles = 0;
  for (auto iter = cont.begin(); iter.isValid(); ++iter) {
    newNumParticles++;
  }

  writeVTKFile("cluster.vtu", newNumParticles, cont);
}
