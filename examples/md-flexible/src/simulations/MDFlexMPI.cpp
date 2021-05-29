/**
 * @file MDFlexMPI.cpp
 * @author J. Körner
 * @date 07.04.2021
 */

#include "MDFlexMPI.h"

#include <algorithm>
#include <iostream>

#include "autopas/utils/ArrayMath.h"

MDFlexMPI::MDFlexMPI(int dimensionCount, int argc, char **argv) {
  this->initialize(dimensionCount, argc, argv);
}

void MDFlexMPI::run() {
  // @todo: make variable part of MDFlexConfig
  int iterationsPerSuperstep = 10;
  int remainingIterations = _configuration->iterations.value;

  for (int i = 0; i < _configuration->iterations.value; i += iterationsPerSuperstep) {
    executeSuperstep(iterationsPerSuperstep);
  }
}

void MDFlexMPI::executeSuperstep(const int iterationsPerSuperstep) {
  _domainDecomposition->exchangeHaloParticles(_autoPasContainer);

  for (int i = 1; i < iterationsPerSuperstep; ++i) {
    updateParticles();
  }

  _domainDecomposition->exchangeMigratingParticles(_autoPasContainer);
}

void MDFlexMPI::updateParticles() {
  updatePositions();
  updateForces();
  updateVelocities();
  updateThermostat();
}

void MDFlexMPI::initializeDomainDecomposition(int &dimensionCount) {
  std::vector<double> boxMin(_configuration->boxMin.value.begin(), _configuration->boxMin.value.end());
  std::vector<double> boxMax(_configuration->boxMax.value.begin(), _configuration->boxMax.value.end());

  _domainDecomposition = std::make_shared<RegularGrid>(dimensionCount, boxMin, boxMax);

  _domainDecomposition->setHaloWidth(_configuration->cutoff.value + _configuration->verletSkinRadius.value);

  std::vector<double> localBoxMin = _domainDecomposition->getLocalBoxMin();
  std::vector<double> localBoxMax = _domainDecomposition->getLocalBoxMax();
  for (int i = 0; i < localBoxMin.size(); ++i) {
    _configuration->boxMin.value[i] = localBoxMin[i];
    _configuration->boxMax.value[i] = localBoxMax[i];
  }

}
