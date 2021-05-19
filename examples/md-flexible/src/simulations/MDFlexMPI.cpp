/**
 * @file MDFlexMPI.cpp
 * @author J. Körner
 * @date 07.04.2021
 */

#include "MDFlexMPI.h"

#include "autopas/utils/ArrayMath.h"

#include <algorithm>
#include <iostream>

MDFlexMPI::MDFlexMPI(int dimensionCount, int argc, char **argv) {
	MDFlexSimulation::initialize(dimensionCount, argc, argv);
}

void MDFlexMPI::run(){
	// @todo: make variable part of MDFlexConfig
	int iterationsPerSuperstep = 10;
	int remainingIterations = _configuration->iterations.value;

	for (int i = 0; i < _configuration->iterations.value; i+=iterationsPerSuperstep){
		executeSuperstep(iterationsPerSuperstep);
	}
}

void MDFlexMPI::executeSuperstep(const int iterationsPerSuperstep){
	_domainDecomposition->exchangeHaloParticles(_autoPasContainer);

	//updateParticles(iterationsPerSuperstep);

	// todo: Not sure if this synchronization is required. Check performance with and without it.
	_domainDecomposition->synchronizeDomains();

	_domainDecomposition->exchangeMigratingParticles(_autoPasContainer);
}

void MDFlexMPI::updateParticles(const int iterationsPerSuperstep){
 	// @todo: update positions
	// I don't think we need to update boundary conditions as with a single rank

	// @todo: update forces

	// @todo: update velocities
}

void MDFlexMPI::initializeDomainDecomposition(int &dimensionCount){
	std::cout << "InitializeDomainDecomposition" << std::endl;
	std::vector<double> boxMin(_configuration->boxMin.value.begin(), _configuration->boxMin.value.end());
	std::vector<double> boxMax(_configuration->boxMax.value.begin(), _configuration->boxMax.value.end());
	
	_domainDecomposition = std::make_shared<RegularGrid>(_argc, _argv, dimensionCount, boxMin, boxMax);

	std::vector<double> localBoxMin = _domainDecomposition->getLocalBoxMin();
	std::vector<double> localBoxMax = _domainDecomposition->getLocalBoxMax();
	for (int i = 0; i < localBoxMin.size(); ++i){
		_configuration->boxMin.value[i] = localBoxMin[i];
		_configuration->boxMax.value[i] = localBoxMax[i];
	}
	std::cout << "Initialized Local Box" << std::endl;
	std::cout << autopas::utils::ArrayUtils::to_string(localBoxMin) << std::endl;
	std::cout << autopas::utils::ArrayUtils::to_string(localBoxMax) << std::endl;
}

