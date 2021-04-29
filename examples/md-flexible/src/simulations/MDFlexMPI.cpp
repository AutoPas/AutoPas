/**
 * @file MDFlexMPI.cpp
 * @author J. Körner
 * @date 07.04.2021
 */

#include "MDFlexMPI.h"

#include "../domainDecomposition/RegularGridDecomposition.h"

#include <array>
#include <iostream>

MDFlexMPI::MDFlexMPI(int argc, char** argv) : MDFlexSimulation(argc, argv){
  	MPI_Init(&argc, &argv);

		int processorCount;
		MPI_Comm_size(MPI_COMM_WORLD, &processorCount);
  	MPI_Comm_rank(MPI_COMM_WORLD, &_rank);

		int dimensionCount = _configuration->boxMin.value.size();
		double* globalBoxMin = &(_configuration->boxMin.value[0]);
		double* globalBoxMax = &(_configuration->boxMax.value[0]);
		RegularGridDecomposition decomposition(processorCount, dimensionCount, _rank, globalBoxMax, globalBoxMax);

}

MDFlexMPI::~MDFlexMPI() {
	MPI_Finalize();
}

void MDFlexMPI::run(){
	autopas::AutoPas<Simulation::ParticleType> autopas(std::cout);	
	_simulation->initialize(*_configuration, autopas);

	_simulation->simulate(autopas);

	_simulation->printStatistics(autopas);
}

