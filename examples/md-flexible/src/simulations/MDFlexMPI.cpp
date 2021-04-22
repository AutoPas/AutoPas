/**
 * @file MDFlexMPI.cpp
 * @author J. Körner
 * @date 07.04.2021
 */

#include "MDFlexMPI.h"

#include "../domainDecomposition/RegularGridDecomposition.h"

#include <array>
#include <mpi.h>

#include <iostream>

MDFlexMPI::MDFlexMPI(int argc, char** argv) : MDFlexSimulation(argc, argv){
  	MPI_Init(&argc, &argv);
  	MPI_Comm_rank(MPI_COMM_WORLD, &_rank);

		int processorCount;
		MPI_Comm_size(MPI_COMM_WORLD, &processorCount);
}

void MDFlexMPI::run(){
	autopas::AutoPas<Simulation::ParticleType> autopas();	
	_simulation->initialize(*_configuration, autopas);

	_simulation->simulate(autopas);

	_simulation->printStatistics(autopas);
}

