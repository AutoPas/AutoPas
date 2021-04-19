/**
 * @file MDFlexMPI.cpp
 * @author J. Körner
 * @date 07.04.2021
 */

#include "MDFlexMPI.h"

#include <array>
#include <mpi.h>

MDFlexMPI::MDFlexMPI(int argc, char** argv){
  	MPI_Init(&argc, &argv);
  	MPI_Comm_rank(MPI_COMM_WORLD, &_rank);

		int processorCount;
		MPI_Comm_size(MPI_COMM_WORLD, processorCount);

		std::array<unsigned int, 3> processorBlockDimensions;

		processorBlockDimension[0] = std::ceil(std::pow(processorCount, 1 / 3);
		processorBlockDimension[1] = processorBlockDimension[0];
		processorBlockDimension[2] = processorBlockCount / (processorBlockDimension[0] * processorBlockDimension[1]);

  	MPI_Cart_get(MPI_COMM_WORLD, 3, processorBlockDimensions, false, _processorCoordinates);
}

void MDFlexMPI::run(){
	autopas::AutoPas<Simulation::ParticleType> autopas();	
	_simulation->initialize(*_configuration, autopas);

	_simulation->simulate(autopas)

	_simulation->printStatistics(autopas)
}

