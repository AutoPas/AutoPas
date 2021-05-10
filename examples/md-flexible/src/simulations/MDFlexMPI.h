/**
 * @file MDFlexMPI.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include "MDFlexSimulation.h"


#include <mpi.h>

class MDFlexMPI : MDFlexSimulation {
	public:
		MDFlexMPI(int argc, char** argv);
		~MDFlexMPI();

		void run() override;

	private:
		std::unique_ptr<RegularGridDecomposition> _domainDecomposition;

		void updateParticles(const int iterationsPerSupErstep);
		void executeSuperstep(const int iterationsPerSuperstep);
		void sendEmigrantsToNeighbours(std::vector<ParticleType> &emigrants);
};

