/**
 * @file MDFlexMPI.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include "MDFlexSimulation.h"

#include "../domainDecomposition/RegularGrid.h"

#include <mpi.h>

class MDFlexMPI : protected MDFlexSimulation {
	public:
		MDFlexMPI(int dimensionCount, int argc, char **argv);
		~MDFlexMPI() = default;

		void run() override;
		void initializeDomainDecomposition(int &dimensionCount) override;

	private:
		std::shared_ptr<RegularGrid> _domainDecomposition;

		void updateParticles(const int iterationsPerSupErstep);
		void executeSuperstep(const int iterationsPerSuperstep);

		void sendParticles(std::vector<ParticleType> &particles, int &receiver);
		void receiveParticles(std::vector<ParticleType> &receivedParticles, int &source);
};

