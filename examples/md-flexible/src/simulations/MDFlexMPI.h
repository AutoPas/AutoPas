/**
 * @file MDFlexMPI.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include "MDFlexSimulation.h"

#include "../domainDecomposition/RegularGridDecomposition.h"

#include <mpi.h>

class MDFlexMPI : MDFlexSimulation {
	public:
		MDFlexMPI(int argc, char** argv);
		~MDFlexMPI();

		void run() override;

	private:
		std::unique_ptr<autopas::AutoPas<ParticleType>> _autoPasContainer;
  	std::unique_ptr<ParticlePropertiesLibraryType> _particlePropertiesLibrary;

		std::unique_ptr<RegularGridDecomposition> _domainDecomposition;
		std::array<unsigned int, 3> _processorCoordinates;

		void updateParticles(const int iterationsPerSupErstep);
		void executeSuperstep(const int iterationsPerSuperstep);
		void sendEmigrantsToNeighbours(std::vector<ParticleType> &emigrants);
};

