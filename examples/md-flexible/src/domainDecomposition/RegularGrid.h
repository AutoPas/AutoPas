/**
 * @file RegularGrid.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#include "DomainDecomposition.h"

#include "../TypeDefinitions.h"
#include "autopas/AutoPas.h"

#include "mpi.h"
#include <list>
#include <memory>

class RegularGrid final : protected DomainDecomposition {
	public:
		RegularGrid(int argc, char** argv, const int &dimensionCount, const std::vector<double> &globalBoxMin,
			const std::vector<double> &globalBoxMax);
		~RegularGrid();
	
		using SharedAutoPasContainer = std::shared_ptr<autopas::AutoPas<ParticleType>>;

	 	void update() override;
		const int getDimensionCount() override { return _dimensionCount; }
		std::vector<double> getLocalBoxMin() { return _localBoxMin; }	
		std::vector<double> getLocalBoxMax() { return _localBoxMax; }	

		int convertIdToIndex(const std::vector<int> &domainIndex);
		void sendDataToNeighbour(std::vector<char> sendBuffer, const int &neighbour);
		void receiveDataFromNeighbour(const int &neighbour, std::vector<char> &dataBuffer);
		void waitForSendRequests();
		void synchronizeDomains();
		void exchangeHaloParticles(SharedAutoPasContainer &autoPasContainer);
		void exchangeMigratingParticles(SharedAutoPasContainer &autoPasContainer);

	private:
		// Global data
		int _dimensionCount;
		int _subdomainCount;
		std::vector<double> _globalBoxMin;
		std::vector<double> _globalBoxMax;
		std::vector<int> _decomposition;
		MPI_Comm _communicator;

		// Domain specific data
		int _domainIndex;
		std::vector<int> _domainId;
		std::vector<int> _neighbourDomainIndices;
  	std::vector<double> _localBoxMin;
  	std::vector<double> _localBoxMax;
		std::list<MPI_Request> _sendRequests;
		std::list<std::vector<char>> _sendBuffers;

		void initializeDecomposition();
		void initializeMPICommunicator();
		void initializeLocalDomain();
		void initializeGlobalBox(const std::vector<double> &globalBoxMin, const std::vector<double> &globalBoxMax);
		void initializeLocalBox();
		void initializeNeighbourIds();

		void updateLocalBox();
		void sendParticles(std::vector<ParticleType> &particles, int &receiver);
		void receiveParticles(std::vector<ParticleType> &receivedParticles, int &source);
};
