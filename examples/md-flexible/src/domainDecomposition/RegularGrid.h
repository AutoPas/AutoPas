/**
 * @file RegularGrid.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#include "DomainDecomposition.h"

#include "mpi.h"

class RegularGrid : protected DomainDecomposition {
	public:
		RegularGrid(int argc, char** argv, const int &dimensionCount, const std::vector<double> &globalBoxMin,
			const std::vector<double> &globalBoxMax);
		~RegularGrid();
	
	 	void update() override;
		void exchangeHaloData() override;
		const int getDimensionCount() override { return _dimensionCount; }
		std::vector<double> getLocalBoxMin() { return _localBoxMin; }	
		std::vector<double> getLocalBoxMax() { return _localBoxMax; }	

		const MPI_Comm getCommunicator() { return _communicator; }

		int convertIdToIndex(const std::vector<int> &domainIndex);

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

		void initializeDecomposition();
		void initializeMPICommunicator();
		void initializeLocalDomain();
		void initializeGlobalBox(const std::vector<double> &globalBoxMin, const std::vector<double> &globalBoxMax);
		void initializeLocalBox();
		void initializeNeighbourIds();

		void updateLocalBox();
};
