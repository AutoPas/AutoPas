/**
 * @file RegularGridDecomposition.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#include "DomainDecomposition.h"

class RegularGridDecomposition : protected DomainDecomposition<std::vector<int>, std::vector<int>> {
	public:
		RegularGridDecomposition(const unsigned int &dimensionCount, const double* globalBoxMin, const double* golbalBoxMax);
		~RegularGridDecomposition() = default;

	 	void update() override;

		const int getDimensionCount() override { return _dimensionCount; }
		const MPI_Comm getCommunicator() override { return _communicator; }

		std::vector<double> getLocalBoxMin() { return _localBoxMin; }	
		std::vector<double> getLocalBoxMax() { return _localBoxMax; }	

		int convertIdToIndex(const std::vector<int> &domainIndex);


	private:
		// Global data
		std::vector<double> _globalBoxMin;
		std::vector<double> _globalBoxMax;

		// Domain specific data
		int _domainIndex;
  	std::vector<double> _localBoxMin;
  	std::vector<double> _localBoxMax;

		void initializeDecomposition();
		void initializeMPICommunicator();
		void initializeLocalDomain();
		void initializeGlobalBox(const double* globalBoxMin, const double* globalBoxMax);
		void initializeLocalBox();
		void initializeNeighbourIds();

		void updateLocalBox();
};
