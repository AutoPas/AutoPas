/**
 * @file RegularGridDecomposition.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#include "DomainDecomposition.h"

class RegularGridDecomposition : protected DomainDecomposition<std::vector<int>, std::vector<int>> {
	public:
		RegularGridDecomposition(const unsigned int &subdomainCount, const unsigned int &dimensionCount, const unsigned int &domainIndex, const double* globalBoxMin, const double* golbalBoxMax);
		~RegularGridDecomposition() = default;

	 	void update() override;

	private:
		// Global data
		std::vector<double> _globalBoxMin;
		std::vector<double> _globalBoxMax;

		// Domain specific data
  	std::vector<double> _localBoxMin;
  	std::vector<double> _localBoxMax;

		void initializeDecomposition();
		void initializeMPICommunicator();
		void initializeLocalDomain();
		void initializeGlobalBox(const double* globalBoxMin, const double* globalBoxMax);
		void initializeLocalBox();
		void initializeNeighbourIndices();

		int convertIdToIndex(const std::vector<int> &domainIndex);
		void updateLocalBox();
};
