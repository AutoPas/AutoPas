/**
 * @file RegularGridDecomposition.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#include "DomainDecomposition.h"

class RegularGridDecomposition : protected DomainDecomposition<std::vector<int>, std::vector<int>> {
	public:
		RegularGridDecomposition(const unsigned int &subdomainCount, const unsigned int &dimensionCount);
		~RegularGridDecomposition() = default;

	 void update() override;

	private:
		unsigned int convertIdToRank(const std::vector<int> &processorId);
		void computeNeighbourRanks();
};
