/**
 * @file RegularGridDecomposition.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#include "DomainDecomposition.h"
#include "storageClass/Vector.h"

class RegularGridDecomposition : DomainDecomposition<std::vector<unsigned int>, std::vector<unsigned int>> {
	public:
		RegularGridDecomposition(const unsigned int &numberOfSubdomains, const unsigned int &numberOfDimensions);
		~RegularGridDecomposition() = default;

	private:
		void convertIdToRank(const std::vector<unsigned int> &processorId);
		void computeNeighbourRanks();
};
