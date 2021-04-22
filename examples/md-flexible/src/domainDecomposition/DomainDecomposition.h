/**
 * @file DomainDecomposition.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#include <vector>

template<class DecompositionStorageType, class DomainIdType>
class DomainDecomposition {
	public:
		virtual void update() = 0;

	protected:
		DomainDecomposition();
		~DomainDecomposition() = default;

		// Global data
		unsigned int _subdomainCount;
		unsigned int _dimensionCount;
		DecompositionStorageType _decomposition;

		// Processor specific data
		DomainIdType _processorId;
		std::vector<unsigned int> _neighbourRanks;
};
