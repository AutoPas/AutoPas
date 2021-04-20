/**
 * @file DomainDecomposition.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#include "storageClass/StorageClass.h"

template<class DecompositionStorageType, DomainIdType>
class DomainDecomposition {
	public:
		virtual void update() = 0;

	protected:
		DomainDecomposition();
		~DomainDecomposition();

		// Global data
		const unsigned int _subdomainCount;
		const unsigned int _dimensionCount;
		DecompositionStorageType _decomposition;

		// Processor specific data
		DomainIdType _processorId;
		std::vector<unsigned int> _neighbourRanks;

} 
