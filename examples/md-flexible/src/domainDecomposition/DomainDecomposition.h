/**
 * @file DomainDecomposition.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#include "mpi.h"

#include <vector>

template<class DecompositionStorageType, class DomainIdType>
class DomainDecomposition {
	public:
		virtual void update() = 0;
		const std::vector<unsigned int> getNeighbourDomainIndices() const {
			return _neighbourDomainIndices;
		};

	protected:
		DomainDecomposition() = default;
		~DomainDecomposition() = default;

		// Global data
		int _subdomainCount;
		int _dimensionCount;
		DecompositionStorageType _decomposition;

		// @todo Try to remove MPI from the decomposition logic
		MPI_Comm _communicator;

		// Domain specific data
		int _domainIndex;
		DomainIdType _domainId;
		std::vector<unsigned int> _neighbourDomainIndices;
};
