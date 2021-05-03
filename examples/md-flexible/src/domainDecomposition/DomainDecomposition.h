/**
 * @file DomainDecomposition.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#include "mpi.h"

#include <vector>
#include <memory>

template<class DecompositionStorageType, class DomainIdType>
class DomainDecomposition {
	public:
		virtual void update() = 0;

		virtual const int getDimensionCount() = 0;
		virtual const MPI_Comm getCommunicator() = 0;

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
		
		// @todo: location of neighbour ideally is implicitly defined by index in _neighbourDOmainIndices. Therefore, this member may be redundant.
		std::vector<DomainIdType> _neighbourDomainIds;
};
