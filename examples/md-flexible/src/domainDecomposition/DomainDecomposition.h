/**
 * @file DomainDecomposition.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#include "StorageClass.h"

class DomainDecomposition {
	public:
		virtual void generate() = 0;
		virtual void update() = 0;
		virtual void getSubdomainId(const unsigned int subdomainIndex) = 0;
		virtual void getNeighboursOfSubdomain(/* Use processor id to identifiy neighbours */) = 0;

	protected:
		DomainDecomposition();
		~DomainDecomposition();

	private:
		const unsigned int _numberOfDimensions;
		StorageClass _decomposition;
		std::vector<unsigned int> _neighbourProcessors;
} 
