/**
 * @file Vector.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#include "StorageClass.h"

template<class T>
class Tree : StorageClass<T> {
	public:
		Tree() = default;
		~Tree() = default;

		void getNeighboursOfSubdomain(DomainId subdomainIndex, std::vector<DomainId> neighbours) override;
}

