/**
 * @file Vector.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#include "StorageClass.h"

#include <vector>

class Vector : StorageClass<unsigned int> {
	public:
		Vector() = default;	
		~Vector() = default;	

		void getNeighboursOfSubdomain(DomainId& subdomainIndex, std::vector<DomainId>& neighbours) override;
}
