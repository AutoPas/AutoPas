/**
 * @file CartesianGridDecomposition.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#include "DomainDecomposition.h"
#include "storageClass/Vector.h"

class CartesianGridDecomposition : DomainDecomposition {
	public:
		CartesianGridDecomposition() = default;
		~CartesianGridDecomposition() = default;

		void generate()	override;
		void update() override;
		void getNeighboursOfProcessor(const unsigned int processorIndex) override;
};
