/**
 * @file MDFlexMPI.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#include "CartesianGridDecomposition.h"

class DomainDecompositionFactory {
	public: 
		DomainDecompositionFactory() = default;
		~DomainDecompositionFactory() = default;

		CartesianGridDecomposition GetCartesianGridDecomposition(const unsigned int& numberOfAvailableProcessors, const unsigned int& numberOfGridDimensions);
};

