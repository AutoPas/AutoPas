/**
 * @file DomainDecomposition.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#include <vector>

class DomainDecomposition {
	public:
		virtual void update() = 0;
		virtual const int getDimensionCount() = 0;
		virtual std::vector<double> getGlobalBoxMin() = 0;
		virtual std::vector<double> getGlobalBoxMax() = 0;
		virtual std::vector<double> getLocalBoxMin() = 0;
		virtual std::vector<double> getLocalBoxMax() = 0;
		virtual bool isInsideLocalDomain(std::vector<double> coordinates) = 0;
};
