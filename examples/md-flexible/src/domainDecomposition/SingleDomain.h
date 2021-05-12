/**
 * @file SingleDomain.h
 * @author J. Körner
 * @date 06.05.2021
 */
#pragma once

#include "DomainDecomposition.h"

class SingleDomain final : DomainDecomposition {
	public:
		SingleDomain(int argc, char** argv, const int &dimensionCount, const std::vector<double> &globalBoxMin,
			const std::vector<double> &globalBoxMax);
		~SingleDomain();

		void update() override;
		const int getDimensionCount() override { return _dimensionCount; }
		std::vector<double> getLocalBoxMin() override { return _globalBoxMin; }
		std::vector<double> getLocalBoxMax() override { return _globalBoxMax; }

	private:
		int _dimensionCount;
		std::vector<double> _globalBoxMin;
		std::vector<double> _globalBoxMax;

		void initializeGlobalBox(const std::vector<double> &globalBoxMin, const std::vector<double> &globalBoxMax);
};

