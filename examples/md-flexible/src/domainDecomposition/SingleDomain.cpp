/**
 * @file SingleDomain.cpp
 * @author J. Körner
 * @date 06.05.2021
 */

#include "SingleDomain.h"

#ifdef AUTOPAS_INTERNODE_TUNING
#include "mpi.h"
#endif

SingleDomain::SingleDomain(int argc, char** argv, const int &dimensionCount,
	const std::vector<double> &globalBoxMin, const std::vector<double> &globalBoxMax) {
	#ifdef AUTOPAS_INTERNODE_TUNING
	MPI_Init(&argc, &argv); 
	#endif

	_dimensionCount = dimensionCount;
	initializeGlobalBox(globalBoxMin, globalBoxMax);
}

SingleDomain::~SingleDomain(){
#ifdef AUTOPAS_INTERNODE_TUNING
	MPI_Finalize(); 
#endif
}

void SingleDomain::update(){
	// Do nothing
}

void SingleDomain::exchangeHaloData() {
	// Do nothing
}

void SingleDomain::initializeGlobalBox(const std::vector<double> &globalBoxMin,
	const std::vector<double> &globalBoxMax){
	_globalBoxMin.resize(_dimensionCount);
	_globalBoxMax.resize(_dimensionCount);
	for (int i = 0; i < _dimensionCount; ++i){
		_globalBoxMin[i] = globalBoxMin[i];
		_globalBoxMax[i] = globalBoxMax[i];
	}
}
