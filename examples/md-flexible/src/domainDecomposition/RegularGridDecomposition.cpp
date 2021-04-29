/**
 * @file RegularGridDecomposition.cpp
 * @author J. Körner
 * @date 19.04.2021
 */
#include "RegularGridDecomposition.h"

#include <algorithm>
#include <functional>
#include <list>
#include <math.h>
#include <numeric>

namespace {
	void calculatePrimeFactors(unsigned int number, std::list<unsigned int>& oPrimeFactors){
		while (number%2 == 0)
		{
			oPrimeFactors.push_back(2);
			number = number / 2;
		}

		for (unsigned int i = 3; i <= number; i = i+2)
		{
			while (number%i == 0)
			{
				oPrimeFactors.push_back(i);
				number = number / i;
			}
		}
	}

}

RegularGridDecomposition::RegularGridDecomposition(const unsigned int &subdomainCount, const unsigned int &dimensionCount, const unsigned int &domainIndex, const double* globalBoxMin, const double* globalBoxMax){
	_subdomainCount = subdomainCount;
	_dimensionCount = dimensionCount;

	initializeDecomposition();

	initializeMPICommunicator();

	initializeLocalDomain();

	initializeGlobalBox(globalBoxMin, globalBoxMax);

	initializeLocalBox();

	initializeNeighbourIndices();
}

void RegularGridDecomposition::update(){
	updateLocalBox();
}

void RegularGridDecomposition::initializeDecomposition(){
	std::list<unsigned int> primeFactors;
	calculatePrimeFactors(_subdomainCount, primeFactors);

	while (primeFactors.size() > _dimensionCount)
	{
		primeFactors.sort();
		auto firstElement = primeFactors.front();
		primeFactors.pop_front();
		primeFactors.front() *= firstElement; 
	}

	_decomposition.resize(_dimensionCount);

	for (auto& dimensionSize : _decomposition)
	{
		if (primeFactors.size() > 0) {
			dimensionSize = primeFactors.front();
			primeFactors.pop_front();
		}
		else {
			dimensionSize = 1;
		}
	}
}

void RegularGridDecomposition::initializeMPICommunicator(){
	std::vector<int> periods(_dimensionCount, 1);
  MPI_Cart_create(MPI_COMM_WORLD, _dimensionCount, _decomposition.data(), periods.data(), true, &_communicator);
}

void RegularGridDecomposition::initializeLocalDomain(){
	_domainId.resize(_dimensionCount);
  MPI_Comm_rank(_communicator, &_domainIndex);

	std::vector<int> periods(_dimensionCount, 1);
	MPI_Cart_get(_communicator, _dimensionCount, _decomposition.data(), periods.data(), _domainId.data());
}


void RegularGridDecomposition::initializeNeighbourIndices(){
	_neighbourDomainIndices.resize(_dimensionCount * 2);

	for (unsigned int i = 0; i < _dimensionCount; ++i){
		auto neighbourIndex = i * 2;
		std::vector<int> preceedingNeighbourId = _domainId;
		preceedingNeighbourId[i] =
			(--preceedingNeighbourId[i] + _decomposition[i])%_decomposition[i];
		_neighbourDomainIndices[neighbourIndex] = convertIdToIndex(preceedingNeighbourId);

		++neighbourIndex;
		std::vector<int> succeedingNeighbourId = _domainId;
		succeedingNeighbourId[i] =
			(++succeedingNeighbourId[i] + _decomposition[i])%_decomposition[i];
		_neighbourDomainIndices[neighbourIndex] = convertIdToIndex(succeedingNeighbourId);
	}
}

void RegularGridDecomposition::initializeGlobalBox(const double* globalBoxMin, const double* globalBoxMax){
	_globalBoxMin.resize(_dimensionCount);
	_globalBoxMax.resize(_dimensionCount);
	for (int i = 0; i < _dimensionCount; ++i){
		_globalBoxMin[i] = globalBoxMin[i];
		_globalBoxMax[i] = globalBoxMax[i];
	}
}

void RegularGridDecomposition::initializeLocalBox(){
	_localBoxMin.resize(_dimensionCount);
	_localBoxMax.resize(_dimensionCount);
	updateLocalBox();
}

int RegularGridDecomposition::convertIdToIndex(const std::vector<int> &domainId){
	int neighbourDomainIndex = 0;

	for (unsigned int i = 0; i < _dimensionCount; ++i){
		int accumulatedTail = 1;

		if (i < _decomposition.size()-1) {
			accumulatedTail = std::accumulate(_decomposition.begin()+i+1, _decomposition.end(),
				1, std::multiplies<int>());
		}

		neighbourDomainIndex += accumulatedTail * domainId[i];
	}

	return neighbourDomainIndex;
}

void RegularGridDecomposition::updateLocalBox(){
  for (int i = 0; i < _dimensionCount; ++i) {
		double _globalBoxWidth = _globalBoxMax[i] - _globalBoxMin[i];
		double localBoxWidth = _globalBoxWidth / _decomposition[i];

    _localBoxMin[i] = _domainId[i] * localBoxWidth + _globalBoxMin[i];
    _localBoxMax[i] = (_domainId[i] + 1) * localBoxWidth + _globalBoxMin[i];

    if (_domainId[i] == 0) {
      _localBoxMin[i] = _globalBoxMin[i];
    } else if (_domainId[i] == _decomposition[i] - 1) {
      _localBoxMax[i] = _globalBoxMax[i];
    }
  }
}
