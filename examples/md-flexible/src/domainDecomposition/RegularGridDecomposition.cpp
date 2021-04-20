/**
 * @file RegularGridDecomposition.cpp
 * @author J. Körner
 * @date 19.04.2021
 */
#include "RegularGridDecomposition.h"

#include <algorithm>
#include <list>
#include <math>

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

	void generateRegularGridDecomposition(const unsigned int subdomainCount, const unsigned int dimensionCount, std::vector<unsigned int>& oGridDimensions){
		std::list<unsigned int> primeFactors;
		calculatePrimeFactors(subdomainCount, primeFactors);

		while (primeFactors.size() > dimensionCount)
		{
			primeFactors.sort();
			auto firstElement = primeFactors.front();
			primeFactors.pop_front();
			primeFactors.front() *= firstElement; 
		}

		oGridDimensions.resize(subdomainCount);

		for (auto& dimensionSize : oGridDimensions)
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

}

RegularGridDecomposition::RegularGridDecomposition(const unsigned int subdomainCount, const unsigned int dimensionCount) : _subdomainCount(subdomainCount), _dimensionCount(dimensionCount) {
	_neighbourRanks.resize(_dimensionCount * 2);

	update();
}

void RegularGridDecomposition::update(){
	generateRegularGridDecomposition(_subdomainCount, _dimensionCount, _decomposition);

	std::vector<int> periods(_dimensionCount, 0);
	MPI_Cart_get(MPI_COMM_WORLD, _dimensionCount, _decomposition, periods, _processorId);

	computeNeighbourRanks();
}

unsigned int RegularGridDecomposition::convertIdToRank(std::vector<unsigned int> processorId){
	unsigned int neighbourRank = 0;
	for (auto &dimension : _dimensions){
		//@todo: calaculate neighbour rank
	}
}

void RegularGridDecomposition::computeNeighbourRanks(){
	for (unsigned int i = 0; i < _dimensionCount; ++i){
		auto neighbourIndex = i * 2;
		std::vector<unsigned int> preceedingNeighbourId = _processorId;
		preceedingNeighbourId[i] = --preceedingNeighbourId[i]%_decomposition[i];
		_neighbourRanks[neighbourIndex] = convertIdToRank(preceedingNeighbourId);

		++neighbourIndex;
		std::vector<unsigned int> succeedingNeighbourId = _processorId;
		succeedingNeighbourId[i] = ++succeedingNeighbourId[i]%_decomposition[i];
		_neighbourRanks[neighbourIndex] = convertIdToRank(succeedingNeighbourId);
	}
}
