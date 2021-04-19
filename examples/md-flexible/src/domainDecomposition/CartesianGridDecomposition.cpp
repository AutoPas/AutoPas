/**
 * @file CartesianGridDecomposition.cpp
 * @author J. Körner
 * @date 19.04.2021
 */
#include "CartesianGridDecomposition.h"

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

	void generateGridDecomposition(const unsigned int subdomainCount, const unsigned int dimensionCount, std::vector<unsigned int>& oGridDimensions){
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

CartesianGridDecomposition(const unsigned int numberOfSubdomains, const unsigned int numberOfDimensions){
	 generateGridDecomposition(subdomainCount, dimensionCount, _decomposition);
}

void CartesianGridDecomposition::generate(){
		
}

void CartesianGridDecomposition::getSubdomainId(const unsigned int subdomainIndex){

}

void CartesianGridDecomposition::update(){

}

void CartesianGridDecomposition::getNeighboursOfSubdomain(/* Use processor id to identifiy neighbours */){
}


