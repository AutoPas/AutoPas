/**
 * @file DomainTools.cpp
 * @author J. Körner
 * @date 13.05.2021
 */
#include "DomainTools.h"

#include <list>

namespace DomainTools {
bool isInsideDomain(const std::vector<double> &coordinates, std::vector<double> &boxMin, std::vector<double> &boxMax) {
  bool isInsideLocalDomain = true;
  for (int i = 0; i < coordinates.size(); ++i) {
    if (!isInsideLocalDomain) {
      break;
    }
    isInsideLocalDomain = coordinates[i] >= boxMin[i] && coordinates[i] < boxMax[i];
  }
  return isInsideLocalDomain;
}

void generateDecomposition(unsigned int subdomainCount, int dimensionCount, std::vector<int> &oDecomposition) {
  std::list<int> primeFactors;
  while (subdomainCount % 2 == 0) {
    primeFactors.push_back(2);
    subdomainCount = subdomainCount / 2;
  }

  for (unsigned int i = 3; i <= subdomainCount; i = i + 2) {
    while (subdomainCount % i == 0) {
      primeFactors.push_back(i);
      subdomainCount = subdomainCount / i;
    }
  }

  while (primeFactors.size() > dimensionCount) {
    primeFactors.sort();
    auto firstElement = primeFactors.front();
    primeFactors.pop_front();
    primeFactors.front() *= firstElement;
  }

  oDecomposition.resize(dimensionCount);

  for (auto &dimensionSize : oDecomposition) {
    if (primeFactors.size() > 0) {
      dimensionSize = primeFactors.front();
      primeFactors.pop_front();
    } else {
      dimensionSize = 1;
    }
  }
}
}
