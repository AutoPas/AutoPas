/**
 * @file DomainTools.cpp
 * @author J. Körner
 * @date 13.05.2021
 */
#include "DomainTools.h"

#include <algorithm>
#include <list>

#include "autopas/utils/ArrayMath.h"

namespace DomainTools {
bool isInsideDomain(const std::array<double, 3> &coordinates, std::array<double, 3> &boxMin,
                    std::array<double, 3> &boxMax) {
  bool isInsideLocalDomain = true;
  for (int i = 0; i < coordinates.size(); ++i) {
    if (not isInsideLocalDomain) {
      break;
    }
    isInsideLocalDomain = coordinates[i] >= boxMin[i] and coordinates[i] < boxMax[i];
  }
  return isInsideLocalDomain;
}

double getDistanceToDomain(const std::array<double, 3> &coordinates, std::array<double, 3> &boxMin,
                           std::array<double, 3> &boxMax) {
  std::array<double, 3> differences = {0, 0, 0};
  for (int i = 0; i < 3; ++i) {
    differences[i] = std::clamp(coordinates[i], boxMin[i], boxMax[i]);
  }

  return autopas::utils::ArrayMath::L2Norm(differences);
}

void generateDecomposition(unsigned int subdomainCount, std::array<bool, 3> subdivideDimension,
                           std::array<int, 3> &decomposition) {
  std::list<int> primeFactors;
  // Add 2 to prime factorization as many times as subdomainCount is dividable by 2.
  while (subdomainCount % 2 == 0) {
    primeFactors.push_back(2);
    subdomainCount = subdomainCount / 2;
  }

  // Add every uneven number smaller than subdomainCount to the prime factors
  // as long as subdomainCount is dividable by this number.
  // Uneven numbers which are not a prime number will be ignored, because subdomainCount can not be divided
  // by those numbers at the point they are being checked.
  for (int i = 3; i <= subdomainCount; i = i + 2) {
    while (subdomainCount % i == 0) {
      primeFactors.push_back(i);
      subdomainCount = subdomainCount / i;
    }
  }

  // Determine number of dimensions which have to be subdivided.
  size_t numberOfDimensionsToSubdivide = 0;
  for (auto element : subdivideDimension) {
    numberOfDimensionsToSubdivide += element;
  }

  // Reduces the primeFactors to 3 elements, one for each dimension of the domain.
  // It multiplies the smalles to factors and stores it in the second factor.
  while (primeFactors.size() > numberOfDimensionsToSubdivide) {
    primeFactors.sort();
    auto firstElement = primeFactors.front();
    primeFactors.pop_front();
    primeFactors.front() *= firstElement;
  }

  // If the prime factorization ends up having less factors than dimensions in the domain,
  // fill those dimensions with 1.
  for (int i = 0; i < 3; ++i) {
    if (not primeFactors.empty() and subdivideDimension[i]) {
      decomposition[i] = primeFactors.front();
      primeFactors.pop_front();
    } else {
      decomposition[i] = 1;
    }
  }
}
}  // namespace DomainTools
