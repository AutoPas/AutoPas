/**
 * @file DomainTools.cpp
 * @author J. Körner
 * @date 13.05.2021
 */
#include "DomainTools.h"

#include <algorithm>
#include <cmath>
#include <list>

#include "autopas/utils/ArrayMath.h"

namespace DomainTools {
bool isInsideDomain(const std::array<double, 3> &coordinates, std::array<double, 3> &boxMin,
                    std::array<double, 3> &boxMax) {
  bool isInsideLocalDomain = true;
  for (int i = 0; i < coordinates.size(); ++i) {
    if (!isInsideLocalDomain) {
      break;
    }
    isInsideLocalDomain = coordinates[i] >= boxMin[i] && coordinates[i] < boxMax[i];
  }
  return isInsideLocalDomain;
}

double getDistanceToDomain(const std::array<double, 3> &coordinates, std::array<double, 3> &boxMin,
                           std::array<double, 3> &boxMax) {
  if (coordinates.size() == boxMin.size() && coordinates.size() == boxMax.size()) {
    std::array<double, 3> differences = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {
      differences[i] = std::clamp(coordinates[i], boxMin[i], boxMax[i]);
    }

    return autopas::utils::ArrayMath::L2Norm(differences);
  }
  return -1;
}

void generateDecomposition(unsigned int subdomainCount, std::array<int, 3> &decomposition) {
  std::list<int> primeFactors;
  while (subdomainCount % 2 == 0) {
    primeFactors.push_back(2);
    subdomainCount = subdomainCount / 2;
  }

  for (int i = 3; i <= subdomainCount; i = i + 2) {
    while (subdomainCount % i == 0) {
      primeFactors.push_back(i);
      subdomainCount = subdomainCount / i;
    }
  }

  while (primeFactors.size() > 3) {
    primeFactors.sort();
    auto firstElement = primeFactors.front();
    primeFactors.pop_front();
    primeFactors.front() *= firstElement;
  }

  for (auto &dimensionSize : decomposition) {
    if (not primeFactors.empty()) {
      dimensionSize = primeFactors.front();
      primeFactors.pop_front();
    } else {
      dimensionSize = 1;
    }
  }
}
}  // namespace DomainTools
