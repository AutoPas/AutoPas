/**
 * @file DomainTools.cpp
 * @author J. Körner
 * @date 13.05.2021
 */
#include "DomainTools.h"

namespace DomainTools {
bool isInsideDomain(const std::vector<double> &coordinates, std::vector<double> &boxMin, std::vector<double> &boxMax) {
  bool isInsideLocalDomain = true;
  for (int i = 0; i < coordinates.size(); ++i) {
    if (!isInsideLocalDomain) {
      break;
    }
    isInsideLocalDomain = coordinates[i] > boxMin[i] && coordinates[i] < boxMax[i];
  }
  return isInsideLocalDomain;
}
}
