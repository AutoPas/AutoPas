/**
 * @file SingleDomain.cpp
 * @author J. Körner
 * @date 06.05.2021
 */
#include "SingleDomain.h"

#include "DomainTools.h"

SingleDomain::SingleDomain(int argc, char **argv, const int &dimensionCount, const std::vector<double> &globalBoxMin,
                           const std::vector<double> &globalBoxMax) {
  _dimensionCount = dimensionCount;
  initializeGlobalBox(globalBoxMin, globalBoxMax);
}

void SingleDomain::update() {
  // Do nothing
}

void SingleDomain::initializeGlobalBox(const std::vector<double> &globalBoxMin,
                                       const std::vector<double> &globalBoxMax) {
  _globalBoxMin.resize(_dimensionCount);
  _globalBoxMax.resize(_dimensionCount);
  for (int i = 0; i < _dimensionCount; ++i) {
    _globalBoxMin[i] = globalBoxMin[i];
    _globalBoxMax[i] = globalBoxMax[i];
  }
}

bool SingleDomain::isInsideLocalDomain(const std::vector<double> &coordinates) {
  return DomainTools::isInsideDomain(coordinates, _globalBoxMin, _globalBoxMax);
}

bool SingleDomain::isInsideLocalDomain(const std::array<double, 3> &coordinates) {
  std::vector<double> coordinatesVector;
  coordinatesVector.insert(coordinatesVector.begin(), coordinates.begin(), coordinates.end());

  return isInsideLocalDomain(coordinatesVector);
}
