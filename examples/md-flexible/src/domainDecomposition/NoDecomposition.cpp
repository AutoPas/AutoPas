/**
 * @file KdTree.cpp
 * @author J. Körner
 * @date 25.05.2021
 */
#include "NoDecomposition.h"

#include "DomainTools.h"

NoDecomposition::NoDecomposition(int argc, char **argv, const int &dimensionCount, const std::vector<double> &globalBoxMin,
             const std::vector<double> &globalBoxMax) {
  _dimensionCount = dimensionCount;
  initializeGlobalBox(globalBoxMin, globalBoxMax);
}

NoDecomposition::~NoDecomposition() { }

void NoDecomposition::update() { }

bool NoDecomposition::isInsideLocalDomain(const std::vector<double> &coordinates) {
  return DomainTools::isInsideDomain(coordinates, _globalBoxMin, _globalBoxMax);
}

void NoDecomposition::initializeGlobalBox(const std::vector<double> &globalBoxMin, const std::vector<double> &globalBoxMax) {
  _globalBoxMin.resize(_dimensionCount);
  _globalBoxMax.resize(_dimensionCount);
  for (int i = 0; i < _dimensionCount; ++i) {
    _globalBoxMin[i] = globalBoxMin[i];
    _globalBoxMax[i] = globalBoxMax[i];
  }
}
