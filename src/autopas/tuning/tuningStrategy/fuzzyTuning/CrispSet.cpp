/**
 * @file CrispSet.cpp
 * @author Manuel Lerchner
 * @date 17.04.24
 */

#include "CrispSet.h"

namespace autopas::fuzzy_logic {

CrispSet::CrispSet(const std::map<std::string, std::pair<double, double>> &dimensions) : _dimensions(dimensions) {}

CrispSet CrispSet::operator*(const CrispSet &rhs) const {
  std::map<std::string, std::pair<double, double>> newDimensions;

  newDimensions.insert(_dimensions.begin(), _dimensions.end());
  newDimensions.insert(rhs._dimensions.begin(), rhs._dimensions.end());
  return CrispSet(newDimensions);
}

const std::map<std::string, std::pair<double, double>> &CrispSet::getDimensions() const { return _dimensions; }

}  // namespace autopas::fuzzy_logic