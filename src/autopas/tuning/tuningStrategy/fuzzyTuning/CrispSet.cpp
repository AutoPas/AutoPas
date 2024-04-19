/**
 * @file CrispSet.cpp
 * @author Manuel Lerchner
 * @date 17.04.24
 */

#include "CrispSet.h"

namespace autopas::fuzzy_logic {

CrispSet::CrispSet(const std::string &name, const std::pair<double, double> &range) : _dimensions({{name, range}}) {}

std::shared_ptr<CrispSet> CrispSet::operator*(CrispSet &rhs) const {
  auto cartesian = std::make_shared<CrispSet>();

  cartesian->getDimensions().insert(_dimensions.begin(), _dimensions.end());
  cartesian->getDimensions().insert(rhs.getDimensions().begin(), rhs.getDimensions().end());

  return cartesian;
}

std::map<std::string, std::pair<double, double>> &CrispSet::getDimensions() { return _dimensions; }

}  // namespace autopas::fuzzy_logic