/**
 * @file CrispSet.cpp
 * @author Manuel Lerchner
 * @date 17.04.24
 */

#include "CrispSet.h"

#include <numeric>

namespace autopas::FuzzyLogic {

CrispSet::CrispSet(const std::string &name, const std::pair<double, double> &range) : _dimensions({{name, range}}) {}

std::shared_ptr<CrispSet> CrispSet::operator*(CrispSet &rhs) const {
  auto newDimensions = std::make_shared<CrispSet>();

  newDimensions->getDimensions().insert(_dimensions.begin(), _dimensions.end());
  newDimensions->getDimensions().insert(rhs.getDimensions().begin(), rhs.getDimensions().end());

  return newDimensions;
}

std::map<std::string, std::pair<double, double>> &CrispSet::getDimensions() { return _dimensions; }

CrispSet::operator std::string() const {
  return std::accumulate(_dimensions.begin(), _dimensions.end(), std::string("CrispSet: {"),
                         [](const std::string &acc, const std::pair<const std::string, std::pair<double, double>> &b) {
                           return acc + b.first + ": [" + std::to_string(b.second.first) + ", " +
                                  std::to_string(b.second.second) + "], ";
                         }) +
         "}";
}

}  // namespace autopas::FuzzyLogic