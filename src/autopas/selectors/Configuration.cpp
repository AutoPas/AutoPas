/**
 * @file Configuration.h
 * @author F. Gratl
 * @date 17. Nov 2022
 */

#include "Configuration.h"

#include "autopas/utils/StringUtils.h"

std::string autopas::Configuration::toString() const {
  return "{Container: " + container.to_string() + " , CellSizeFactor: " + std::to_string(cellSizeFactor) +
         " , Traversal: " + traversal.to_string() + " , Load Estimator: " + loadEstimator.to_string() +
         " , Data Layout: " + dataLayout.to_string() + " , Newton 3: " + newton3.to_string() + "}";
}

std::string autopas::Configuration::getCSVHeader() const { return getCSVRepresentation(true); }

std::string autopas::Configuration::getCSVLine() const { return getCSVRepresentation(false); }

bool autopas::Configuration::hasValidValues() const {
  return container != ContainerOption() and cellSizeFactor != -1 and traversal != TraversalOption() and
         loadEstimator != LoadEstimatorOption() and dataLayout != DataLayoutOption() and newton3 != Newton3Option();
}

std::string autopas::Configuration::getCSVRepresentation(bool returnHeaderOnly) const {
  auto rgx = returnHeaderOnly ?
                              // match any sequence before a colon and drop any spaces, comma or brackets before it
                 std::regex("[{, ]+([^:]+):[^,]*")
                              :
                              // match any sequence after a colon and drop any spaces, comma or brackets around it
                 std::regex(": ([^,]+)(?: ,|})");
  auto searchString = toString();
  std::sregex_iterator matchIter(searchString.begin(), searchString.end(), rgx);
  std::sregex_iterator end;
  std::stringstream retStream;

  while (matchIter != end) {
    // first submatch is the match of the capture group
    retStream << matchIter->str(1) << ",";
    ++matchIter;
  }
  auto retString = retStream.str();
  // drop trailing ','
  retString.pop_back();
  return retString;
}

std::ostream &autopas::operator<<(std::ostream &os, const autopas::Configuration &configuration) {
  return os << configuration.toString();
}

bool autopas::operator==(const autopas::Configuration &lhs, const autopas::Configuration &rhs) {
  return lhs.container == rhs.container and lhs.cellSizeFactor == rhs.cellSizeFactor and
         lhs.traversal == rhs.traversal and lhs.loadEstimator == rhs.loadEstimator and
         lhs.dataLayout == rhs.dataLayout and lhs.newton3 == rhs.newton3;
}

bool autopas::operator!=(const autopas::Configuration &lhs, const autopas::Configuration &rhs) {
  return not(lhs == rhs);
}

bool autopas::operator<(const autopas::Configuration &lhs, const autopas::Configuration &rhs) {
  return std::tie(lhs.container, lhs.cellSizeFactor, lhs.traversal, lhs.loadEstimator, lhs.dataLayout, lhs.newton3) <
         std::tie(rhs.container, rhs.cellSizeFactor, rhs.traversal, rhs.loadEstimator, rhs.dataLayout, rhs.newton3);
}
