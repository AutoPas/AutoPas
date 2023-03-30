/**
 * @file Configuration.h
 * @author F. Gratl
 * @date 2/1/19
 */

#pragma once

#include <tuple>

#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/utils/StringUtils.h"

namespace autopas {

/**
 * Class containing multiple options that form an algorithm configuration for the pairwise iteration.
 */
class Configuration {
 public:
  /**
   * Constructor
   * @param _container
   * @param _traversal
   * @param _loadEstimator
   * @param _dataLayout
   * @param _newton3
   * @param _cellSizeFactor
   * @param _verletRebuildFrequency
   */
  constexpr Configuration(ContainerOption _container, double _cellSizeFactor, TraversalOption _traversal,
                          LoadEstimatorOption _loadEstimator, DataLayoutOption _dataLayout, Newton3Option _newton3, int _verletRebuildFrequency)
      : container(_container),
        traversal(_traversal),
        loadEstimator(_loadEstimator),
        dataLayout(_dataLayout),
        newton3(_newton3),
        cellSizeFactor(_cellSizeFactor),
        verletRebuildFrequency(_verletRebuildFrequency) {}

  /**
   * Constructor taking no arguments. Initializes all properties to an invalid choice or false.
   */
  constexpr Configuration() : container(), traversal(), loadEstimator(), dataLayout(), newton3(), cellSizeFactor(-1.), verletRebuildFrequency(-1) {}




  /**
   * Returns string representation in JSON style of the configuration object.
   * @return String representation.
   */
  [[nodiscard]] std::string toString() const {
    return "{Container: " + container.to_string() + " , CellSizeFactor: " + std::to_string(cellSizeFactor) +
           " , Traversal: " + traversal.to_string() + " , Load Estimator: " + loadEstimator.to_string() +
           " , Data Layout: " + dataLayout.to_string() + " , Newton 3: " + newton3.to_string() +
           " , VerletRebuildFrequency: " + std::to_string(verletRebuildFrequency) + "}";
  }

  /**
   * Generate a csv header containing all keys from the toString() method.
   * @return Contains the header.
   */
  [[nodiscard]] std::string getCSVHeader() const { return getCSVRepresentation(true); }

  /**
   * Generate a csv representation containing all values from the toString() method.
   * @return String representing the current configuration.
   */
  [[nodiscard]] std::string getCSVLine() const { return getCSVRepresentation(false); }

  /**
   * Returns whether the configuration has been initialized with valid values or as an invalid one.
   * Does not return false if it has valid values whose combination is invalid (e.g. when the container and traversal do
   * not fit).
   * @return
   */
  [[nodiscard]] bool hasValidValues() const {
    return container != ContainerOption() and cellSizeFactor != -1 and traversal != TraversalOption() and
           loadEstimator != LoadEstimatorOption() and dataLayout != DataLayoutOption() and newton3 != Newton3Option() and
           verletRebuildFrequency != -1;
  }

  /**
   * Container option.
   */
  ContainerOption container;
  /**
   * Traversal option.
   */
  TraversalOption traversal;
  /**
   * Load Estimator option.
   */
  LoadEstimatorOption loadEstimator;
  /**
   * Data Layout option.
   */
  DataLayoutOption dataLayout;
  /**
   * Newton 3 option.
   */
  Newton3Option newton3;
  /**
   * CellSizeFactor
   */
  double cellSizeFactor;
  /**
   * VerletRebuildFrequency
   */
  int verletRebuildFrequency;

 private:
  /**
   * Helper function to return a csv representation of the current object.
   * @param returnHeaderOnly Switch to return the header or content.
   * @return
   */
  [[nodiscard]] std::string getCSVRepresentation(bool returnHeaderOnly) const {
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

};

/**
 * Stream insertion operator.
 * @param os
 * @param configuration
 * @return
 */
inline std::ostream &operator<<(std::ostream &os, const Configuration &configuration) {
  return os << configuration.toString();
}

/**
 * Equals operator for Configuration objects.
 *
 * @remark removing "inline" here leads to multiple definition errors.
 *
 * @param lhs
 * @param rhs
 * @return true iff all members are equal.
 */
inline bool operator==(const Configuration &lhs, const Configuration &rhs) {
  return lhs.container == rhs.container and lhs.cellSizeFactor == rhs.cellSizeFactor and
         lhs.traversal == rhs.traversal and lhs.loadEstimator == rhs.loadEstimator and
         lhs.dataLayout == rhs.dataLayout and lhs.newton3 == rhs.newton3 and
         lhs.verletRebuildFrequency == rhs.verletRebuildFrequency;
}

/**
 * Not-Equals operator for Configuration objects.
 *
 * @remark removing "inline" here leads to multiple definition errors.
 *
 * @param lhs
 * @param rhs
 * @return true iff at least one member is different.
 */
inline bool operator!=(const Configuration &lhs, const Configuration &rhs) { return not(lhs == rhs); }

/**
 * Comparison operator for Configuration objects. This is mainly used for configurations to have a sane ordering in e.g.
 * sets.
 *
 * Configurations are compared member wise in the order: container, cellSizeFactor, traversal, loadEstimator,
 * dataLayout, newton3.
 *
 * @param lhs
 * @param rhs
 * @return
 */
inline bool operator<(const Configuration &lhs, const Configuration &rhs) {
  return std::tie(lhs.container, lhs.cellSizeFactor, lhs.traversal, lhs.loadEstimator, lhs.dataLayout, lhs.newton3, lhs.verletRebuildFrequency) <
         std::tie(rhs.container, rhs.cellSizeFactor, rhs.traversal, rhs.loadEstimator, rhs.dataLayout, rhs.newton3, rhs.verletRebuildFrequency);
}

/**
 * Hash function for Configuration objects to be used in e.g. unordered maps.
 */
struct ConfigHash {
  /**
   * Hash Function operator
   * @param configuration
   * @return
   */
  std::size_t operator()(Configuration configuration) const {
    std::size_t enumHash = static_cast<std::size_t>(configuration.newton3) +
                           static_cast<std::size_t>(configuration.dataLayout) * 10 +
                           static_cast<std::size_t>(configuration.loadEstimator) * 100 +
                           static_cast<std::size_t>(configuration.traversal) * 1000 +
                           static_cast<std::size_t>(configuration.container) * 10000;
    std::size_t doubleHash = std::hash<double>{}(configuration.cellSizeFactor);
    std::size_t intHash = std::hash<int>{}(configuration.verletRebuildFrequency);

    return enumHash ^ doubleHash ^ intHash;
  }
};



}  // namespace autopas
