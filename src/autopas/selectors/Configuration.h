/**
 * @file Configuration.h
 * @author F. Gratl
 * @date 1 Feb. 2019
 */

#pragma once

#include <tuple>

#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"

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
   *
   * @note needs constexpr (hence inline) constructor to be a literal.
   */
  constexpr Configuration(ContainerOption _container, double _cellSizeFactor, TraversalOption _traversal,
                          LoadEstimatorOption _loadEstimator, DataLayoutOption _dataLayout, Newton3Option _newton3)
      : container(_container),
        traversal(_traversal),
        loadEstimator(_loadEstimator),
        dataLayout(_dataLayout),
        newton3(_newton3),
        cellSizeFactor(_cellSizeFactor) {}

  /**
   * Constructor taking no arguments. Initializes all properties to an invalid choice or false.
   * @note needs constexpr (hence inline) constructor to be a literal.
   */
  constexpr Configuration() : container(), traversal(), loadEstimator(), dataLayout(), newton3(), cellSizeFactor(-1.) {}

  /**
   * Returns string representation in JSON style of the configuration object.
   * @return String representation.
   */
  [[nodiscard]] std::string toString() const;

  /**
   * Returns a short string representation of the configuration object, suitable for tabular output.
   * @return A short string representation.
   */
  [[nodiscard]] std::string toShortString() const {
    return "{" + container.to_string(true) + " , " + std::to_string(cellSizeFactor) + " , " +
           traversal.to_string(true) + " , " + loadEstimator.to_string(true) + " , " + dataLayout.to_string(true) +
           " , " + newton3.to_string(true) + "}";
  }

  /**
   * Generate a csv header containing all keys from the toString() method.
   * @return Contains the header.
   */
  [[nodiscard]] std::string getCSVHeader() const;

  /**
   * Generate a csv representation containing all values from the toString() method.
   * @return String representing the current configuration.
   */
  [[nodiscard]] std::string getCSVLine() const;

  /**
   * Returns whether the configuration has been initialized with valid values or as an invalid one.
   * Does not return false if it has valid values whose combination is invalid (e.g. when the container and traversal do
   * not fit).
   * @return
   */
  [[nodiscard]] bool hasValidValues() const;

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

 private:
  /**
   * Helper function to return a csv representation of the current object.
   * @param returnHeaderOnly Switch to return the header or content.
   * @return
   */
  [[nodiscard]] std::string getCSVRepresentation(bool returnHeaderOnly) const;
};

/**
 * Stream insertion operator.
 * @param os
 * @param configuration
 * @return
 */
std::ostream &operator<<(std::ostream &os, const Configuration &configuration);

/**
 * Stream extraction operator.
 * @param in
 * @param configuration
 * @return
 */
inline std::istream &operator>>(std::istream &in, Configuration &configuration) {
  constexpr auto max = std::numeric_limits<std::streamsize>::max();
  in.ignore(max, ':');
  in >> configuration.container;
  in.ignore(max, ':');
  in >> configuration.cellSizeFactor;
  in.ignore(max, ':');
  in >> configuration.traversal;
  in.ignore(max, ':');
  in >> configuration.loadEstimator;
  in.ignore(max, ':');
  in >> configuration.dataLayout;
  in.ignore(max, ':');
  in >> configuration.newton3;
  return in;
}

/**
 * Equals operator for Configuration objects.
 * @param lhs
 * @param rhs
 * @return true iff all members are equal.
 */
bool operator==(const Configuration &lhs, const Configuration &rhs);

/**
 * Not-Equals operator for Configuration objects.
 * @param lhs
 * @param rhs
 * @return true iff at least one member is different.
 */
bool operator!=(const Configuration &lhs, const Configuration &rhs);

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
bool operator<(const Configuration &lhs, const Configuration &rhs);

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

    return enumHash ^ doubleHash;
  }
};

}  // namespace autopas
