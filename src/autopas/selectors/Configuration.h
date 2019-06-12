/**
 * @file Configuration.h
 * @author F. Gratl
 * @date 2/1/19
 */

#pragma once

#include <tuple>
#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
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
   * @param _dataLayout
   * @param _newton3
   */
  Configuration(ContainerOption _container, TraversalOption _traversal, DataLayoutOption _dataLayout,
                Newton3Option _newton3)
      : _container(_container), _traversal(_traversal), _dataLayout(_dataLayout), _newton3(_newton3) {}

  /**
   * Constructor taking no arguments. Initializes all properties to an invalid choice or false.
   */
  Configuration()
      : _container(ContainerOption(-1)),
        _traversal(TraversalOption(-1)),
        _dataLayout(DataLayoutOption(-1)),
        _newton3(Newton3Option(-1)) {}

  /**
   * Returns string representation in JSON style of the configuration object.
   * @return String representation.
   */
  std::string toString() const {
    return "{Container: " + utils::StringUtils::to_string(_container) +
           " , Traversal: " + utils::StringUtils::to_string(_traversal) +
           " , Data Layout: " + utils::StringUtils::to_string(_dataLayout) +
           " , Newton 3: " + utils::StringUtils::to_string(_newton3) + "}";
  }

  /**
   * Container option.
   */
  ContainerOption _container;
  /**
   * Traversal option.
   */
  TraversalOption _traversal;
  /**
   * Data Layout option.
   */
  DataLayoutOption _dataLayout;
  /**
   * Newton 3 option.
   */
  Newton3Option _newton3;
};

/**
 * Stream insertion operator.
 * @param os
 * @param configuration
 * @return
 */
inline std::ostream& operator<<(std::ostream& os, const Configuration& configuration) {
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
inline bool operator==(const Configuration& lhs, const Configuration& rhs) {
  return lhs._container == rhs._container and lhs._traversal == rhs._traversal and
         lhs._dataLayout == rhs._dataLayout and lhs._newton3 == rhs._newton3;
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
inline bool operator!=(const Configuration& lhs, const Configuration& rhs) { return not(lhs == rhs); }

/**
 * Comparison operator for Configuration objects. This is mainly used for configurations to have a sane ordering in e.g.
 * sets.
 *
 * Configurations are compared member wise in the order: _container, _traversal, _dataLayout, _newton3.
 *
 * @param lhs
 * @param rhs
 * @return
 */
inline bool operator<(const Configuration& lhs, const Configuration& rhs) {
  return std::tie(lhs._container, lhs._traversal, lhs._dataLayout, lhs._newton3) <
         std::tie(rhs._container, rhs._traversal, rhs._dataLayout, rhs._newton3);
}

/**
 * Hash function for Configuration objects to be used in e.g. unordered maps.
 * Aims to place integer representations of members in one large number s.th. they never overlap.
 */
struct ConfigHash {
  /**
   * Hash Function operator
   * @param configuration
   * @return
   */
  std::size_t operator()(Configuration configuration) const {
    return static_cast<std::size_t>(configuration._newton3) + static_cast<std::size_t>(configuration._dataLayout) * 10 +
           static_cast<std::size_t>(configuration._traversal) * 100 +
           static_cast<std::size_t>(configuration._container) * 10000;
  }
};

}  // namespace autopas
