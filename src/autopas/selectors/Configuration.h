/**
 * @file Configuration.h
 * @author F. Gratl
 * @date 2/1/19
 */

#pragma once

#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/utils/StringUtils.h"

namespace autopas {

class Configuration {
 public:
  /**
   * Constructor
   * @param _container
   * @param _traversal
   * @param _dataLayout
   * @param _newton3
   */
  Configuration(ContainerOption _container, TraversalOption _traversal, DataLayoutOption _dataLayout, bool _newton3)
      : _container(_container), _traversal(_traversal), _dataLayout(_dataLayout), _newton3(_newton3) {}

  /**
   * Constructor taking no arguments. Initializes all properties to an invalid choice or false.
   */
  Configuration()
      : _container(ContainerOption(-1)),
        _traversal(TraversalOption(-1)),
        _dataLayout(DataLayoutOption(-1)),
        _newton3(false) {}

  /**
   * Returns string representation in JSON style of the configuration object.
   * @return String representation.
   */
  std::string toString() const {
    return "{Container : " + utils::StringUtils::to_string(_container) +
           " , Traversal : " + utils::StringUtils::to_string(_traversal) +
           " , Data Layout : " + utils::StringUtils::to_string(_dataLayout) +
           " , Newton 3 : " + (_newton3 ? "On " : "Off") + "}";
  }

  ContainerOption _container;
  TraversalOption _traversal;
  DataLayoutOption _dataLayout;
  bool _newton3;
};

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
 * Hash function for Configuration objects to be used in e.g. unordered maps.
 * Aims to place integer representations of members in one large number s.th. they never overlap.
 */
struct ConfigHash {
  std::size_t operator()(Configuration configuration) const {
    return static_cast<std::size_t>(configuration._newton3) + static_cast<std::size_t>(configuration._dataLayout) * 10 +
           static_cast<std::size_t>(configuration._traversal) * 100 +
           static_cast<std::size_t>(configuration._container) * 10000;
  }
};

}  // namespace autopas