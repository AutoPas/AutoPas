/**
 * @file Configuration.h
 * @author F. Gratl
 * @date 1 Feb. 2019
 */

#pragma once

#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>

#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/CompatibleVectorizationPattern.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/OpenMPKindOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/options/VectorizationPatternOption.h"
#include "autopas/utils/HashCombine.h"

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
   * @param _ompKind
   * @param _ompChunkSize
   * @param _interactionType
   * @param _vecPattern
   *
   * @note needs constexpr (hence inline) constructor to be a literal.
   */
  constexpr Configuration(ContainerOption _container, double _cellSizeFactor, TraversalOption _traversal,
                          LoadEstimatorOption _loadEstimator, DataLayoutOption _dataLayout, Newton3Option _newton3,
                          OpenMPKindOption _ompKind, size_t _ompChunkSize, InteractionTypeOption _interactionType,
                          VectorizationPatternOption _vecPattern = VectorizationPatternOption::p1xVec)
      : container(_container),
        traversal(_traversal),
        vecPattern(_vecPattern),
        loadEstimator(_loadEstimator),
        dataLayout(_dataLayout),
        newton3(_newton3),
        cellSizeFactor(_cellSizeFactor),
        ompKind(_ompKind),
        ompChunkSize(_ompChunkSize),
        interactionType(_interactionType) {}

  /**
   * Constructor taking no arguments. Initializes all properties to an invalid choice or false.
   * @note needs constexpr (hence inline) constructor to be a literal.
   */
  constexpr Configuration()
      : container(),
        traversal(),
        loadEstimator(),
        dataLayout(),
        newton3(),
        cellSizeFactor(-1.),
        ompKind(),
        ompChunkSize(0),
        interactionType() {}

  /**
   * Returns string representation in JSON style of the configuration object.
   * @return String representation.
   */
  [[nodiscard]] std::string toString() const;

  /**
   * Returns a short string representation of the configuration object, suitable for tabular output or test name.
   * @param fixedLength See Option::to_string().
   * @param forParameterizedTestName if true, creates a string representation that is safe for use as a test name.
   * @return Short string representation.
   */
  [[nodiscard]] std::string toShortString(bool fixedLength = true, bool forParameterizedTestName = false) const {
    const std::string delimiter = forParameterizedTestName ? "_" : " , ";
    auto result = (forParameterizedTestName ? "" : "{") + interactionType.to_string() + delimiter +
                  container.to_string(fixedLength) + delimiter + std::to_string(cellSizeFactor) + delimiter +
                  traversal.to_string(fixedLength) + delimiter + loadEstimator.to_string(fixedLength) + delimiter +
                  dataLayout.to_string(fixedLength) + delimiter + newton3.to_string(fixedLength) + delimiter +
                  ompKind.to_string(fixedLength) + delimiter + std::to_string(ompChunkSize) + delimiter +
                  vecPattern.to_string(fixedLength) + (forParameterizedTestName ? "" : "}");

    // For parameterized test names, no punctuation is allowed except "_"
    if (forParameterizedTestName) {
      std::ranges::replace(result, '.', '_');
      std::ranges::replace(result, '-', '_');
      std::ranges::replace(result, '/', '_');
    }

    return result;
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
   * Checks if any of the configuration values are incompatible with each other.
   * @return True if all options are compatible to each other.
   */
  [[nodiscard]] bool hasCompatibleValues() const;

  /**
   * A tuple (std::tie) of references to all members, in a fixed order.
   *
   * Useful for treating the configuration as a tuple and applying it in variadic templates, or in tuple comparison
   * operators, which improves the maintainability of configuration ordering, comparison, hashing, and
   * (de)serialization, as these features do not need to be modified when new configuration components are added, only
   * this tie operator.
   *
   * @return Tuple of const references to all members.
   */
  [[nodiscard]] auto tie() const {
    return std::tie(container, cellSizeFactor, traversal, loadEstimator, dataLayout, newton3, ompKind, ompChunkSize,
                    interactionType, vecPattern);
  }

  /**
   * Non-const overload of tie()
   * @return Tuple of mutable references to all members.
   */
  [[nodiscard]] auto tie() {
    return std::tie(container, cellSizeFactor, traversal, loadEstimator, dataLayout, newton3, ompKind, ompChunkSize,
                    interactionType, vecPattern);
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
   * Vectorization Pattern option
   */
  VectorizationPatternOption vecPattern;
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
   * OpenMP (Schedule) Kind Option, e.g. static, dynamic.
   */
  OpenMPKindOption ompKind;
  /**
   * OpenMP Chunk Size.
   */
  size_t ompChunkSize;
  /**
   * Interaction type of the configuration.
   */
  InteractionTypeOption interactionType;

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
std::istream &operator>>(std::istream &in, Configuration &configuration);

/**
 * Equals operator for Configuration objects.
 * @param lhs
 * @param rhs
 * @return true iff all components are equal.
 */
bool operator==(const Configuration &lhs, const Configuration &rhs);

/**
 * Not-Equals operator for Configuration objects.
 * @param lhs
 * @param rhs
 * @return true iff at least one component is different.
 */
bool operator!=(const Configuration &lhs, const Configuration &rhs);

/**
 * Comparison operator for Configuration objects. This is mainly used for configurations to have a sane ordering in e.g.
 * sets.
 *
 * Configurations are compared member wise in the order given by Configuration::tie().
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
  std::size_t operator()(const Configuration &configuration) const {
    return std::apply([](const auto &...members) { return utils::hashCombine(members...); }, configuration.tie());
  }
};

namespace internal {
/**
 * Sum of the sizes of all members referenced by a Configuration::tie() tuple.
 * @tparam Tuple_T Type of the tie tuple.
 * @tparam Is Index sequence over the tuple elements.
 * @return Sum of sizeof over all referenced members.
 */
template <class Tuple_T, std::size_t... Is>
constexpr std::size_t sumOfMemberSizes(std::index_sequence<Is...>) {
  return (sizeof(std::remove_reference_t<std::tuple_element_t<Is, Tuple_T>>) + ...);
}
}  // namespace internal

/**
 * Type of the tuple returned by Configuration::tie().
 */
using ConfigurationTie = decltype(std::declval<const Configuration &>().tie());

/**
 * Number of bytes a Configuration occupies when every member is serialized by value.
 *
 * Derived from Configuration::tie(), so it cannot drift out of sync with the members.
 */
inline constexpr std::size_t serializedConfigurationSize =
    internal::sumOfMemberSizes<ConfigurationTie>(std::make_index_sequence<std::tuple_size_v<ConfigurationTie>>{});

}  // namespace autopas
