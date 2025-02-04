/**
 * @file Configuration.h
 * @author F. Gratl
 * @date 1 Feb. 2019
 */

#pragma once

#include <tuple>

#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/baseFunctors/Functor.h"
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
   * @param _interactionType
   *
   * @note needs constexpr (hence inline) constructor to be a literal.
   */
  constexpr Configuration(ContainerOption _container, double _cellSizeFactor, TraversalOption _traversal,
                          LoadEstimatorOption _loadEstimator, DataLayoutOption _dataLayout, Newton3Option _newton3,
                          InteractionTypeOption _interactionType)
      : container(_container),
        traversal(_traversal),
        loadEstimator(_loadEstimator),
        dataLayout(_dataLayout),
        newton3(_newton3),
        cellSizeFactor(_cellSizeFactor),
        interactionType(_interactionType) {}

  /**
   * Constructor taking no arguments. Initializes all properties to an invalid choice or false.
   * @note needs constexpr (hence inline) constructor to be a literal.
   */
  constexpr Configuration()
      : container(), traversal(), loadEstimator(), dataLayout(), newton3(), cellSizeFactor(-1.), interactionType() {}

  /**
   * Returns string representation in JSON style of the configuration object.
   * @return String representation.
   */
  [[nodiscard]] std::string toString() const;

  /**
   * Returns a short string representation of the configuration object, suitable for tabular output.
   * @param fixedLength See Option::to_string().
   * @return A short string representation.
   */
  [[nodiscard]] std::string toShortString(bool fixedLength = true) const {
    return "{" + interactionType.to_string(interactionType) + " , " + container.to_string(fixedLength) + " , " +
           std::to_string(cellSizeFactor) + " , " + traversal.to_string(fixedLength) + " , " +
           loadEstimator.to_string(fixedLength) + " , " + dataLayout.to_string(fixedLength) + " , " +
           newton3.to_string(fixedLength) + "}";
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
   * @warning This function checks by using CompatibleTraversals and CompatibleLoadEstimators. What is not checked by
   * this function is whether some limitations are imposed during the construction of the traversal nor if the
   * desired functor imposes some find of restriction on newton3.
   * @param silent if false, WARN-level logs are printed if a configuration is not compatible. This is desirable during
   * a simulation as it suggests something with the tuning strategy or elsewhere is wrong. Before the simulation, this
   * function is used to reduce the search space to only compatible values. Here incompatible values do not indicate
   * that something is wrong and so silent==true can be used.
   * @return True if all options are compatible to each other.
   */
  [[nodiscard]] bool hasCompatibleValues(bool silent = false) const;

  /**
   * Checks if any of the configuration values are incompatible with each other. Uses the non-functor variant of
   * hasCompatibleValues to check compatibility using CompatibleTraversals and CompatibleLoadEstimators, and then
   * checks that the configuration is compatible with the functor (i.e. the functor supports the newton3 type) and that
   * the traversal can be constructed with the chosen container and configuration (which may fail due to additional
   * constraints that are only forced in the generation of the traversal.)
   * @tparam Functor_T Functor type
   * @tparam Particle_T Particle type
   * @param functor functor intended for use with configuration
   * @param builtContainer the container with which the traversal is attempted to be generated
   * @param silent if false, WARN-level logs are printed if a configuration is not compatible.
   * @return True if the configuration is compatible.
   */
  template <class Functor_T, class Particle_T>
  [[nodiscard]] bool hasCompatibleValues(Functor_T &functor, ParticleContainerInterface<Particle_T> &builtContainer, bool silent = false) const;

  /**
   * Check if all discrete options of the given configuration are equal to this'.
   * @param rhs
   * @return
   */
  bool equalsDiscreteOptions(const Configuration &rhs) const;

  /**
   * Check if all continuous options of the given configuration are equal to this configuration.
   * @param rhs configuration compared against.
   * @param epsilon Maximal allowed absolute difference between two continuous values to be considered equal.
   * @return
   */
  bool equalsContinuousOptions(const autopas::Configuration &rhs, double epsilon = 1e-12) const;

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
    std::size_t enumHash = static_cast<std::size_t>(configuration.interactionType) +
                           static_cast<std::size_t>(configuration.newton3) * 10 +
                           static_cast<std::size_t>(configuration.dataLayout) * 100 +
                           static_cast<std::size_t>(configuration.loadEstimator) * 1000 +
                           static_cast<std::size_t>(configuration.traversal) * 10000 +
                           static_cast<std::size_t>(configuration.container) * 100000;
    std::size_t doubleHash = std::hash<double>{}(configuration.cellSizeFactor);

    return enumHash ^ doubleHash;
  }
};

}  // namespace autopas
