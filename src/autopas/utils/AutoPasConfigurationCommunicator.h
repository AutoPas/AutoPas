/**
 * @file AutoPasConfigurationCommunicator.h
 * @author W. Thieme
 * @date 30.04.2020
 */

#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "WrapMPI.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/selectors/Configuration.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/NumberSet.h"

/**
 * Provides several functions for handling configurations among mpi ranks.
 * This includes functionality for (de)serialization of configurations, splitting up search spaces based on ranks,
 * and finding the globally optimal configuration given time measurements.
 */

namespace autopas::AutoPasConfigurationCommunicator {

/** type definition for the serialization of configurations. A serialized config is an array of 12 bytes */
using SerializedConfiguration = std::array<std::byte, 12>;

/**
 * Simply a shorter way of static_casting from Option to std::byte
 * @tparam TOption
 * @param option
 * @return
 */
template <typename TOption>
inline std::byte castToByte(TOption option) {
  return static_cast<std::byte>(static_cast<typename TOption::Value>(option));
}

/**
 * Holds functionality needed to iterate through ranks and configurations simultaneously in an evenly distributed manner
 * Also provides useful information for dealing with infinite NumberSets as cellSizeFactors
 */
struct IteratorHandler {
  /**
   * Constructor for IteratorHandler
   * @param containerOptions
   * @param cellSizeFactors
   * @param traversalOptions
   * @param dataLayoutOptions
   * @param newton3Options
   * @param numConfigs
   * @param commSize
   */
  IteratorHandler(std::set<ContainerOption> &containerOptions, std::set<double> &cellSizeFactors,
                  std::set<TraversalOption> &traversalOptions, std::set<DataLayoutOption> &dataLayoutOptions,
                  std::set<Newton3Option> &newton3Options, const int numConfigs, const int commSize)
      : _containerOptions(&containerOptions),
        _cellSizeFactors(&cellSizeFactors),
        _traversalOptions(&traversalOptions),
        _dataLayoutOptions(&dataLayoutOptions),
        _newton3Options(&newton3Options),
        _containerIt(containerOptions.begin()),
        _cellSizeFactorIt(cellSizeFactors.begin()),
        _traversalIt(traversalOptions.begin()),
        _dataLayoutIt(dataLayoutOptions.begin()),
        _newton3It(newton3Options.begin()),
        _rankIterator(0),
        _remainingBlockSize(commSize >= numConfigs ? commSize / numConfigs : numConfigs / commSize),
        _remainder(commSize >= numConfigs ? commSize % numConfigs : numConfigs % commSize),
        _infiniteCellSizeFactorsOffset(0),
        _infiniteCellSizeFactorsBlockSize(commSize >= numConfigs ? commSize / numConfigs : numConfigs / commSize) {}

  /**
   * Advances the rankIterator (getRankIterator()) and/or the Option iterators for a single step such that repeated
   * execution of this function ends up in both reaching their respective ends simultaneously
   * @param numConfigs
   * @param commSize
   */
  void advanceIterators(const int numConfigs, const int commSize);

  /**
   * Alternative getter for all Configuration iterators
   * @param containerIt out
   * @param cellSizeFactorIt out
   * @param traversalIt out
   * @param dataLayoutIt out
   * @param newton3It out
   */
  inline void getConfigIterators(std::set<ContainerOption>::iterator &containerIt,
                                 std::set<double>::iterator &cellSizeFactorIt,
                                 std::set<TraversalOption>::iterator &traversalIt,
                                 std::set<DataLayoutOption>::iterator &dataLayoutIt,
                                 std::set<Newton3Option>::iterator &newton3It) {
    containerIt = _containerIt;
    cellSizeFactorIt = _cellSizeFactorIt;
    traversalIt = _traversalIt;
    dataLayoutIt = _dataLayoutIt;
    newton3It = _newton3It;
  }

  /**
   * Getter for the rankIterator. The value will correspond to the rank that holds the Options that the other iterators
   * point to
   * @return
   */
  inline int getRankIterator() const { return _rankIterator; }

  /**
   * Getter for the number of ranks smaller than getRankIterator that have the exact same configs assigned to them.
   * @return
   */
  inline int getInfiniteCellSizeFactorsOffset() const { return _infiniteCellSizeFactorsOffset; }

  /**
   * Getter for the number of ranks in total that have the exact same configs assigned to them.
   * @return
   */
  inline int getInfiniteCellSizeFactorsBlockSize() const { return _infiniteCellSizeFactorsBlockSize; }

  /**
   * Getter for containerIterator
   * @return
   */
  inline std::set<ContainerOption>::iterator getContainerIterator() const { return _containerIt; }

  /**
   * Getter for the CellSizeFactorIterator
   * @return
   */
  inline std::set<double>::iterator getCellSizeFactorIterator() const { return _cellSizeFactorIt; }

  /**
   * Getter for the TraversalIterator
   * @return
   */
  inline std::set<TraversalOption>::iterator getTraversalIterator() const { return _traversalIt; }

  /**
   * Getter for the DataLayoutIterator
   * @return
   */
  inline std::set<DataLayoutOption>::iterator getDataLayoutIterator() const { return _dataLayoutIt; }

  /**
   * Getter for the Newton3Iterator
   * @return
   */
  inline std::set<Newton3Option>::iterator getNewton3Iterator() const { return _newton3It; }

 private:
  /**
   * Selects the next config from containers X cellSizeFactors X traversals X dataLayouts X newton3Options
   * Later named options are changed first, because they are easier to switch between simulation steps
   */
  inline void advanceConfigIterators();

  const std::set<ContainerOption> * const _containerOptions;
  const std::set<double> * const _cellSizeFactors;
  const std::set<TraversalOption> * const _traversalOptions;
  const std::set<DataLayoutOption> * const _dataLayoutOptions;
  const std::set<Newton3Option> * const _newton3Options;
  std::set<ContainerOption>::iterator _containerIt;
  std::set<double>::iterator _cellSizeFactorIt;
  std::set<TraversalOption>::iterator _traversalIt;
  std::set<DataLayoutOption>::iterator _dataLayoutIt;
  std::set<Newton3Option>::iterator _newton3It;
  int _rankIterator;
  int _remainingBlockSize;
  int _remainder;
  int _infiniteCellSizeFactorsOffset;
  int _infiniteCellSizeFactorsBlockSize;
};

/**
 * Distributes the provided configurations globally for equal work loads.
 * All parameters' values (except for comm) are only relevant at the root node (0).
 * All parameters' values (except for comm) will be changed by this function
 * @param containerOptions
 * @param cellSizeFactors
 * @param traversalOptions
 * @param dataLayoutOptions
 * @param newton3Options
 * @param comm
 */
void distributeConfigurations(std::set<ContainerOption> &containerOptions, NumberSet<double> &cellSizeFactors,
                              std::set<TraversalOption> &traversalOptions,
                              std::set<DataLayoutOption> &dataLayoutOptions, std::set<Newton3Option> &newton3Options,
                              AutoPas_MPI_Comm comm);

/**
 * Serializes a configuration object for communication via MPI
 * @param configuration: the configuration to be sent
 * @return The serialization
 */
SerializedConfiguration serializeConfiguration(Configuration configuration);

/**
 * Recreates a Configuration object from the object obtained by _serializeConfiguration
 * @param config: The SerializedConfiguration objects returned by _serializeConfiguration
 * @return The deserialized Configuration object
 */
Configuration deserializeConfiguration(SerializedConfiguration config);

/**
 * Handles communication to select the globally best configuration.
 * @param comm: The communicator used for sending and receiving the optimal configuration
 * @param localOptimalConfig: The locally optimal configuration to be compared with others
 * @param localOptimalTime: The time measured for localOptimalConfig
 * @return The globally optimal configuration
 */
Configuration optimizeConfiguration(AutoPas_MPI_Comm comm, Configuration localOptimalConfig, size_t localOptimalTime);

}  // namespace autopas::AutoPasConfigurationCommunicator
