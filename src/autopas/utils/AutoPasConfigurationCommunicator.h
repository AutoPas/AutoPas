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
 * Distributes the provided configurations globally for equal work loads.
 * All parameters' values (except for comm) are only relevant at the root node (0).
 * All parameters' values (except for comm) will be changed by this function
 * @param containerOptions
 * @param cellSizeFactors
 * @param traversalOptions
 * @param datalayouts
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
Configuration deserializeConfiguration(SerializedConfiguration config) {
  double cellSizeFactor;
  std::memcpy(&cellSizeFactor, &config[4], sizeof(double));
  return Configuration(static_cast<ContainerOption::Value>(config[0]), cellSizeFactor,
                       static_cast<TraversalOption::Value>(config[1]), static_cast<DataLayoutOption::Value>(config[2]),
                       static_cast<Newton3Option::Value>(config[3]));
}

/**
 * Handles communication to select the globally best configuration.
 * @param comm: The communicator used for sending and receiving the optimal configuration
 * @param localOptimalConfig: The locally optimal configuration to be compared with others
 * @param localOptimalTime: The time measured for localOptimalConfig
 * @return The globally optimal configuration
 */
Configuration optimizeConfiguration(AutoPas_MPI_Comm comm, Configuration localOptimalConfig, size_t localOptimalTime);

}  // namespace autopas::AutoPasConfigurationCommunicator
