/**
 * @file AutoPasConfigurationCommunicator.h
 * @author W. Thieme
 * @date 30.04.2020
 */

#pragma once

#include <array>
#include <vector>

#include "WrapMPI.h"
#include "autopas/selectors/Configuration.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/containers/CompatibleTraversals.h"

namespace autopas::AutoPasConfigurationCommunicator {

/** type definition for the serialization of configurations. A serialized config is an array of 12 bytes */
using SerializedConfiguration = std::array<std::byte, 12>;

template <typename SerializedType>
using OptionsPerRank = std::unique_ptr<std::set<SerializedType>[]>;

template <typename TOption>
std::byte castToByte(TOption option) {
  return static_cast<std::byte>(static_cast<typename TOption::Value>(option));
}

template <typename SerializedType>
void prepareForCommunication(const OptionsPerRank<SerializedType> &optionsPerRank, const int commSize,
                             std::unique_ptr<int[]> &displs, std::unique_ptr<int[]> &sendCounts, int &size,
                             std::vector<SerializedType> &optionsForAllRanks) {
  int displacementAcc = 0;
  for (int rankIndex = 0; rankIndex < commSize; ++rankIndex) {
    displs[rankIndex] = displacementAcc;
    sendCounts[rankIndex] = optionsPerRank[rankIndex].size();
    displacementAcc += optionsPerRank[rankIndex].size();
    for (auto option : optionsPerRank[rankIndex]) {
      optionsForAllRanks.push_back(option);
    }
  }
  size = displacementAcc;
}

template <typename SerializedType>
AutoPas_MPI_Datatype templateToMPI();
template <>
AutoPas_MPI_Datatype templateToMPI<std::byte>() { return AUTOPAS_MPI_BYTE; }
template <>
AutoPas_MPI_Datatype templateToMPI<double>() { return AUTOPAS_MPI_DOUBLE; }

template <typename TOption, typename SerializedType>
void distributeOption(const AutoPas_MPI_Comm comm, const int commSize, const int rank, std::set<TOption> &localOptions,
                      OptionsPerRank<SerializedType> &optionsPerRank) {
  auto optionsForAllRanks = std::vector<SerializedType>();
  auto displs = std::make_unique<int[]>(commSize);
  auto sendCounts = std::make_unique<int[]>(commSize);
  int arraySize;
  prepareForCommunication<SerializedType>(optionsPerRank, commSize, displs, sendCounts, arraySize, optionsForAllRanks);

  auto recvBuf = std::vector<SerializedType>(sendCounts[rank]);
  MPI_Scatterv(&optionsForAllRanks.data(), sendCounts, displs, templateToMPI<SerializedType>(), &recvBuf.data(),
               sendCounts[rank], templateToMPI<SerializedType>(), 0, comm);

  localOptions.clear();
  for (int i = 0; i < sendCounts[rank]; ++i) {
    localOptions.emplace(recvBuf[i]);
  }
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
void distributeConfigurations(std::set<ContainerOption> &containerOptions, std::set<double> &cellSizeFactors,
                              std::set<TraversalOption> &traversalOptions, std::set<DataLayoutOption> &dataLayoutOptions,
                              std::set<Newton3Option> &newton3Options, AutoPas_MPI_Comm comm);

/**
 * Serializes a configuration object for communication via MPI
 * @param configuration: the configuration to be sent
 * @return The serialization
 */
SerializedConfiguration serializeConfiguration(Configuration configuration) {
  // @todo maybe consider endianness for different processors
  SerializedConfiguration config;
  config[0] = static_cast<std::byte>(static_cast<ContainerOption::Value>(configuration.container));
  config[1] = static_cast<std::byte>(static_cast<TraversalOption::Value>(configuration.traversal));
  config[2] = static_cast<std::byte>(static_cast<DataLayoutOption::Value>(configuration.dataLayout));
  config[3] = static_cast<std::byte>(static_cast<Newton3Option::Value>(configuration.newton3));
  std::memcpy(&config[4], &configuration.cellSizeFactor, sizeof(double));
  return config;
}

/**
 * Recreates a Configuration object from the object obtained by _serializeConfiguration
 * @param config: The SerializedConfiguration objects returned by _serializeConfiguration
 * @return The deserialized Configuration object
 */
static Configuration deserializeConfig(SerializedConfiguration config) {
  double cellSizeFactor;
  std::memcpy(&cellSizeFactor, &config[4], sizeof(double));
  return Configuration(static_cast<ContainerOption::Value>(config[0]), cellSizeFactor,
                       static_cast<TraversalOption::Value>(config[1]),
                       static_cast<DataLayoutOption::Value>(config[2]), static_cast<Newton3Option::Value>(config[3]));
}

/**
 * Helper class for FullSearchMPI.
 * Handles serialization and deserialization of configurations, as well selecting the global optimal via MPI calls
 */
/**
 * Handles communication to select the globally best configuration.
 * @param comm: The communicator used for sending and receiving the optimal configuration
 * @param localOptimalConfig: The locally optimal configuration to be compared with others
 * @param localOptimalTime: The time measured for localOptimalConfig
 * @return The globally optimal configuration
 */
Configuration optimizeConfiguration(AutoPas_MPI_Comm comm, Configuration localOptimalConfig,
                                    size_t localOptimalTime) {
  SerializedConfiguration serializedConfiguration = serializeConfiguration(localOptimalConfig);
  size_t optimalTimeOut;
  int optimalRankIn, optimalRankOut;

  AutoPas_MPI_Allreduce(&localOptimalTime, &optimalTimeOut, 1, AUTOPAS_MPI_UNSIGNED_LONG, AUTOPAS_MPI_MIN, comm);

  // Send own rank if local optimal time is equal to the global optimal time.
  // Send something higher than the highest rank otherwise.
  if (localOptimalTime == optimalTimeOut) {
    AutoPas_MPI_Comm_rank(comm, &optimalRankIn);
  } else {
    AutoPas_MPI_Comm_size(comm, &optimalRankIn);
  }
  AutoPas_MPI_Allreduce(&optimalRankIn, &optimalRankOut, 1, AUTOPAS_MPI_INT, AUTOPAS_MPI_MIN, comm);

  AutoPas_MPI_Bcast(&serializedConfiguration, sizeof(serializedConfiguration), AUTOPAS_MPI_BYTE, optimalRankOut,
                    comm);

  return deserializeConfig(serializedConfiguration);
}

}  // namespace autopas
