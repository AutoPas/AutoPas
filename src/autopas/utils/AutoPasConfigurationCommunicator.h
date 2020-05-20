/**
 * @file AutoPasConfigurationCommunicator.h
 * @author W. Thieme
 * @date 30.04.2020
 */

#pragma once

#include <array>

#include "WrapMPI.h"
#include "autopas/selectors/Configuration.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

/** type definition for the serialization of configurations. A serialized config is an array of 12 bytes */
using SerializedConfiguration = std::array<std::byte, 12>;

/**
 * Helper class for FullSearchMPI.
 * Handles serialization and deserialization of configurations, as well selecting the global optimal via MPI calls
 */
class AutoPasConfigurationCommunicator {
 public:
  /**
   * Handles communication to select the globally best configuration.
   * @param comm: The communicator used for sending and receiving the optimal configuration
   * @param localOptimalConfig: The locally optimal configuration to be compared with others
   * @param localOptimalTime: The time measured for localOptimalConfig
   * @return The globally optimal configuration
   */
  Configuration optimizeConfiguration(AutoPas_MPI_Comm comm, Configuration localOptimalConfig,
                                      size_t localOptimalTime) {
    _serializedConfiguration = serializeConfiguration(localOptimalConfig);
    _optimalTimeIn = localOptimalTime;
    AutoPas_MPI_Allreduce(&_optimalTimeIn, &_optimalTimeOut, 1, AUTOPAS_MPI_UNSIGNED_LONG, AUTOPAS_MPI_MIN, comm);

    // Send own rank if local optimal time is equal to the global optimal time.
    // Send something higher than the highest rank otherwise.
    if (_optimalTimeIn == _optimalTimeOut) {
      AutoPas_MPI_Comm_rank(comm, &_optimalRankIn);
    } else {
      AutoPas_MPI_Comm_size(comm, &_optimalRankIn);
    }
    AutoPas_MPI_Allreduce(&_optimalRankIn, &_optimalRankOut, 1, AUTOPAS_MPI_INT, AUTOPAS_MPI_MIN, comm);

    AutoPas_MPI_Bcast(&_serializedConfiguration, sizeof(_serializedConfiguration), AUTOPAS_MPI_BYTE, _optimalRankOut,
                      comm);

    return deserializeConfig(_serializedConfiguration);
  }

 private:
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

  SerializedConfiguration _serializedConfiguration;
  size_t _optimalTimeIn, _optimalTimeOut;
  int _optimalRankIn, _optimalRankOut;
};

}  // namespace autopas
