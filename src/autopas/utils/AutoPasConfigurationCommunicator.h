/**
 * @file AutoPasConfigurationCommunicator.h
 * @author W. Thieme
 * @date 04/30/20
 */

#include "WrapMPI.h"
#include "autopas/selectors/Configuration.h"
#include "autopas/utils/ExceptionHandler.h"
#include <array>

namespace autopas {

using SerializedConfiguration = std::array<std::byte, 12>;

class AutoPasConfigurationCommunicator {
public:
  /**
   * Handles communication to select the globally best configuration.
   * @param localOptimalConfig: The locally optimal configuration to be compared with others
   * @param localOptimalTime: The time measured for localOptimalConfig
   * @return The globally optimal configuration
   */
  Configuration optimizeConfiguration(MPI_Comm comm, Configuration localOptimalConfig, size_t localOptimalTime) {
    _serializedConfiguration = serializeConfiguration(localOptimalConfig);
    _optimalTime = localOptimalTime;
    MPI_Allreduce(&_optimalTime, &_optimalTime, 1, MPI_UNSIGNED_LONG, MPI_MIN, comm);

    // Send own rank if local optimal time is equal to the global optimal time.
    // Send something higher than the highest rank otherwise.
    if (localOptimalTime == _optimalTime) {
      MPI_Comm_rank(comm, &_optimalRank);
    } else {
      MPI_Comm_size(comm, &_optimalRank);
    }
    MPI_Allreduce(&_optimalRank, &_optimalRank, 1, MPI_INT, MPI_MIN, comm);

    MPI_Bcast(&_serializedConfiguration, sizeof(_serializedConfiguration), MPI_BYTE, _optimalRank, comm);

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
    return Configuration(static_cast<ContainerOption::Value>(config[0]),
                         cellSizeFactor,
                         static_cast<TraversalOption::Value>(config[1]),
                         static_cast<DataLayoutOption::Value>(config[2]),
                         static_cast<Newton3Option::Value>(config[3]));
  }

  SerializedConfiguration _serializedConfiguration;
  size_t _optimalTime;
  int _optimalRank;
};

}
