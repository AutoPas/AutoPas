/**
 * @file AutoPasConfigurationPasser.h
 * @author W. Thieme
 * @date 04/30/20
 */

#include "WrapMPI.h"
#include "autopas/selectors/Configuration.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

struct TraversalTimeAndRank {
  size_t time;
  int rank;
};

struct Config {
  int container, traversal, dataLayout, newton3;
  double cellSizeFactor;
};

class AutoPasConfigurationPasser {
public:
  /**
   * Constructor for AutoPasConfigurationPasser
   * @param configuration: local optimal configuration
   */
  explicit AutoPasConfigurationPasser(Configuration configuration)
    : _request(nullptr),
      _config(_serializeConfiguration(configuration)),
      _optimalTraversalTimeAndRank({0,-1}),
      _state(0) {

  }

  /**
   * Handles communication to select the globally best configuration.
   * To be called repeatedly.
   * @param configuration: outputs optimal configuration
   * @return if true the function will have to be called again. Else the optimization is finished.
   */
  bool optimizeConfiguration(Configuration  *configuration, size_t traversalTime) {
    int completed;
    int rank;
    AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &rank);
    switch (_state) {
      case 0:
        _findGlobalOptimalTime(traversalTime);
        ++_state;
        return true;
      case 1:
        AutoPas_MPI_Test(&_request, &completed, AUTOPAS_MPI_STATUS_IGNORE);
        if (completed) {
          _findGlobalOptimalRank(traversalTime);
          ++_state;
        }
        return true;
      case 2:
        AutoPas_MPI_Test(&_request, &completed, AUTOPAS_MPI_STATUS_IGNORE);
        if (completed) {
          _broadcastGlobalOptimalConfiguration();
          ++_state;
        }
        return true;
      case 3:
        AutoPas_MPI_Test(&_request, &completed, AUTOPAS_MPI_STATUS_IGNORE);
        if (completed) {
          *configuration = _getGlobalOptimalConfiguration();
          ++_state;
          return false;
        }
        return true;
      default:
        autopas::utils::ExceptionHandler::exception("optimizeConfiguration called after the optimization was already finished");
        return false;
    }
  }

private:
  /**
   * Find the globally optimal configuration based on traversal times provided by each rank.
   * @param traversalTime: time of the locally optimal configuration
   * @param outRank: The rank of the process which has the globally optimal configuration.
   *                 May only be accessed after AutoPas_MPI_Wait or AutoPas_MPI_Test
   * @return handle to a AutoPas_MPI_Request.
   */
  void _findGlobalOptimalTime(size_t traversalTime) {
    int rank;
    AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &rank);

    std::cout << "sending " << traversalTime << " as " << rank << std::endl;
    AutoPas_MPI_Iallreduce(&traversalTime, &_optimalTraversalTimeAndRank.time, 1, MPI_UNSIGNED_LONG, MPI_MIN,
                           AUTOPAS_MPI_COMM_WORLD, &_request);
  }


  void _findGlobalOptimalRank(size_t traversalTime) {
    int rank;

    // Send own rank if local optimal time is equal to the global optimal time.
    // Send something higher than the highest rank otherwise.
    if (traversalTime == _optimalTraversalTimeAndRank.time) {
      AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &rank);
    } else {
      AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &rank);
    }
    std::cout << "the best time was " << _optimalTraversalTimeAndRank.time << std::endl;

    std::cout << "sending rank " << rank << " as optimum" << std::endl;
    AutoPas_MPI_Iallreduce(&rank, &_optimalTraversalTimeAndRank.rank, 1, MPI_INT, MPI_MIN, AUTOPAS_MPI_COMM_WORLD,
                           &_request);
  }

  /**
   * Broadcasts the best configuration to all ranks
   * @param configuration: inout for globally optimal configuration.
   *                       May only be accessed after AutoPas_MPI_Wait or AutoPas_MPI_Test.
   * @param bestRank: The rank containing the optimal configuration.
   * @return handle to a AutoPas_MPI_Request.
   */
  void _broadcastGlobalOptimalConfiguration() {
    int rank;
    AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &rank);

    std::cout << "bcast from " << _optimalTraversalTimeAndRank.rank << std::endl;
    AutoPas_MPI_Ibcast(&_config, sizeof(_config), AUTOPAS_MPI_BYTE, _optimalTraversalTimeAndRank.rank,
                       AUTOPAS_MPI_COMM_WORLD, &_request);
  }

  Configuration _getGlobalOptimalConfiguration() {
    std::cout << "best config arrived: " << _config.container << _config.traversal << _config.dataLayout << _config.newton3 << _config.cellSizeFactor << std::endl;
    return _deserializeConfig(_config);
  }

  /**
   * Serializes a configuration object for communication via MPI
   * @param configuration: the configuration to be sent
   * @return The serialization
   */
  static Config _serializeConfiguration(Configuration configuration) {
    Config config{configuration.container,
                  configuration.traversal,
                  configuration.dataLayout,
                  configuration.newton3,
                  configuration.cellSizeFactor};
    int rank;
    AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &rank);
    std::cout << rank << ":configuration: " << config.container << ", " << config.traversal << ", " << config.dataLayout << ", " << config.newton3 << ", " << config.cellSizeFactor << std::endl;
    return config;
  }

  static Configuration _deserializeConfig(Config config) {
    return Configuration(static_cast<ContainerOption::Value>(config.container),
                         config.cellSizeFactor,
                         static_cast<TraversalOption::Value>(config.traversal),
                         static_cast<DataLayoutOption::Value>(config.dataLayout),
                         static_cast<Newton3Option::Value>(config.newton3));
  }

  AutoPas_MPI_Request _request;
  Config _config;
  TraversalTimeAndRank _optimalTraversalTimeAndRank;
  int _state;
};

}
