/**
 * @file AutoPasConfigurationCommunicator.cpp
 * @author W. Thieme
 * @date 29.05.2020
 */

#include "AutoPasConfigurationCommunicator.h"
#include "autopas/utils/Logger.h"

namespace autopas {

inline void AutoPasConfigurationCommunicator::IteratorHandler::advanceConfigIterators() {
  // advance to the next valid config
  ++_newton3It;
  if (_newton3It != _newton3Options.end()) return;
  _newton3It = _newton3Options.begin();
  ++_dataLayoutIt;
  if (_dataLayoutIt != _dataLayoutOptions.end()) return;
  _dataLayoutIt = _dataLayoutOptions.begin();
  ++_traversalIt;
  if (_traversalIt != _allowedAndApplicableTraversalOptions.end()) return;
  _traversalIt = _allowedAndApplicableTraversalOptions.begin();
  ++_cellSizeFactorIt;
  if (_cellSizeFactorIt != _cellSizeFactors.end()) return;
  _cellSizeFactorIt = _cellSizeFactors.begin();
  ++_containerIt;
  selectTraversalsForCurrentContainer();
}

void AutoPasConfigurationCommunicator::IteratorHandler::advanceIterators(const int numConfigs, const int commSize) {
  if (numConfigs >= commSize or _remainingBlockSize == 0) {
    advanceConfigIterators();
  }

  if (commSize >= numConfigs or _remainingBlockSize == 0) {
    // advance to the next rank
    ++_rankIterator;

    // advance offset to the position relative to the first rank with the same configuration
    ++_infiniteCellSizeFactorsOffset;
  }

  // Set information necessary to compute the next block.
  // Block here means either a block of ranks that all have the same configuration or a set of configuration that all
  // have the same ranks.
  if (_remainingBlockSize == 0) {
    if (numConfigs >= commSize) {
      _remainingBlockSize = numConfigs / commSize;
    } else {
      _remainingBlockSize = commSize / numConfigs;

      _infiniteCellSizeFactorsBlockSize = _remainingBlockSize;
      _infiniteCellSizeFactorsOffset = 0;
    }
    if (_remainder > 0) {
      ++_remainingBlockSize;
      --_remainder;
    }
  }

  --_remainingBlockSize;
}

void AutoPasConfigurationCommunicator::IteratorHandler::selectTraversalsForCurrentContainer() {
  // get all traversals of the container and restrict them to the allowed ones
  const std::set<TraversalOption> &allContainerTraversals =
      compatibleTraversals::allCompatibleTraversals(*_containerIt);
  std::set_intersection(_allowedTraversalOptions.begin(), _allowedTraversalOptions.end(),
                        allContainerTraversals.begin(), allContainerTraversals.end(),
                        std::inserter(_allowedAndApplicableTraversalOptions,
                                             _allowedAndApplicableTraversalOptions.begin()));
  _traversalIt = _allowedAndApplicableTraversalOptions.begin();
}

/**
 * Calculates the maximum number of valid configs from several sets of options.
 * This does not equal the cartesian product as not all containers are compatible with all traversals.
 * @param containerOptions
 * @param cellSizeFactors The size of cellSizeFactors will only be taken into account if the NumberSet is finite
 * @param traversalOptions
 * @param dataLayoutOptions
 * @param newton3Options
 * @return
 */
size_t getSearchSpaceSize(std::set<ContainerOption> &containerOptions, NumberSet<double> &cellSizeFactors,
                          std::set<TraversalOption> &traversalOptions, std::set<DataLayoutOption> &dataLayoutOptions,
                          std::set<Newton3Option> &newton3Options) {
  size_t numConfigs = 0;
  // only take into account finite sets of cellSizeFactors
  size_t cellSizeFactorArraySize;
  if (cellSizeFactors.isFinite()) {
    cellSizeFactorArraySize = cellSizeFactors.size();
  } else {
    cellSizeFactorArraySize = 1;
  }

  for (const auto &containerOption : containerOptions) {
    // get all traversals of the container and restrict them to the allowed ones
    const std::set<TraversalOption> &allContainerTraversals =
        compatibleTraversals::allCompatibleTraversals(containerOption);
    std::set<TraversalOption> allowedAndApplicable;
    std::set_intersection(traversalOptions.begin(), traversalOptions.end(), allContainerTraversals.begin(),
                          allContainerTraversals.end(),
                          std::inserter(allowedAndApplicable, allowedAndApplicable.begin()));
    numConfigs +=
        cellSizeFactorArraySize * allowedAndApplicable.size() * dataLayoutOptions.size() * newton3Options.size();
  }
  return numConfigs;
}

/**
 * Calculates which Options the current rank should handle based on the total number of options and ranks
 * @param numConfigs in
 * @param commSize in
 * @param rank in
 * @param containerOptions inout
 * @param cellSizeFactors inout
 * @param traversalOptions inout
 * @param dataLayoutOptions inout
 * @param newton3Options inout
 */
void generateDistribution(const int numConfigs, const int commSize, const int rank,
                          std::set<ContainerOption> &containerOptions, NumberSet<double> &cellSizeFactors,
                          std::set<TraversalOption> &traversalOptions, std::set<DataLayoutOption> &dataLayoutOptions,
                          std::set<Newton3Option> &newton3Options) {
  // ============== setup ======================================================

  // These will be set to the Options specific to this rank and will overwrite the input sets
  auto newContainerOptions = std::set<ContainerOption>();
  auto newCellSizeFactors = std::set<double>();
  auto newTraversalOptions = std::set<TraversalOption>();
  auto newDataLayoutOptions = std::set<DataLayoutOption>();
  auto newNewton3Options = std::set<Newton3Option>();

  // Distribution works only with finite sets of cellSizeFactors.
  // If the set is infinite a dummy value will be used and replaced later on.
  std::set<double> finiteCellSizeFactors;
  if (cellSizeFactors.isFinite()) {
    finiteCellSizeFactors = cellSizeFactors.getAll();
  } else {
    // Dummy value which makes the code simpler in case the cellSizeFactors are not a finite set.
    finiteCellSizeFactors = std::set<double>{-1};
  }

  // ============== main computation ===========================================

  AutoPasConfigurationCommunicator::IteratorHandler iteratorHandler(containerOptions, finiteCellSizeFactors,
                                                                    traversalOptions, dataLayoutOptions, newton3Options,
                                                                    numConfigs, commSize);

  while (iteratorHandler.getRankIterator() < rank) {
    iteratorHandler.advanceIterators(numConfigs, commSize);
  }

  // Only important for infinite cellSizeFactors if commSize > numConfigs
  int infiniteCellSizeFactorsOffset = iteratorHandler.getInfiniteCellSizeFactorsOffset();
  int infiniteCellSizeFactorsBlockSize = iteratorHandler.getInfiniteCellSizeFactorsBlockSize();

  while (iteratorHandler.getRankIterator() == rank) {
    // std::set handles duplicate elements
    newContainerOptions.emplace(*iteratorHandler.getContainerIterator());
    newCellSizeFactors.emplace(*iteratorHandler.getCellSizeFactorIterator());
    newTraversalOptions.emplace(*iteratorHandler.getTraversalIterator());
    newDataLayoutOptions.emplace(*iteratorHandler.getDataLayoutIterator());
    newNewton3Options.emplace(*iteratorHandler.getNewton3Iterator());

    iteratorHandler.advanceIterators(numConfigs, commSize);
  }

  // ============== assigning to local search space ============================

  containerOptions = newContainerOptions;
  if (not cellSizeFactors.isFinite()) {
    // newCellSizeFactors == {-1}
    double min = cellSizeFactors.getMin();
    double max = cellSizeFactors.getMax();
    double delta = (max - min) / infiniteCellSizeFactorsBlockSize;
    std::set<double> values{min + delta * infiniteCellSizeFactorsOffset,
                            min + delta * (infiniteCellSizeFactorsOffset + 1)};
    cellSizeFactors.resetValues(values);
  } else {
    cellSizeFactors.resetValues(newCellSizeFactors);
  }
  traversalOptions = newTraversalOptions;
  dataLayoutOptions = newDataLayoutOptions;
  newton3Options = newNewton3Options;
}

void AutoPasConfigurationCommunicator::distributeConfigurations(std::set<ContainerOption> &containerOptions,
                                                                NumberSet<double> &cellSizeFactors,
                                                                std::set<TraversalOption> &traversalOptions,
                                                                std::set<DataLayoutOption> &dataLayoutOptions,
                                                                std::set<Newton3Option> &newton3Options,
                                                                AutoPas_MPI_Comm comm) {
  int rank, commSize;
  AutoPas_MPI_Comm_rank(comm, &rank);
  AutoPas_MPI_Comm_size(comm, &commSize);

  int numConfigs =
      getSearchSpaceSize(containerOptions, cellSizeFactors, traversalOptions, dataLayoutOptions, newton3Options);

  if (numConfigs == 0) {
    utils::ExceptionHandler::exception("Could not generate valid configurations, aborting");
    return;
  }

  // Creates a set for each option and each rank containing the serialized versions (std::byte or double) of all
  // options assigned to that rank.
  generateDistribution(numConfigs, commSize, rank, containerOptions, cellSizeFactors, traversalOptions,
                       dataLayoutOptions, newton3Options);

  size_t cellSizeFactorsSize = cellSizeFactors.isFinite() ? cellSizeFactors.size() : 1;
  //AutoPasLog(debug, "After distributing {} containers, {} cellSizeFactors, {} traversals, {} dataLayouts, {} newton3s",
  //           containerOptions.size(), cellSizeFactorsSize, traversalOptions.size(),
  //           dataLayoutOptions.size(), newton3Options.size());
}

Configuration AutoPasConfigurationCommunicator::optimizeConfiguration(AutoPas_MPI_Comm comm,
                                                                      Configuration localOptimalConfig,
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

  AutoPas_MPI_Bcast(&serializedConfiguration, sizeof(serializedConfiguration), AUTOPAS_MPI_BYTE, optimalRankOut, comm);

  Configuration deserializedConfig = deserializeConfiguration(serializedConfiguration);
  //AutoPasLog(debug, "Globally best configuration: {}", deserializedConfig.toString());

  return deserializedConfig;
}

AutoPasConfigurationCommunicator::SerializedConfiguration AutoPasConfigurationCommunicator::serializeConfiguration(
    Configuration configuration) {
  // @todo maybe consider endianness for different processors
  SerializedConfiguration config;
  config[0] = castToByte(configuration.container);
  config[1] = castToByte(configuration.traversal);
  config[2] = castToByte(configuration.dataLayout);
  config[3] = castToByte(configuration.newton3);
  std::memcpy(&config[4], &configuration.cellSizeFactor, sizeof(double));
  return config;
}

Configuration AutoPasConfigurationCommunicator::deserializeConfiguration(SerializedConfiguration config) {
  double cellSizeFactor;
  std::memcpy(&cellSizeFactor, &config[4], sizeof(double));
  return Configuration(static_cast<ContainerOption::Value>(config[0]), cellSizeFactor,
                       static_cast<TraversalOption::Value>(config[1]), static_cast<DataLayoutOption::Value>(config[2]),
                       static_cast<Newton3Option::Value>(config[3]));
}

}  // namespace autopas
