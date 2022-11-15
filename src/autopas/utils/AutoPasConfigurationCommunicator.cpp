/**
 * @file AutoPasConfigurationCommunicator.cpp
 * @author W. Thieme
 * @date 29.05.2020
 */

#include "AutoPasConfigurationCommunicator.h"

#include "ThreeDimensionalMapping.h"
#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/utils/ConfigurationAndRankIteratorHandler.h"
#include "autopas/utils/logging/Logger.h"

namespace autopas::utils::AutoPasConfigurationCommunicator {

size_t getSearchSpaceSize(const std::set<ContainerOption> &containerOptions, const NumberSet<double> &cellSizeFactors,
                          const std::set<TraversalOption> &traversalOptions,
                          const std::set<LoadEstimatorOption> &loadEstimatorOptions,
                          const std::set<DataLayoutOption> &dataLayoutOptions,
                          const std::set<Newton3Option> &newton3Options,
                          const NumberSet<int> &verletRebuildFrequencies) {
  size_t numConfigs = 0;
  // only take into account finite sets of cellSizeFactors.
  size_t cellSizeFactorArraySize;
  if (cellSizeFactors.isFinite()) {
    cellSizeFactorArraySize = cellSizeFactors.size();
  } else {
    cellSizeFactorArraySize = 1;
  }

  size_t verletRebuildFrequencyArraySize;
  if (verletRebuildFrequencies.isFinite()) {
    verletRebuildFrequencyArraySize = verletRebuildFrequencies.size();
  } else {
    verletRebuildFrequencyArraySize = 1;
  }

  for (const auto &containerOption : containerOptions) {
    // get all traversals of the container and restrict them to the allowed ones.
    const std::set<TraversalOption> &allContainerTraversals =
        compatibleTraversals::allCompatibleTraversals(containerOption);
    std::set<TraversalOption> allowedAndApplicableTraversalOptions;
    std::set_intersection(
        traversalOptions.begin(), traversalOptions.end(), allContainerTraversals.begin(), allContainerTraversals.end(),
        std::inserter(allowedAndApplicableTraversalOptions, allowedAndApplicableTraversalOptions.begin()));

    for (const auto &traversalOption : allowedAndApplicableTraversalOptions) {
      // if load estimators are not applicable LoadEstimatorOption::none is returned.
      const std::set<LoadEstimatorOption> allowedAndApplicableLoadEstimators =
          loadEstimators::getApplicableLoadEstimators(containerOption, traversalOption, loadEstimatorOptions);
      numConfigs += cellSizeFactorArraySize * allowedAndApplicableLoadEstimators.size() * dataLayoutOptions.size() *
                    newton3Options.size() * verletRebuildFrequencyArraySize;
    }
  }
  return numConfigs;
}

/**
 * Calculates which Options the current rank should handle based on the total number of options and ranks.
 * @param numConfigs in
 * @param commSize in
 * @param rank in
 * @param containerOptions inout
 * @param cellSizeFactors inout
 * @param traversalOptions inout
 * @param loadEstimatorOptions inout
 * @param dataLayoutOptions inout
 * @param newton3Options inout
 */
void generateDistribution(const int numConfigs, const int commSize, const int rank,
                          std::set<ContainerOption> &containerOptions, NumberSet<double> &cellSizeFactors,
                          std::set<TraversalOption> &traversalOptions,
                          std::set<LoadEstimatorOption> &loadEstimatorOptions,
                          std::set<DataLayoutOption> &dataLayoutOptions, std::set<Newton3Option> &newton3Options,
                          NumberSet<int> &verletRebuildFrequencies) {
  // ============== setup ======================================================

  // These will be set to the Options specific to this rank and will overwrite the input sets.
  auto newContainerOptions = std::set<ContainerOption>();
  auto newCellSizeFactors = std::set<double>();
  auto newTraversalOptions = std::set<TraversalOption>();
  auto newLoadEstimatorOptions = std::set<LoadEstimatorOption>();
  auto newDataLayoutOptions = std::set<DataLayoutOption>();
  auto newNewton3Options = std::set<Newton3Option>();
  auto newVerletRebuildFrequency = std::set<int>();

  // Distribution works only with finite sets of cellSizeFactors.
  // If the set is infinite a dummy value will be used and replaced later on.
  std::set<double> finiteCellSizeFactors;
  if (cellSizeFactors.isFinite()) {
    finiteCellSizeFactors = cellSizeFactors.getAll();
  } else {
    // Dummy value which makes the code simpler in case the cellSizeFactors are not a finite set.
    finiteCellSizeFactors = std::set<double>{-1};
  }

  //TODO verletRebuildFrequency is Finite?

  // ============== main computation ===========================================

  ConfigurationAndRankIteratorHandler iteratorHandler(containerOptions, finiteCellSizeFactors, traversalOptions,
                                                      loadEstimatorOptions, dataLayoutOptions, newton3Options,
                                                      verletRebuildFrequencies.getAll(),
                                                      numConfigs, commSize);

  while (iteratorHandler.getRankIterator() < rank) {
    iteratorHandler.advanceIterators(numConfigs, commSize);
  }

  // Only important for infinite cellSizeFactors if commSize > numConfigs.
  int infiniteCellSizeFactorsOffset = iteratorHandler.getInfiniteCellSizeFactorsOffset();
  int infiniteCellSizeFactorsBlockSize = iteratorHandler.getInfiniteCellSizeFactorsBlockSize();

  while (iteratorHandler.getRankIterator() == rank) {
    // std::set handles duplicate elements.
    newContainerOptions.emplace(*iteratorHandler.getContainerIterator());
    newCellSizeFactors.emplace(*iteratorHandler.getCellSizeFactorIterator());
    newTraversalOptions.emplace(*iteratorHandler.getTraversalIterator());
    newLoadEstimatorOptions.emplace(*iteratorHandler.getLoadEstimatorIterator());
    newDataLayoutOptions.emplace(*iteratorHandler.getDataLayoutIterator());
    newNewton3Options.emplace(*iteratorHandler.getNewton3Iterator());
    newVerletRebuildFrequency.emplace(*iteratorHandler.getVerletRebuildFrequenzyIterator());

    iteratorHandler.advanceIterators(numConfigs, commSize);
  }

  // ============== assigning to local search space ============================

  containerOptions = newContainerOptions;
  if (not cellSizeFactors.isFinite()) {
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
  loadEstimatorOptions = newLoadEstimatorOptions;
  dataLayoutOptions = newDataLayoutOptions;
  newton3Options = newNewton3Options;
}

void distributeConfigurations(std::set<ContainerOption> &containerOptions, NumberSet<double> &cellSizeFactors,
                              std::set<TraversalOption> &traversalOptions,
                              std::set<LoadEstimatorOption> &loadEstimatorOptions,
                              std::set<DataLayoutOption> &dataLayoutOptions, std::set<Newton3Option> &newton3Options,
                              NumberSet<int> &verletRebuildFrequencies,
                              const int rank, const int commSize) {
  int numConfigs = getSearchSpaceSize(containerOptions, cellSizeFactors, traversalOptions, loadEstimatorOptions,
                                      dataLayoutOptions, newton3Options, verletRebuildFrequencies);

  if (numConfigs == 0) {
    utils::ExceptionHandler::exception("Could not generate valid configurations, aborting");
    return;
  }

  // Creates a set for each option and each rank containing the serialized versions (std::byte or double) of all
  // options assigned to that rank.
  generateDistribution(numConfigs, commSize, rank, containerOptions, cellSizeFactors, traversalOptions,
                       loadEstimatorOptions, dataLayoutOptions, newton3Options, verletRebuildFrequencies);

  size_t cellSizeFactorsSize = cellSizeFactors.isFinite() ? cellSizeFactors.size() : 1;
  AutoPasLog(debug,
             "After distributing: {} containers, {} cellSizeFactors, {} traversals, {} dataLayouts, {} newton3s"
             " => {} total configs",
             containerOptions.size(), cellSizeFactorsSize, traversalOptions.size(), dataLayoutOptions.size(),
             newton3Options.size(),
             getSearchSpaceSize(containerOptions, cellSizeFactors, traversalOptions, loadEstimatorOptions,
                                dataLayoutOptions, newton3Options, verletRebuildFrequencies));
}

Configuration optimizeConfiguration(AutoPas_MPI_Comm comm, Configuration localOptimalConfig, size_t localOptimalTime) {
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

  AutoPas_MPI_Bcast(serializedConfiguration.data(), serializedConfiguration.size(), AUTOPAS_MPI_BYTE, optimalRankOut,
                    comm);

  Configuration deserializedConfig = deserializeConfiguration(serializedConfiguration);
  AutoPasLog(debug, "Globally optimal configuration: {}", deserializedConfig.toString());

  return deserializedConfig;
}

SerializedConfiguration serializeConfiguration(Configuration configuration) {
  // @todo maybe consider endianness for different processors
  SerializedConfiguration config;
  config[0] = castToByte(configuration.container);
  config[1] = castToByte(configuration.traversal);
  config[2] = castToByte(configuration.loadEstimator);
  config[3] = castToByte(configuration.dataLayout);
  config[4] = castToByte(configuration.newton3);
  std::memcpy(&config[5], &configuration.cellSizeFactor, sizeof(double));
  std::memcpy(&config[6], &configuration.cellSizeFactor, sizeof(int));
  return config;
}

Configuration deserializeConfiguration(SerializedConfiguration config) {
  double cellSizeFactor;
  int verletRebuildFrequency;
  std::memcpy(&cellSizeFactor, &config[5], sizeof(double));
  std::memcpy(&verletRebuildFrequency, &config[6], sizeof(int));
  return Configuration(static_cast<ContainerOption::Value>(config[0]), cellSizeFactor,
                       static_cast<TraversalOption::Value>(config[1]),
                       static_cast<LoadEstimatorOption::Value>(config[2]),
                       static_cast<DataLayoutOption::Value>(config[3]), static_cast<Newton3Option::Value>(config[4]),
                       verletRebuildFrequency);
}

}  // namespace autopas::utils::AutoPasConfigurationCommunicator
