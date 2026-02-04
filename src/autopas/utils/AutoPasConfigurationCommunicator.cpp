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
                          const InteractionTypeOption &interactionTypeOption,
                          const NumberSetFinite<int> &threadCounts) {
  // only take into account finite sets of cellSizeFactors.
  const size_t cellSizeFactorArraySize = cellSizeFactors.isFinite() ? cellSizeFactors.size() : 1;

  size_t numConfigs{0};
  for (const auto &containerOption : containerOptions) {
    // get all traversals of the container and restrict them to the allowed ones.
    const std::set<TraversalOption> &allContainerTraversals =
        compatibleTraversals::allCompatibleTraversals(containerOption, interactionTypeOption);
    std::set<TraversalOption> allowedAndApplicableTraversalOptions;
    std::set_intersection(
        traversalOptions.begin(), traversalOptions.end(), allContainerTraversals.begin(), allContainerTraversals.end(),
        std::inserter(allowedAndApplicableTraversalOptions, allowedAndApplicableTraversalOptions.begin()));

    for (const auto &traversalOption : allowedAndApplicableTraversalOptions) {
      // if load estimators are not applicable LoadEstimatorOption::none is returned.
      const std::set<LoadEstimatorOption> allowedAndApplicableLoadEstimators =
          loadEstimators::getApplicableLoadEstimators(containerOption, traversalOption, loadEstimatorOptions);
      numConfigs += cellSizeFactorArraySize * allowedAndApplicableLoadEstimators.size() * dataLayoutOptions.size() *
                    newton3Options.size() * threadCounts.size();
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
 * @param interactionType Handle configurations for this interaction type
 */
void generateDistribution(const int numConfigs, const int commSize, const int rank,
                          std::set<ContainerOption> &containerOptions, NumberSet<double> &cellSizeFactors,
                          std::set<TraversalOption> &traversalOptions,
                          std::set<LoadEstimatorOption> &loadEstimatorOptions,
                          std::set<DataLayoutOption> &dataLayoutOptions, std::set<Newton3Option> &newton3Options,
                          const InteractionTypeOption interactionType, NumberSetFinite<int> &threadCounts) {
  // ============== setup ======================================================

  // These will be set to the Options specific to this rank and will overwrite the input sets.
  auto newContainerOptions = std::set<ContainerOption>();
  auto newCellSizeFactors = std::set<double>();
  auto newTraversalOptions = std::set<TraversalOption>();
  auto newLoadEstimatorOptions = std::set<LoadEstimatorOption>();
  auto newDataLayoutOptions = std::set<DataLayoutOption>();
  auto newNewton3Options = std::set<Newton3Option>();
  auto newThreadCounts = std::set<int>();

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

  ConfigurationAndRankIteratorHandler iteratorHandler(containerOptions, finiteCellSizeFactors, traversalOptions,
                                                      loadEstimatorOptions, dataLayoutOptions, newton3Options,
                                                      interactionType, threadCounts.getAll(), numConfigs, commSize);

  while (iteratorHandler.getRankIterator() < rank) {
    iteratorHandler.advanceIterators(numConfigs, commSize);
  }

  // Only important for infinite cellSizeFactors if commSize > numConfigs.
  const int infiniteCellSizeFactorsOffset = iteratorHandler.getInfiniteCellSizeFactorsOffset();
  const int infiniteCellSizeFactorsBlockSize = iteratorHandler.getInfiniteCellSizeFactorsBlockSize();

  while (iteratorHandler.getRankIterator() == rank) {
    // std::set handles duplicate elements.
    newContainerOptions.emplace(*iteratorHandler.getContainerIterator());
    newCellSizeFactors.emplace(*iteratorHandler.getCellSizeFactorIterator());
    newTraversalOptions.emplace(*iteratorHandler.getTraversalIterator());
    newLoadEstimatorOptions.emplace(*iteratorHandler.getLoadEstimatorIterator());
    newDataLayoutOptions.emplace(*iteratorHandler.getDataLayoutIterator());
    newNewton3Options.emplace(*iteratorHandler.getNewton3Iterator());
    newThreadCounts.emplace(*iteratorHandler.getThreadCountIterator());

    iteratorHandler.advanceIterators(numConfigs, commSize);
  }

  // ============== assigning to local search space ============================

  containerOptions = newContainerOptions;
  if (not cellSizeFactors.isFinite()) {
    const double min = cellSizeFactors.getMin();
    const double max = cellSizeFactors.getMax();
    const double delta = (max - min) / infiniteCellSizeFactorsBlockSize;
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
                              NumberSetFinite<int> &threadCounts, InteractionTypeOption interactionType, const int rank,
                              const int commSize) {
  const auto numConfigs =
      static_cast<int>(getSearchSpaceSize(containerOptions, cellSizeFactors, traversalOptions, loadEstimatorOptions,
                                          dataLayoutOptions, newton3Options, interactionType, threadCounts));

  if (numConfigs == 0) {
    utils::ExceptionHandler::exception("Could not generate valid configurations, aborting");
    return;
  }

  // Creates a set for each option and each rank containing the serialized versions (std::byte or double) of all
  // options assigned to that rank.
  generateDistribution(numConfigs, commSize, rank, containerOptions, cellSizeFactors, traversalOptions,
                       loadEstimatorOptions, dataLayoutOptions, newton3Options, interactionType, threadCounts);

  AutoPasLog(DEBUG,
             "After distributing: {} containers, {} cellSizeFactors, {} traversals, {} dataLayouts, {} newton3s, {} "
             "threadCounts"
             " => {} total configs",
             containerOptions.size(), /*cellSizeFactorsSize*/ (cellSizeFactors.isFinite() ? cellSizeFactors.size() : 1),
             traversalOptions.size(), dataLayoutOptions.size(), newton3Options.size(), threadCounts.size(),
             getSearchSpaceSize(containerOptions, cellSizeFactors, traversalOptions, loadEstimatorOptions,
                                dataLayoutOptions, newton3Options, interactionType, threadCounts));
}

Configuration findGloballyBestConfiguration(AutoPas_MPI_Comm comm, Configuration localOptimalConfig,
                                            long localOptimalTime) {
  // local helper struct to pack a long and int
  int myRank{};
  AutoPas_MPI_Comm_rank(comm, &myRank);
  struct LongIntStruct {
    long value;
    int rank;
  };
  const LongIntStruct localOptimum{localOptimalTime, myRank};
  LongIntStruct globalOptimum;
  // find the rank and value with the minimal value in one collective operation
  AutoPas_MPI_Allreduce(&localOptimum, &globalOptimum, 1, AUTOPAS_MPI_LONG_INT, AUTOPAS_MPI_MINLOC, comm);
  // broadcast the best configuration
  SerializedConfiguration serializedConfiguration = serializeConfiguration(localOptimalConfig);
  AutoPas_MPI_Bcast(serializedConfiguration.data(), serializedConfiguration.size(), AUTOPAS_MPI_BYTE,
                    globalOptimum.rank, comm);

  const Configuration deserializedConfig = deserializeConfiguration(serializedConfiguration);
  AutoPasLog(DEBUG, "Globally optimal configuration: {}", deserializedConfig.toString());

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
  config[5] = castToByte(configuration.interactionType);
  // Doubles can't be easily truncated, so store all 8 bytes via memcpy
  std::memcpy(&config[6], &configuration.cellSizeFactor, sizeof(double));
  std::memcpy(&config[6 + sizeof(double)], &configuration.threadCount, sizeof(int));
  return config;
}

std::vector<std::byte> serializeConfigurations(const std::vector<Configuration> &configurations) {
  // The size of one configuration in bytes as defined by the class type.
  constexpr auto serializedConfSize = std::tuple_size_v<SerializedConfiguration>;

  std::vector<std::byte> confsSerialized{};
  confsSerialized.reserve(configurations.size() * serializedConfSize);
  for (const auto &conf : configurations) {
    const auto serializedConf = serializeConfiguration(conf);
    confsSerialized.insert(confsSerialized.end(), serializedConf.begin(), serializedConf.end());
  }
  return confsSerialized;
}

Configuration deserializeConfiguration(SerializedConfiguration config) {
  double cellSizeFactor{0.};
  std::memcpy(&cellSizeFactor, &config[6], sizeof(double));
  int threadCount(autopas_get_max_threads());
  std::memcpy(&cellSizeFactor, &config[6 + sizeof(double)], sizeof(int));
  return {static_cast<ContainerOption::Value>(config[0]),       cellSizeFactor,
          static_cast<TraversalOption::Value>(config[1]),       static_cast<LoadEstimatorOption::Value>(config[2]),
          static_cast<DataLayoutOption::Value>(config[3]),      static_cast<Newton3Option::Value>(config[4]),
          static_cast<InteractionTypeOption::Value>(config[5]), threadCount};
}

std::vector<Configuration> deserializeConfigurations(const std::vector<std::byte> &configurationsSerialized) {
  // The size of one configuration in bytes as defined by the class type.
  constexpr auto serializedConfSize = std::tuple_size_v<SerializedConfiguration>;

  std::vector<Configuration> configurations{};
  configurations.reserve(configurationsSerialized.size() / serializedConfSize);
  for (size_t i = 0; i < configurationsSerialized.size(); i += serializedConfSize) {
    // copy the bytes of one configuration into a dedicated buffer
    SerializedConfiguration serializedConfig{};
    std::copy(configurationsSerialized.begin() + i, configurationsSerialized.begin() + i + serializedConfSize,
              serializedConfig.begin());
    // turn the byte buffer into a config and store it in the return vector
    configurations.push_back(deserializeConfiguration(serializedConfig));
  }
  return configurations;
}

void distributeRanksInBuckets(AutoPas_MPI_Comm comm, AutoPas_MPI_Comm *bucket, double smoothedPDBinStdDevDensity,
                              double smoothedPDBinMaxDensity, double MPITuningMaxDifferenceForBucket,
                              double MPITuningWeightForMaxDensity) {
  int rank{};
  AutoPas_MPI_Comm_rank(comm, &rank);
  int commSize{};
  AutoPas_MPI_Comm_size(comm, &commSize);

  // If invalid values are passed, every rank gets its own bucket.
  if (smoothedPDBinStdDevDensity < 0 or smoothedPDBinMaxDensity < 0) {
    AutoPas_MPI_Comm_split(comm, rank, rank, bucket);
    return;
  }

  std::vector<double> similarityMetrics(commSize);
  double similarityMetric = smoothedPDBinStdDevDensity + MPITuningWeightForMaxDensity * smoothedPDBinMaxDensity;

  // debug print for evaluation
  AutoPasLog(DEBUG, "similarityMetric           of rank {} is: {}", rank, similarityMetric);
  AutoPasLog(DEBUG, "smoothedPDBinStdDevDensity of rank {} is: {}", rank, smoothedPDBinStdDevDensity);
  AutoPasLog(DEBUG, "smoothedPDBinMaxDensity    of rank {} is: {}", rank, smoothedPDBinMaxDensity);

  // get all the similarityMetrics of the other ranks
  AutoPas_MPI_Allgather(&similarityMetric, 1, AUTOPAS_MPI_DOUBLE, similarityMetrics.data(), 1, AUTOPAS_MPI_DOUBLE,
                        comm);

  // sort all values
  std::sort(similarityMetrics.begin(), similarityMetrics.end());

  // calculate absolute differences between neighbouring values
  std::vector<double> differences;
  std::adjacent_difference(similarityMetrics.begin(), similarityMetrics.end(), std::back_inserter(differences));

  // convert differences to percentage changes
  std::transform(differences.begin(), differences.end(), similarityMetrics.begin(), differences.begin(),
                 std::divides<>());

  int current_bucket = 0;
  int my_bucket = 0;

  for (int i = 0; (size_t)i < similarityMetrics.size(); i++) {
    // if a difference exceeds MPITuningMaxDifferenceForBucket, start a new bucket
    if (differences[i] > MPITuningMaxDifferenceForBucket) current_bucket++;

    // debug print for evaluation
    AutoPasLog(DEBUG, "I am rank {} bucket {} new value {} ", rank, current_bucket, similarityMetrics[i]);
    if (similarityMetrics[i] == similarityMetric) my_bucket = current_bucket;
  }
  // split MPI_Comm in as many new communications as there are groups with similar scenarios
  AutoPas_MPI_Comm_split(comm, my_bucket, rank, bucket);
}

std::vector<Configuration> gatherConfigurations(AutoPas_MPI_Comm comm,
                                                const std::vector<Configuration> &localConfigurations, int root) {
  int rank{};
  AutoPas_MPI_Comm_rank(comm, &rank);
  int numRanks{};
  AutoPas_MPI_Comm_size(comm, &numRanks);

  // turn local vector of configurations into a vector of bits that we can send
  std::vector<std::byte> localConfsSerialized = serializeConfigurations(localConfigurations);

  const auto localNumBytes = static_cast<int>(localConfsSerialized.size());
  std::vector<int> globalNumBytes(numRanks);
  AutoPas_MPI_Gather(&localNumBytes, 1, AUTOPAS_MPI_INT, globalNumBytes.data(), 1, AUTOPAS_MPI_INT, root, comm);

  const auto totalNumBytes = std::reduce(globalNumBytes.begin(), globalNumBytes.end());
  std::vector<std::byte> globalConfsSerialized(totalNumBytes);

  // offsets of each ranks data in the receive-buffer
  std::vector<int> displacements(numRanks);
  for (int i = 1; i < numRanks; ++i) {
    // the displacement is the sum of all previous numBytes
    displacements[i] = displacements[i - 1] + globalNumBytes[i - 1];
  }

  AutoPas_MPI_Gatherv(localConfsSerialized.data(), localNumBytes, AUTOPAS_MPI_BYTE, globalConfsSerialized.data(),
                      globalNumBytes.data(), displacements.data(), AUTOPAS_MPI_BYTE, root, comm);

  return deserializeConfigurations(globalConfsSerialized);
}

}  // namespace autopas::utils::AutoPasConfigurationCommunicator
