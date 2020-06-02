/**
 * @file AutoPasConfigurationCommunicator.cpp
 * @author W. Thieme
 * @date 29.05.2020
 */

#include "AutoPasConfigurationCommunicator.h"

namespace autopas {

int getSearchSpaceSize(std::set<ContainerOption> &containerOptions,
                       std::unique_ptr<NumberSet<double>> &cellSizeFactors,
                       std::set<TraversalOption> &traversalOptions,
                       std::set<DataLayoutOption> &dataLayoutOptions,
                       std::set<Newton3Option> &newton3Options) {
  int numConfigs = 0;
  int cellSizeFactorArraySize;
  if (cellSizeFactors->isFinite()) {
    cellSizeFactorArraySize = cellSizeFactors->size();
  } else {
    cellSizeFactorArraySize = 1;
  }
  for (const auto &containerOption : containerOptions) {
    // get all traversals of the container and restrict them to the allowed ones
    const std::set<TraversalOption> &allContainerTraversals =
            compatibleTraversals::allCompatibleTraversals(containerOption);
    std::set<TraversalOption> allowedAndApplicable;
    std::set_intersection(traversalOptions.begin(), traversalOptions.end(),
                          allContainerTraversals.begin(), allContainerTraversals.end(),
                          std::inserter(allowedAndApplicable, allowedAndApplicable.begin()));
    numConfigs += cellSizeFactorArraySize * allowedAndApplicable.size() * dataLayoutOptions.size()
                  * newton3Options.size();
  }
  return numConfigs;
}

inline void advanceIterators(std::set<ContainerOption>::iterator &containerIT,
                             std::set<double>::iterator &cellSizeFactorsIT,
                             const std::set<double> &cellSizeFactors,
                             std::set<TraversalOption>::iterator &traversalIT,
                             const std::set<TraversalOption> &traversalOptions,
                             std::set<DataLayoutOption>::iterator &dataLayoutIT,
                             const std::set<DataLayoutOption> &dataLayouts,
                             std::set<Newton3Option>::iterator &newton3OptionIT,
                             const std::set<Newton3Option> &newton3Options,
                             int &rankIterator, int &remainingBlockSize, int &remainder,
                             int &infiniteCellSizeFactorsOffset, int &infiniteCellSizeFactorsBlockSize,
                             const int numConfigs, const int commSize) {
  if (numConfigs >= commSize or remainingBlockSize == 0) {
    // advance to the next valid config
    ++newton3OptionIT;
    if (newton3OptionIT != newton3Options.end())
      return;
    newton3OptionIT = newton3Options.begin();
    ++dataLayoutIT;
    if (dataLayoutIT != dataLayouts.end())
      return;
    dataLayoutIT = dataLayouts.begin();
    // @todo handle invalid combinations of containers and traversals
    ++traversalIT;
    if (traversalIT != traversalOptions.end())
      return;
    traversalIT = traversalOptions.begin();
    ++cellSizeFactorsIT;
    if (cellSizeFactorsIT != cellSizeFactors.end())
      return;
    cellSizeFactorsIT = cellSizeFactors.begin();
    ++containerIT;
  }

  if (commSize >= numConfigs or remainingBlockSize == 0) {
    // advance to the next rank
    ++rankIterator;

    ++infiniteCellSizeFactorsOffset;
  }

  if (remainingBlockSize == 0) {
    if (numConfigs >= commSize) {
      remainingBlockSize = numConfigs / commSize;
    } else {
      remainingBlockSize = commSize / numConfigs;

      infiniteCellSizeFactorsBlockSize = remainingBlockSize;
      infiniteCellSizeFactorsOffset = 0;
    }
    if (remainder > 0) {
      ++remainingBlockSize;
      --remainder;
    }

  }
}

void generateDistribution(const int numConfigs, const int commSize, const int rank,
                          std::set<ContainerOption> &containerOptions,
                          std::unique_ptr<autopas::NumberSet<double>> &cellSizeFactors,
                          std::set<TraversalOption> &traversalOptions, std::set<DataLayoutOption> &dataLayoutOptions,
                          std::set<Newton3Option> &newton3Options) {

  // ============== setup ======================================================

  auto newContainerOptions = std::set<ContainerOption>();
  auto newCellSizeFactors = std::set<double>();
  auto newTraversalOptions = std::set<TraversalOption>();
  auto newDataLayoutOptions = std::set<DataLayoutOption>();
  auto newNewton3Options = std::set<Newton3Option>();
  // If there are more configurations than ranks, the ranks partition the configurations.
  // If there are more ranks than configurations, the configurations partition the ranks.

  int remainingBlockSize = 0;
  int rankIterator = 0;
  int remainder;
  if (commSize >= numConfigs) {
    remainder = commSize % numConfigs;
  } else {
    remainder = numConfigs % commSize;
  }

  std::set<double> finiteCellSizeFactors;
  if (cellSizeFactors->isFinite()) {
    finiteCellSizeFactors = cellSizeFactors->getAll();
  } else {
    // Dummy value which makes the code simpler in case the cellSizeFactors are not a finite set.
    finiteCellSizeFactors = std::set<double>{-1};
  }

  auto containerIt = containerOptions.begin();
  auto cellSizeIt = finiteCellSizeFactors.begin();
  auto traversalIt = traversalOptions.begin();
  auto dataLayoutIt = dataLayoutOptions.begin();
  auto newton3It = newton3Options.begin();

  int infiniteCellSizeFactorsOffset = 0;
  int infiniteCellSizeFactorsBlockSize = 1;

  // ============== main computation ===========================================

  while (rankIterator < rank) {
    advanceIterators(containerIt, cellSizeIt, finiteCellSizeFactors, traversalIt, traversalOptions, dataLayoutIt,
                     dataLayoutOptions, newton3It, newton3Options, rankIterator, remainingBlockSize, remainder,
                     infiniteCellSizeFactorsOffset, infiniteCellSizeFactorsBlockSize, numConfigs, commSize);
  }

  while (rankIterator == rank) {
    // std::set handles duplicate elements
    newContainerOptions.emplace(*containerIt);
    newCellSizeFactors.emplace(*cellSizeIt);
    newTraversalOptions.emplace(*traversalIt);
    newDataLayoutOptions.emplace(*dataLayoutIt);
    newNewton3Options.emplace(*newton3It);

    int dummy;
    advanceIterators(containerIt, cellSizeIt, finiteCellSizeFactors, traversalIt, traversalOptions, dataLayoutIt,
                     dataLayoutOptions, newton3It, newton3Options, rankIterator, remainingBlockSize, remainder,
                     dummy, dummy, numConfigs, commSize);
  }


  // ============== assigning to local search space ============================

  containerOptions = newContainerOptions;
  if (not cellSizeFactors->isFinite()) {
    // newCellSizeFactors == {-1}
    double min = cellSizeFactors->getMin();
    double max = cellSizeFactors->getMax();
    double delta = (max - min)/infiniteCellSizeFactorsBlockSize;
    cellSizeFactors = std::make_unique<NumberInterval<double>>(
                                       min + delta*infiniteCellSizeFactorsOffset,
                                       min + delta*(infiniteCellSizeFactorsOffset+1));
  } else {
    cellSizeFactors = std::make_unique<NumberSetFinite<double>>(newCellSizeFactors);
  }
  traversalOptions = newTraversalOptions;
  dataLayoutOptions = newDataLayoutOptions;
  newton3Options = newNewton3Options;
}

void AutoPasConfigurationCommunicator::distributeConfigurations(
        std::set<ContainerOption> &containerOptions, std::unique_ptr<NumberSet<double>> &cellSizeFactors,
        std::set<TraversalOption> &traversalOptions, std::set<DataLayoutOption> &dataLayoutOptions,
        std::set<Newton3Option> &newton3Options, AutoPas_MPI_Comm comm) {
  int rank, commSize;
  AutoPas_MPI_Comm_rank(comm, &rank);
  AutoPas_MPI_Comm_size(comm, &commSize);
  int numConfigs = getSearchSpaceSize(containerOptions, cellSizeFactors, traversalOptions, dataLayoutOptions,
                                      newton3Options);

  // Creates a set for each option and each rank containing the serialized versions (std::byte or double) of all
  // options assigned to that rank.
  generateDistribution(numConfigs, commSize, rank, containerOptions, cellSizeFactors, traversalOptions,
                       dataLayoutOptions, newton3Options);

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

  AutoPas_MPI_Bcast(&serializedConfiguration, sizeof(serializedConfiguration), AUTOPAS_MPI_BYTE, optimalRankOut,
                    comm);

  return deserializeConfig(serializedConfiguration);
}


} // namespace autopas
