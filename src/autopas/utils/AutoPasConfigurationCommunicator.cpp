/**
 * @file AutoPasConfigurationCommunicator.cpp
 * @author W. Thieme
 * @date 29.05.2020
 */

#include "AutoPasConfigurationCommunicator.h"

namespace autopas::AutoPasConfigurationCommunicator {

int getSearchSpaceSize(std::set<ContainerOption> &containerOptions,
                       std::set<double> &cellSizeFactors,
                       std::set<TraversalOption> &traversalOptions,
                       std::set<DataLayoutOption> &dataLayoutOptions,
                       std::set<Newton3Option> &newton3Options) {
  int numConfigs = 0;
  for (const auto &containerOption : containerOptions) {
    // get all traversals of the container and restrict them to the allowed ones
    const std::set<TraversalOption> &allContainerTraversals =
            compatibleTraversals::allCompatibleTraversals(containerOption);
    std::set<TraversalOption> allowedAndApplicable;
    std::set_intersection(traversalOptions.begin(), traversalOptions.end(),
                          allContainerTraversals.begin(), allContainerTraversals.end(),
                          std::inserter(allowedAndApplicable, allowedAndApplicable.begin()));
    numConfigs += cellSizeFactors.size() * allowedAndApplicable.size() * dataLayoutOptions.size()
                  * newton3Options.size();
  }
  return numConfigs;
}

inline void advanceIterators(std::set<ContainerOption>::iterator &containerIT,
                             std::set<double>::iterator &cellSizeFactorsIT,
                             const std::set<double>::iterator cellSizeFactorsStart,
                             const std::set<double>::iterator cellSizeFactorsEnd,
                             std::set<TraversalOption>::iterator &traversalIT,
                             const std::set<TraversalOption>::iterator traversalStart,
                             const std::set<TraversalOption>::iterator traversalEnd,
                             std::set<DataLayoutOption>::iterator &dataLayoutIT,
                             const std::set<DataLayoutOption>::iterator dataLayoutStart,
                             const std::set<DataLayoutOption>::iterator dataLayoutEnd,
                             std::set<Newton3Option>::iterator &newton3OptionIT,
                             const std::set<Newton3Option>::iterator newton3OptionStart,
                             const std::set<Newton3Option>::iterator newton3OptionEnd) {
  ++newton3OptionIT;
  if (newton3OptionIT != newton3OptionEnd)
    return;
  newton3OptionIT = newton3OptionStart;
  ++dataLayoutIT;
  if (dataLayoutIT != dataLayoutEnd)
    return;
  dataLayoutIT = dataLayoutStart;
  ++traversalIT;
  if (traversalIT != traversalEnd)
    return;
  traversalIT = traversalStart;
  ++cellSizeFactorsIT;
  if (cellSizeFactorsIT != cellSizeFactorsEnd)
    return;
  cellSizeFactorsIT = cellSizeFactorsStart;
  ++containerIT;
}

void generateDistribution(const int numConfigs, const int commSize, const std::set<ContainerOption> &containerOptions,
                          const std::set<double> &cellSizeFactors, const std::set<TraversalOption> &traversalOptions,
                          const std::set<DataLayoutOption> &dataLayoutOptions,
                          const std::set<Newton3Option> &newton3Options,
                          OptionsPerRank<std::byte> &containerOptionsPerRank,
                          OptionsPerRank<double> &cellSizeFactorsPerRank,
                          OptionsPerRank<std::byte> &traversalOptionsPerRank,
                          OptionsPerRank<std::byte> &dataLayoutOptionsPerRank,
                          OptionsPerRank<std::byte> &newton3OptionsPerRank) {
  // If there are more configurations than ranks, the ranks partition the configurations.
  // If there are more ranks than configurations, the configurations partition the ranks.
  if (commSize <= numConfigs) {
    int configsLeftForRank = 0;
    int remainder = numConfigs % commSize;
    // holds the value of the rank that currently receives new configurations.
    // Start at -1 to be 0 in the first iteration
    int rankIndex = -1;
    // Hold the first option that was assigned to each rank;
    std::set<ContainerOption>::iterator firstContainer;
    std::set<double>::iterator firstCellSizeFactor;
    std::set<TraversalOption>::iterator firstTraversalOption;
    std::set<DataLayoutOption>::iterator firstDataLayout;
    std::set<Newton3Option>::iterator firstNewton3Option;

    for (const auto &container : containerOptions) {
      // get all traversals of the container and restrict them to the allowed ones
      const std::set<TraversalOption> &allContainerTraversals =
              compatibleTraversals::allCompatibleTraversals(container);
      std::set<TraversalOption> allowedAndApplicable;
      std::set_intersection(traversalOptions.begin(), traversalOptions.end(),
                            allContainerTraversals.begin(), allContainerTraversals.end(),
                            std::inserter(allowedAndApplicable, allowedAndApplicable.begin()));
      for (const auto &cellSizeFactor : cellSizeFactors)
        for (const auto &traversal : traversalOptions) {
          for (const auto &dataLayout : dataLayoutOptions) {
            for (const auto &newton3Option : newton3Options) {

              // Advance the rank index whenever one block of configurations is finished
              // Reset new.. and first.. values for the new rankIndex
              if (configsLeftForRank == 0) {
                ++rankIndex;
                configsLeftForRank = numConfigs / commSize;

                // Assign all remaining configs to the first remainder-many ranks
                if (remainder > 0) {
                  ++configsLeftForRank;
                  --remainder;
                }
              }

              containerOptionsPerRank[rankIndex].emplace(castToByte(container));
              cellSizeFactorsPerRank[rankIndex].emplace(cellSizeFactor);
              traversalOptionsPerRank[rankIndex].emplace(castToByte(traversal));
              dataLayoutOptionsPerRank[rankIndex].emplace(castToByte(dataLayout));
              newton3OptionsPerRank[rankIndex].emplace(castToByte(newton3Option));

              --configsLeftForRank;
            }
          }
        }
    }

  } else {
    int ranksLeftForConfig = 0;
    int remainder = commSize % numConfigs;

    auto containerIT = containerOptions.begin();
    auto cellSizeFactorIT = cellSizeFactors.begin();
    auto traversalIT = traversalOptions.begin();
    auto dataLayoutIT = dataLayoutOptions.begin();
    auto newton3OptionIT = newton3Options.begin();

    for (int rankIndex = 0; rankIndex < commSize; ++rankIndex) {
      if (ranksLeftForConfig == 0) {
        advanceIterators(containerIT,
                         cellSizeFactorIT, cellSizeFactors.begin(), cellSizeFactors.end(),
                         traversalIT, traversalOptions.begin(), traversalOptions.end(),
                         dataLayoutIT, dataLayoutOptions.begin(), dataLayoutOptions.end(),
                         newton3OptionIT, newton3Options.begin(), newton3Options.end());
        ranksLeftForConfig = commSize / numConfigs;
        if (remainder > 0) {
          ++ranksLeftForConfig;
          --remainder;
        }
      }

      containerOptionsPerRank[rankIndex].emplace(castToByte(*containerIT));
      cellSizeFactorsPerRank[rankIndex].emplace(*cellSizeFactorIT);
      traversalOptionsPerRank[rankIndex].emplace(castToByte(*traversalIT));
      dataLayoutOptionsPerRank[rankIndex].emplace(castToByte(*dataLayoutIT));
      newton3OptionsPerRank[rankIndex].emplace(castToByte(*newton3OptionIT));

      --ranksLeftForConfig;
    }
  }
}

void distributeConfigurations(std::set<ContainerOption> &containerOptions, std::set<double> &cellSizeFactors,
                              std::set<TraversalOption> &traversalOptions,
                              std::set<DataLayoutOption> &dataLayoutOptions, std::set<Newton3Option> &newton3Options,
                              AutoPas_MPI_Comm comm) {
  int rank, commSize;
  AutoPas_MPI_Comm_rank(comm, &rank);
  AutoPas_MPI_Comm_size(comm, &commSize);
  int numConfigs = getSearchSpaceSize(containerOptions, cellSizeFactors, traversalOptions, dataLayoutOptions,
                                      newton3Options);

  // Creates a set for each option and each rank containing the serialized versions (std::byte or double) of all
  // options assigned to that rank.
  auto containerOptionsPerRank = std::make_unique<std::set<std::byte>[]>(commSize);
  auto cellSizeFactorsPerRank = std::make_unique<std::set<double>[]>(commSize);
  auto traversalOptionsPerRank = std::make_unique<std::set<std::byte>[]>(commSize);
  auto dataLayoutOptionsPerRank = std::make_unique<std::set<std::byte>[]>(commSize);
  auto newton3OptionsPerRank = std::make_unique<std::set<std::byte>[]>(commSize);
  generateDistribution(numConfigs, commSize, containerOptions, cellSizeFactors, traversalOptions, dataLayoutOptions,
                       newton3Options, containerOptionsPerRank, cellSizeFactorsPerRank, traversalOptionsPerRank,
                       dataLayoutOptionsPerRank, newton3OptionsPerRank);

  distributeOption(comm, commSize, rank, containerOptions, containerOptionsPerRank);
  distributeOption(comm, commSize, rank, cellSizeFactors, cellSizeFactorsPerRank);
  distributeOption(comm, commSize, rank, traversalOptions, traversalOptionsPerRank);
  distributeOption(comm, commSize, rank, dataLayoutOptions, dataLayoutOptionsPerRank);
  distributeOption(comm, commSize, rank, newton3Options, newton3OptionsPerRank);
}

} // namespace autopas
