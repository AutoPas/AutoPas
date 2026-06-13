/**
 * @file SearchSpaceGenerators.cpp
 * @author F. Gratl
 * @date 23.06.23
 */

#include "autopas/tuning/utils/SearchSpaceGenerators.h"

#include <memory>
#include <set>

#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/utils/NumberSet.h"

namespace autopas {

std::set<Configuration> SearchSpaceGenerators::cartesianProduct(
    const std::set<ContainerOption> &allowedContainerOptions, const std::set<TraversalOption> &allowedTraversalOptions,
    const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
    const std::set<DataLayoutOption> &allowedDataLayoutOptions, const std::set<Newton3Option> &allowedNewton3Options,
    const NumberSet<double> *allowedCellSizeFactors,
    const std::set<VectorizationPatternOption> &allowedVecPatternOptions,
    const InteractionTypeOption &interactionType) {
  if (allowedCellSizeFactors->isInterval()) {
    utils::ExceptionHandler::exception("Cross product does not work with continuous cell size factors!");
  }
  const auto cellSizeFactors = allowedCellSizeFactors->getAll();

  std::set<Configuration> searchSet;
  // generate all potential configs
  for (const auto &containerOption : allowedContainerOptions) {
    // get all traversals of the container and restrict them to the allowed ones
    const std::set<TraversalOption> &allContainerTraversals =
        compatibleTraversals::allCompatibleTraversals(containerOption, interactionType);
    std::set<TraversalOption> allowedAndApplicable;
    std::set_intersection(allowedTraversalOptions.begin(), allowedTraversalOptions.end(),
                          allContainerTraversals.begin(), allContainerTraversals.end(),
                          std::inserter(allowedAndApplicable, allowedAndApplicable.begin()));

    for (const auto &traversalOption : allowedAndApplicable) {
      // if load estimators are not applicable LoadEstimatorOption::none is returned.
      const std::set<LoadEstimatorOption> allowedAndApplicableLoadEstimators =
          loadEstimators::getApplicableLoadEstimators(containerOption, traversalOption, allowedLoadEstimatorOptions);
      for (const auto csf : cellSizeFactors) {
        for (const auto &loadEstimatorOption : allowedAndApplicableLoadEstimators) {
          for (const auto &dataLayoutOption : allowedDataLayoutOptions) {
            for (const auto &newton3Option : allowedNewton3Options) {
              for (const auto &vecPatternOption : allowedVecPatternOptions) {
                const Configuration configuration{containerOption,  csf,           traversalOption, loadEstimatorOption,
                                                  dataLayoutOption, newton3Option, interactionType, vecPatternOption};
                if (configuration.hasCompatibleValues()) {
                  searchSet.insert(configuration);
                }
              }
            }
          }
        }
      }
    }
  }

  if (searchSet.empty()) {
    utils::ExceptionHandler::exception("No valid configurations could be created.");
  }
  return searchSet;
}

SearchSpaceGenerators::OptionSpace SearchSpaceGenerators::inferOptionDimensions(
    const std::set<Configuration> &searchSet) {
  OptionSpace optionSpace;
  for (const auto &[container, traversal, vecPattern, loadEst, dataLayout, newton3, csf, interactT] : searchSet) {
    optionSpace.containerOptions.insert(container);
    optionSpace.traversalOptions.insert(traversal);
    optionSpace.loadEstimatorOptions.insert(loadEst);
    optionSpace.dataLayoutOptions.insert(dataLayout);
    optionSpace.newton3Options.insert(newton3);
    optionSpace.cellSizeFactors.insert(csf);
    optionSpace.vecPatternOptions.insert(vecPattern);
  }
  return optionSpace;
}

std::set<double> SearchSpaceGenerators::calculateRelevantCsfs(const NumberInterval<double> &numberInterval,
                                                              double interactionLength, double domainLengthX) {
  // helper function for readability
  auto calcCsf = [](double domainLength, double interactionLength, unsigned int numberOfCells) {
    return domainLength / (static_cast<double>(numberOfCells) * interactionLength);
  };

  // ceil to stay within the given interval
  const auto numCellsMin =
      static_cast<unsigned int>(std::ceil(domainLengthX / (interactionLength * numberInterval.getMax())));
  // floor to stay within the given interval
  const auto numCellsMax =
      static_cast<unsigned int>(std::floor(domainLengthX / (interactionLength * numberInterval.getMin())));

  // find the CSFs for all possible amounts of cells in the interval
  std::set<double> relevantCsfs{};
  for (auto numCells = numCellsMin; numCells <= numCellsMax; ++numCells) {
    relevantCsfs.insert(calcCsf(domainLengthX, interactionLength, numCells));
  }

  return relevantCsfs;
}

SearchSpaceGenerators::OptionSpace::OptionSpace() = default;

SearchSpaceGenerators::OptionSpace::~OptionSpace() noexcept = default;
}  // namespace autopas
