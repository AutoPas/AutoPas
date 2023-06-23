/**
 * @file SearchSpaceGenerators.cpp
 * @author F. Gratl
 * @date 23.06.23
 */

#include "SearchSpaceGenerators.h"

#include <map>
#include <set>
#include <vector>

#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/tuning/utils/Evidence.h"

std::map<autopas::Configuration, std::vector<autopas::Evidence>> autopas::SearchSpaceGenerators::optionCrossProduct(
    const std::set<ContainerOption> &allowedContainerOptions, const std::set<double> &allowedCellSizeFactors,
    const std::set<TraversalOption> &allowedTraversalOptions,
    const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
    const std::set<DataLayoutOption> &allowedDataLayoutOptions, const std::set<Newton3Option> &allowedNewton3Options) {
  std::map<autopas::Configuration, std::vector<autopas::Evidence>> searchSpace;

  // generate all potential configs
  for (const auto &containerOption : allowedContainerOptions) {
    // get all traversals of the container and restrict them to the allowed ones
    const std::set<TraversalOption> &allContainerTraversals =
        compatibleTraversals::allCompatibleTraversals(containerOption);
    std::set<TraversalOption> allowedAndApplicable;
    std::set_intersection(allowedTraversalOptions.begin(), allowedTraversalOptions.end(),
                          allContainerTraversals.begin(), allContainerTraversals.end(),
                          std::inserter(allowedAndApplicable, allowedAndApplicable.begin()));

    for (const auto &cellSizeFactor : allowedCellSizeFactors)
      for (const auto &traversalOption : allowedAndApplicable) {
        // if load estimators are not applicable LoadEstimatorOption::none is returned.
        const std::set<LoadEstimatorOption> allowedAndApplicableLoadEstimators =
            loadEstimators::getApplicableLoadEstimators(containerOption, traversalOption, allowedLoadEstimatorOptions);
        for (const auto &loadEstimatorOption : allowedAndApplicableLoadEstimators) {
          for (const auto &dataLayoutOption : allowedDataLayoutOptions) {
            for (const auto &newton3Option : allowedNewton3Options) {
              const Configuration configuration{containerOption,     cellSizeFactor,   traversalOption,
                                                loadEstimatorOption, dataLayoutOption, newton3Option};
              if (configuration.isValid()) {
                searchSpace.emplace(configuration, std::vector<autopas::Evidence>{});
              }
            }
          }
        }
      }
  }

  if (searchSpace.empty()) {
    autopas::utils::ExceptionHandler::exception("No valid configurations could be created.");
  }
  return searchSpace;
}
