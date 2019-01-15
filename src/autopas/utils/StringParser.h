/**
 * @file StringParser.h
 * @author F. Gratl
 * @date 1/14/19
 */

#pragma once

#include <string>
#include <vector>
#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/containers/cellPairTraversals/TraversalInterface.h"
#include "autopas/selectors/AutoTuner.h"

namespace autopas {
namespace utils {
/**
 * Some functions to parse enums from (input-) strings.
 */
namespace StringParser {

/**
 * Converts a string of options to a vector of enums. The options are expected to be lower case.
 *
 * Possible options: c01, c08, c18, direct, sliced, v01, v18, vsl
 *
 * @param traversalOptionsString String containing traversal options.
 * @return Vector of TraversalOption enums. If no valid option was found the empty vector is returned.
 */
static std::vector<autopas::TraversalOptions> parseTraversalOptions(std::string &traversalOptionsString) {
  std::vector<autopas::TraversalOptions> traversalOptions;

  // @TODO: rewrite such that sli and slicedVerlet don't crash
  if (traversalOptionsString.find("c01") != std::string::npos) {
    traversalOptions.push_back(autopas::TraversalOptions::c01);
  }
  if (traversalOptionsString.find("c08") != std::string::npos) {
    traversalOptions.push_back(autopas::TraversalOptions::c08);
  }
  if (traversalOptionsString.find("c18") != std::string::npos) {
    traversalOptions.push_back(autopas::TraversalOptions::c18);
  }
  if (traversalOptionsString.find("dir") != std::string::npos) {
    traversalOptions.push_back(autopas::TraversalOptions::directSumTraversal);
  }
  if (traversalOptionsString.find("sli") != std::string::npos) {
    traversalOptions.push_back(autopas::TraversalOptions::sliced);
  }
  if (traversalOptionsString.find("v01") != std::string::npos) {
    traversalOptions.push_back(autopas::TraversalOptions::c01Verlet);
  }
  if (traversalOptionsString.find("v18") != std::string::npos) {
    traversalOptions.push_back(autopas::TraversalOptions::c18Verlet);
  }
  if (traversalOptionsString.find("vsl") != std::string::npos) {
    traversalOptions.push_back(autopas::TraversalOptions::slicedVerlet);
  }
  return traversalOptions;
}

/**
 * Converts a string of options to a vector of enums. The options are expected to be lower case.
 *
 * Possible options: directSum, linkedCells, verletLists, vcells, vcluster
 *
 * @param containerOptionsString String containing container options.
 * @return Vector of ContainerOption enums. If no valid option was found the empty vector is returned.
 */
static std::vector<autopas::ContainerOptions> parseContainerOptions(std::string &containerOptionsString) {
  std::vector<autopas::ContainerOptions> containerOptions;

  // @TODO: rewrite such that verlet stuff doesn't crash
  if (containerOptionsString.find("direct") != std::string::npos or
      containerOptionsString.find("ds") != std::string::npos) {
    containerOptions.push_back(autopas::directSum);
  }
  if (containerOptionsString.find("linked") != std::string::npos or
      containerOptionsString.find("lc") != std::string::npos) {
    containerOptions.push_back(autopas::linkedCells);
  }
  if (containerOptionsString.find("verlet") != std::string::npos or
      containerOptionsString.find("vl") != std::string::npos) {
    containerOptions.push_back(autopas::verletLists);
  }
  if (containerOptionsString.find("vcell") != std::string::npos) {
    containerOptions.push_back(autopas::verletListsCells);
  }
  if (containerOptionsString.find("vcluster") != std::string::npos) {
    containerOptions.push_back(autopas::verletClusterLists);
  }
  return containerOptions;
}

/**
 * Converts a string containing a selector strategy to an enum. The option is expected to be lower case.
 *
 * Possible options: absolute, mean, median
 *
 * @param selectorStrategyString String containing the selector option.
 * @return An enum representing the selector Strategy. If no valid option was found 'autopas::SelectorStrategy(-1)' is
 * returned.
 */
static autopas::SelectorStrategy parseSelectorStrategy(std::string &selectorStrategyString) {
  // hack to initialize the enum out of range as an error value.
  auto selectorStrategy(autopas::SelectorStrategy(-1));
  if (selectorStrategyString.find("abs") != std::string::npos) {
    selectorStrategy = autopas::SelectorStrategy::fastestAbs;
  } else if (selectorStrategyString.find("mea") != std::string::npos) {
    selectorStrategy = autopas::SelectorStrategy::fastestMean;
  } else if (selectorStrategyString.find("med") != std::string::npos) {
    selectorStrategy = autopas::SelectorStrategy::fastestMedian;
  }
  return selectorStrategy;
}

/**
 * Converts a string containing a data layout to an enum. The option is expected to be lower case.
 *
 * Possible options: aos, soa
 *
 * @param dataLayoutSting String containing the data layout option.
 * @return An enum representing the data layout. If no valid option was found 'autopas::DataLayoutOption(-1)' is
 * returned.
 */
static autopas::DataLayoutOption parseDataLayout(std::string &dataLayoutSting) {
  // hack to initialize the enum out of range as an error value.
  auto dataLayout(autopas::DataLayoutOption(-1));
  if (dataLayoutSting.find("aos") != std::string::npos) {
    dataLayout = autopas::DataLayoutOption::aos;
  } else if (dataLayoutSting.find("soa") != std::string::npos) {
    dataLayout = autopas::DataLayoutOption::soa;
  }
  return dataLayout;
}
}  // namespace StringParser
}  // namespace utils
}  // namespace autopas