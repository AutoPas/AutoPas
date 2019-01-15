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
 * All accepted delimiters to split input strings.
 */
constexpr char delimiters[] = " ,;|/";

/**
 * Splits a string by multiple delimiters.
 * @param searchString
 * @param delimiters
 * @return Vector of substrings.
 */
static std::vector<std::string> tokenize(std::string &searchString, std::string delimiters) {
  std::vector<std::string> wordVector;

  std::size_t prev = 0, pos;
  while ((pos = searchString.find_first_of(delimiters, prev)) != std::string::npos) {
    if (pos > prev) wordVector.push_back(searchString.substr(prev, pos - prev));
    prev = pos + 1;
  }
  if (prev < searchString.length()) wordVector.push_back(searchString.substr(prev, std::string::npos));

  return wordVector;
}

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

  auto words = tokenize(traversalOptionsString, delimiters);

  for (auto &word : words) {
    if (word.find("01") != std::string::npos) {
      if (word.find('v') != std::string::npos)
        traversalOptions.push_back(autopas::TraversalOptions::c01Verlet);
      else
        traversalOptions.push_back(autopas::TraversalOptions::c01);
    } else if (word.find("c08") != std::string::npos) {
      traversalOptions.push_back(autopas::TraversalOptions::c08);
    } else if (word.find("18") != std::string::npos) {
      if (word.find('v') != std::string::npos)
        traversalOptions.push_back(autopas::TraversalOptions::c18Verlet);
      else
        traversalOptions.push_back(autopas::TraversalOptions::c18);
    } else if (word.find("dir") != std::string::npos) {
      traversalOptions.push_back(autopas::TraversalOptions::directSumTraversal);
    } else if (word.find("sli") != std::string::npos) {
      if (word.find('v') != std::string::npos)
        traversalOptions.push_back(autopas::TraversalOptions::slicedVerlet);
      else
        traversalOptions.push_back(autopas::TraversalOptions::sliced);
    }
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

  auto words = tokenize(containerOptionsString, delimiters);

  for (auto &word : words) {
    if (word.find("dir") != std::string::npos or word.find("ds") != std::string::npos) {
      containerOptions.push_back(autopas::ContainerOptions::directSum);
    } else if (word.find("linked") != std::string::npos or word.find("lc") != std::string::npos) {
      containerOptions.push_back(autopas::ContainerOptions::linkedCells);
    } else if (word.find('v') != std::string::npos) {
      if (word.find("cl") != std::string::npos) {
        containerOptions.push_back(autopas::ContainerOptions::verletClusterLists);
      } else if (word.find("cel") != std::string::npos) {
        containerOptions.push_back(autopas::ContainerOptions::verletListsCells);
      } else {
        containerOptions.push_back(autopas::ContainerOptions::verletLists);
      }
    }
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