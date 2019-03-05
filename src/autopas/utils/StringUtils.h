/**
 * @file StringUtils.h
 * @author F. Gratl
 * @date 1/14/19
 */

#pragma once

#include <string>
#include <vector>
#include "autopas/options/ContainerOptions.h"
#include "autopas/options/DataLayoutOptions.h"
#include "autopas/options/SelectorStrategies.h"
#include "autopas/options/TraversalOptions.h"

namespace autopas {
namespace utils {
/**
 * Some functions to parse enums from (input-) strings.
 */
namespace StringUtils {

/**
 * Converts a SelectorStrategy to its respective string representation.
 * @param option
 * @return The string representation or "Unknown option (<IntValue>)".
 */
inline std::string to_string(autopas::SelectorStrategy option) {
  switch (option) {
    case autopas::SelectorStrategy::fastestAbs: {
      return "Fastest-Absolute-Value";
    }
    case autopas::SelectorStrategy::fastestMean: {
      return "Fastest-Mean-Value";
    }
    case autopas::SelectorStrategy::fastestMedian: {
      return "Fastest-Median-Value";
    }
  }
  // do not implement default case to provoke compiler warnings if new options are introduced.
  return "Unknown option (" + std::to_string(option) + ")";
}

/**
 * Converts a DataLayoutOption to its respective string representation.
 * @param option
 * @return The string representation or "Unknown option (<IntValue>)".
 */
inline std::string to_string(DataLayoutOption option) {
  switch (option) {
    case autopas::DataLayoutOption::aos: {
      return "Array-of-Structures";
    }
    case autopas::DataLayoutOption::soa: {
      return "Structure-of-Arrays";
    }
    case autopas::DataLayoutOption::cuda: {
      return "Structure-of-Arrays on Cuda capable device";
    }
  }
  // do not implement default case to provoke compiler warnings if new options are introduced.
  return "Unknown option (" + std::to_string(option) + ")";
}

/**
 * Converts a ContainerOption to its respective string representation.
 * @param option
 * @return The string representation or "Unknown option (<IntValue>)".
 */
inline std::string to_string(ContainerOptions option) {
  switch (option) {
    case autopas::ContainerOptions::directSum: {
      return "DirectSum";
    }
    case autopas::ContainerOptions::linkedCells: {
      return "LinkedCells";
    }
    case autopas::ContainerOptions::verletLists: {
      return "VerletLists";
    }
    case autopas::ContainerOptions::verletListsCells: {
      return "VerletListsCells";
    }
    case autopas::ContainerOptions::verletClusterLists: {
      return "VerletClusterLists";
    }
  }
  // do not implement default case to provoke compiler warnings if new options are introduced.
  return "Unknown option (" + std::to_string(option) + ")";
}

/**
 * Converts a TraversalOption to its respective string representation.
 * @param option
 * @return The string representation or "Unknown option (<IntValue>)".
 */
inline std::string to_string(TraversalOptions option) {
  switch (option) {
    case autopas::TraversalOptions::dummyTraversal: {
      return "dummyTraversal";
    }
    case autopas::TraversalOptions::c01: {
      return "c01";
    }
    case autopas::TraversalOptions::c08: {
      return "c08";
    }
    case autopas::TraversalOptions::c18: {
      return "c18";
    }
    case autopas::TraversalOptions::sliced: {
      return "sliced";
    }
    case autopas::TraversalOptions::directSumTraversal: {
      return "directSum";
    }
    case autopas::TraversalOptions::c01Verlet: {
      return "verlet-c01";
    }
    case autopas::TraversalOptions::c18Verlet: {
      return "verlet-c18";
    }
    case autopas::TraversalOptions::slicedVerlet: {
      return "verlet-sliced";
    }
  }
  // do not implement default case to provoke compiler warnings if new options are introduced.
  return "Unknown option (" + std::to_string(option) + ")";
}

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
inline std::vector<std::string> tokenize(const std::string &searchString, const std::string &delimiters) {
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
 * Allowed delimiters can be found in autopas::utils::StringUtils::delimiters
 *
 * Possible options: c01, c08, c18, direct, sliced, verlet01, verlet18, verlet-sliced
 *
 * @param traversalOptionsString String containing traversal options.
 * @param ignoreUnknownOptions If set to false, a 'autopas::TraversalOptions(-1)' will be inserted in the return vector
 * for each not parsable word.
 * @return Vector of TraversalOption enums. If no valid option was found the empty vector is returned.
 */
inline std::vector<autopas::TraversalOptions> parseTraversalOptions(const std::string &traversalOptionsString,
                                                                    bool ignoreUnknownOptions = true) {
  std::vector<autopas::TraversalOptions> traversalOptions;

  auto words = tokenize(traversalOptionsString, delimiters);

  for (auto &word : words) {
    if (word.find("01") != std::string::npos) {
      if (word.find('v') != std::string::npos)
        traversalOptions.emplace_back(autopas::TraversalOptions::c01Verlet);
      else
        traversalOptions.emplace_back(autopas::TraversalOptions::c01);
    } else if (word.find("c08") != std::string::npos) {
      traversalOptions.emplace_back(autopas::TraversalOptions::c08);
    } else if (word.find("18") != std::string::npos) {
      if (word.find('v') != std::string::npos)
        traversalOptions.emplace_back(autopas::TraversalOptions::c18Verlet);
      else
        traversalOptions.emplace_back(autopas::TraversalOptions::c18);
    } else if (word.find("dir") != std::string::npos) {
      traversalOptions.emplace_back(autopas::TraversalOptions::directSumTraversal);
    } else if (word.find("sli") != std::string::npos) {
      if (word.find('v') != std::string::npos)
        traversalOptions.emplace_back(autopas::TraversalOptions::slicedVerlet);
      else
        traversalOptions.emplace_back(autopas::TraversalOptions::sliced);
    } else if (not ignoreUnknownOptions) {
      traversalOptions.emplace_back(autopas::TraversalOptions(-1));
    }
  }
  return traversalOptions;
}

/**
 * Converts a string of options to a vector of enums. The options are expected to be lower case.
 * Allowed delimiters can be found in autopas::utils::StringUtils::delimiters
 *
 * Possible options: directSum, linkedCells, verletLists, vcells, vcluster
 *
 * @param containerOptionsString String containing container options.
 * @param ignoreUnknownOptions If set to false, a 'autopas::ContainerOptions(-1)' will be inserted in the return vector
 * for each not parsable word.
 * @return Vector of ContainerOption enums. If no valid option was found the empty vector is returned.
 */
inline std::vector<autopas::ContainerOptions> parseContainerOptions(const std::string &containerOptionsString,
                                                                    bool ignoreUnknownOptions = true) {
  std::vector<autopas::ContainerOptions> containerOptions;

  auto words = tokenize(containerOptionsString, delimiters);

  for (auto &word : words) {
    if (word.find("dir") != std::string::npos or word.find("ds") != std::string::npos) {
      containerOptions.emplace_back(autopas::ContainerOptions::directSum);
    } else if (word.find("linked") != std::string::npos or word.find("lc") != std::string::npos) {
      containerOptions.emplace_back(autopas::ContainerOptions::linkedCells);
    } else if (word.find('v') != std::string::npos) {
      if (word.find("cl") != std::string::npos) {
        containerOptions.emplace_back(autopas::ContainerOptions::verletClusterLists);
      } else if (word.find("cel") != std::string::npos) {
        containerOptions.emplace_back(autopas::ContainerOptions::verletListsCells);
      } else {
        containerOptions.emplace_back(autopas::ContainerOptions::verletLists);
      }
    } else if (not ignoreUnknownOptions) {
      containerOptions.emplace_back(autopas::ContainerOptions(-1));
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
inline autopas::SelectorStrategy parseSelectorStrategy(const std::string &selectorStrategyString) {
  // hack to initialize the enum out of range as an error value.
  auto selectorStrategy(autopas::SelectorStrategy(-1));
  if (selectorStrategyString.find("abs") != std::string::npos or
      selectorStrategyString.find("min") != std::string::npos) {
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
 * @param dataLayoutString String containing the data layout option.
 * @return An enum representing the data layout. If no valid option was found 'autopas::DataLayoutOption(-1)' is
 * returned.
 */
inline autopas::DataLayoutOption parseDataLayout(const std::string &dataLayoutString) {
  // hack to initialize the enum out of range as an error value.
  auto dataLayout(autopas::DataLayoutOption(-1));
  if (dataLayoutString.find("aos") != std::string::npos or
      dataLayoutString.find("array-of-struct") != std::string::npos) {
    dataLayout = autopas::DataLayoutOption::aos;
  } else if (dataLayoutString.find("soa") != std::string::npos or
             dataLayoutString.find("-of-array") != std::string::npos) {
    dataLayout = autopas::DataLayoutOption::soa;
  } else if (dataLayoutString.find("cuda") != std::string::npos) {
    dataLayout = autopas::DataLayoutOption::cuda;
  }
  return dataLayout;
}
}  // namespace StringUtils
}  // namespace utils
}  // namespace autopas
