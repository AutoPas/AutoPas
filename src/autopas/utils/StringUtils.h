/**
 * @file StringUtils.h
 * @author F. Gratl
 * @date 1/14/19
 */

#pragma once

#include <autopas/options/TuningStrategyOption.h>
#include <cmath>
#include <regex>
#include <set>
#include <string>
#include <vector>
#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/SelectorStrategyOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/utils/NumberSet.h"

namespace autopas {
namespace utils {
/**
 * Some functions to parse enums from (input-) strings.
 */
namespace StringUtils {

/**
 * Converts a Newton3Option to its respective string representation.
 * @param option
 * @return The string representation or "Unknown option (<IntValue>)".
 */
inline std::string to_string(const Newton3Option &option) {
  switch (option) {
    case autopas::Newton3Option::enabled: {
      return "enabled";
    }
    case autopas::Newton3Option::disabled: {
      return "disabled";
    }
  }
  // do not implement default case to provoke compiler warnings if new options are introduced.
  return "Unknown option (" + std::to_string(option) + ")";
}

/**
 * Converts a SelectorStrategy to its respective string representation.
 * @param option
 * @return The string representation or "Unknown option (<IntValue>)".
 */
inline std::string to_string(const autopas::SelectorStrategyOption &option) {
  switch (option) {
    case autopas::SelectorStrategyOption::fastestAbs: {
      return "Fastest-Absolute-Value";
    }
    case autopas::SelectorStrategyOption::fastestMean: {
      return "Fastest-Mean-Value";
    }
    case autopas::SelectorStrategyOption::fastestMedian: {
      return "Fastest-Median-Value";
    }
  }
  // do not implement default case to provoke compiler warnings if new options are introduced.
  return "Unknown SelectorStrategyOption (" + std::to_string(option) + ")";
}

/**
 * Converts a DataLayoutOption to its respective string representation.
 * @param option
 * @return The string representation or "Unknown option (<IntValue>)".
 */
inline std::string to_string(const DataLayoutOption &option) {
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
  return "Unknown DataLayoutOption (" + std::to_string(option) + ")";
}

/**
 * Converts a ContainerOption to its respective string representation.
 * @param option
 * @return The string representation or "Unknown option (<IntValue>)".
 */
inline std::string to_string(const ContainerOption &option) {
  switch (option) {
    case autopas::ContainerOption::directSum: {
      return "DirectSum";
    }
    case autopas::ContainerOption::linkedCells: {
      return "LinkedCells";
    }
    case autopas::ContainerOption::verletLists: {
      return "VerletLists";
    }
    case autopas::ContainerOption::verletListsCells: {
      return "VerletListsCells";
    }
    case autopas::ContainerOption::verletClusterLists: {
      return "VerletClusterLists";
    }
  }
  // do not implement default case to provoke compiler warnings if new options are introduced.
  return "Unknown ContainerOption (" + std::to_string(option) + ")";
}

/**
 * Converts a TraversalOption to its respective string representation.
 * @param option
 * @return The string representation or "Unknown option (<IntValue>)".
 */
inline std::string to_string(const TraversalOption &option) {
  switch (option) {
    case autopas::TraversalOption::dummyTraversal: {
      return "dummyTraversal";
    }
    case autopas::TraversalOption::c01: {
      return "c01";
    }
    case autopas::TraversalOption::c04SoA: {
      return "c04SoA";
    }
    case autopas::TraversalOption::c08: {
      return "c08";
    }
    case autopas::TraversalOption::c18: {
      return "c18";
    }
    case autopas::TraversalOption::sliced: {
      return "sliced";
    }
    case autopas::TraversalOption::directSumTraversal: {
      return "directSum";
    }
    case autopas::TraversalOption::c01Verlet: {
      return "verlet-c01";
    }
    case autopas::TraversalOption::c18Verlet: {
      return "verlet-c18";
    }
    case autopas::TraversalOption::slicedVerlet: {
      return "verlet-sliced";
    }
    case autopas::TraversalOption::c01Cuda: {
      return "cuda-c01";
    }
    case autopas::TraversalOption::verletTraversal: {
      return "verlet-lists";
    }
    case autopas::TraversalOption::c01CombinedSoA: {
      return "c01-combined-SoA";
    }
  }
  // do not implement default case to provoke compiler warnings if new options are introduced.
  return "Unknown TraversalOption (" + std::to_string(option) + ")";
}

/**
 * Converts a TuningStrategyOption to its respective string representation.
 * @param option
 * @return The string representation or "Unknown option (<IntValue>)".
 */
inline std::string to_string(const TuningStrategyOption &option) {
  switch (option) {
    case autopas::TuningStrategyOption::fullSearch: {
      return "full-Search";
    }
  }
  // do not implement default case to provoke compiler warnings if new options are introduced.
  return "Unknown TuningStrategyOption (" + std::to_string(option) + ")";
}
/**
 * Converts a double to its respective string representation.
 * @param value
 * @return The string representation.
 */
inline std::string to_string(const double &value) { return std::to_string(value); }
/**
 * All accepted delimiters to split input strings.
 */
constexpr char delimiters[] = " ,;|/";
/**
 * Regex for all delimiters to split input strings.
 */
constexpr char delimitersRgx[] = "[\\s,;|/]";
/**
 * Regex for all but delimiters to split input strings as regex.
 */
constexpr char delimitersRgxInv[] = "[^\\s,;|/]";

/**
 * Splits a string by multiple delimiters.
 * @param searchString
 * @param delimiters
 * @return Set of substrings.
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
 * Converts a string of options to a set of enums. The options are expected to be lower case.
 * Allowed delimiters can be found in autopas::utils::StringUtils::delimiters
 *
 * Possible options: enabled, disabled
 *
 * @param newton3OptionsString String containing newton3 options.
 * @param ignoreUnknownOptions If set to false, a 'autopas::Newton3Option(-1)' will be inserted in the return set
 * for each not parsable word.
 * @return Set of Newton3Option enums. If no valid option was found and unknown options are ignored
 * the empty set is returned.
 */
inline std::set<autopas::Newton3Option> parseNewton3Options(const std::string &newton3OptionsString,
                                                            bool ignoreUnknownOptions = true) {
  std::set<autopas::Newton3Option> newton3Options;

  auto words = tokenize(newton3OptionsString, delimiters);

  for (auto &word : words) {
    if (word.find("on") != std::string::npos or word.find("true") != std::string::npos or
        word.find("enable") != std::string::npos) {
      newton3Options.insert(Newton3Option::enabled);
    } else if (word.find("of") != std::string::npos or word.find("false") != std::string::npos or
               word.find("disable") != std::string::npos) {
      newton3Options.insert(Newton3Option::disabled);
    } else if (not ignoreUnknownOptions) {
      newton3Options.insert(autopas::Newton3Option(-1));
    }
  }

  return newton3Options;
}

/**
 * Converts a string of options to a set of enums. The options are expected to be lower case.
 * Allowed delimiters can be found in autopas::utils::StringUtils::delimiters
 *
 * Possible options: c01, c08, c18, direct, sliced, verlet01, verlet18, verlet-sliced, c01-combined-SoA
 *
 * @param traversalOptionsString String containing traversal options.
 * @param ignoreUnknownOptions If set to false, a 'autopas::TraversalOption(-1)' will be inserted in the return set
 * for each not parsable word.
 * @return Set of TraversalOption enums. If no valid option was found and unknown options are ignored the empty
 * Set is returned.
 */
inline std::set<autopas::TraversalOption> parseTraversalOptions(const std::string &traversalOptionsString,
                                                                bool ignoreUnknownOptions = true) {
  std::set<autopas::TraversalOption> traversalOptions;

  auto words = tokenize(traversalOptionsString, delimiters);

  for (auto &word : words) {
    if (word.find("verlet-lists") != std::string::npos) {
      traversalOptions.insert(autopas::TraversalOption::verletTraversal);
    } else if (word.find("01") != std::string::npos) {
      if (word.find("cuda") != std::string::npos) {
        traversalOptions.insert(autopas::TraversalOption::c01Cuda);
      } else if (word.find("com") != std::string::npos) {
        traversalOptions.insert(autopas::TraversalOption::c01CombinedSoA);
      } else if (word.find('v') != std::string::npos) {
        traversalOptions.insert(autopas::TraversalOption::c01Verlet);
      } else {
        traversalOptions.insert(autopas::TraversalOption::c01);
      }
    } else if (word.find("c08") != std::string::npos) {
      traversalOptions.insert(autopas::TraversalOption::c08);
    } else if (word.find("c04s") != std::string::npos) {
      traversalOptions.insert(autopas::TraversalOption::c04SoA);
    } else if (word.find("18") != std::string::npos) {
      if (word.find('v') != std::string::npos)
        traversalOptions.insert(autopas::TraversalOption::c18Verlet);
      else
        traversalOptions.insert(autopas::TraversalOption::c18);
    } else if (word.find("dir") != std::string::npos) {
      traversalOptions.insert(autopas::TraversalOption::directSumTraversal);
    } else if (word.find("sli") != std::string::npos) {
      if (word.find('v') != std::string::npos)
        traversalOptions.insert(autopas::TraversalOption::slicedVerlet);
      else
        traversalOptions.insert(autopas::TraversalOption::sliced);
    } else if (not ignoreUnknownOptions) {
      traversalOptions.insert(autopas::TraversalOption(-1));
    }
  }
  return traversalOptions;
}

/**
 * Converts a string of options to a set of enums. The options are expected to be lower case.
 * Allowed delimiters can be found in autopas::utils::StringUtils::delimiters
 *
 * Possible options: directSum, linkedCells, verletLists, vcells, vcluster
 *
 * @param containerOptionsString String containing container options.
 * @param ignoreUnknownOptions If set to false, a 'autopas::ContainerOption(-1)' will be inserted in the return set
 * for each not parsable word.
 * @return Set of ContainerOption enums. If no valid option was found and unknown options are ignored the empty
 * set is returned.
 */
inline std::set<autopas::ContainerOption> parseContainerOptions(const std::string &containerOptionsString,
                                                                bool ignoreUnknownOptions = true) {
  std::set<autopas::ContainerOption> containerOptions;

  auto words = tokenize(containerOptionsString, delimiters);

  for (auto &word : words) {
    if (word.find("dir") != std::string::npos or word.find("ds") != std::string::npos) {
      containerOptions.insert(autopas::ContainerOption::directSum);
    } else if (word.find("linked") != std::string::npos or word.find("lc") != std::string::npos) {
      containerOptions.insert(autopas::ContainerOption::linkedCells);
    } else if (word.find('v') != std::string::npos) {
      if (word.find("cl") != std::string::npos) {
        containerOptions.insert(autopas::ContainerOption::verletClusterLists);
      } else if (word.find("cel") != std::string::npos) {
        containerOptions.insert(autopas::ContainerOption::verletListsCells);
      } else {
        containerOptions.insert(autopas::ContainerOption::verletLists);
      }
    } else if (not ignoreUnknownOptions) {
      containerOptions.insert(autopas::ContainerOption(-1));
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
inline autopas::SelectorStrategyOption parseSelectorStrategy(const std::string &selectorStrategyString) {
  // hack to initialize the enum out of range as an error value.
  auto selectorStrategy(autopas::SelectorStrategyOption(-1));
  if (selectorStrategyString.find("abs") != std::string::npos or
      selectorStrategyString.find("min") != std::string::npos) {
    selectorStrategy = autopas::SelectorStrategyOption::fastestAbs;
  } else if (selectorStrategyString.find("mea") != std::string::npos) {
    selectorStrategy = autopas::SelectorStrategyOption::fastestMean;
  } else if (selectorStrategyString.find("med") != std::string::npos) {
    selectorStrategy = autopas::SelectorStrategyOption::fastestMedian;
  }
  return selectorStrategy;
}

/**
 * Converts a string containing a data layout to an enum. The option is expected to be lower case.
 *
 * Possible options: aos, soa
 *

 * @param dataLayoutsSting String containing the data layout option.
 * @param ignoreUnknownOptions If set to false, a 'autopas::DataLayoutOption(-1)' will be inserted in the return set
 * for each not parsable word.
 * @return An enum representing the data layout. If no valid option was found and unknown options are ignored the empty
 * set is returned.
 */
inline std::set<autopas::DataLayoutOption> parseDataLayout(const std::string &dataLayoutsSting,
                                                           bool ignoreUnknownOptions = true) {
  auto words = tokenize(dataLayoutsSting, delimiters);

  std::set<autopas::DataLayoutOption> dataLayouts;

  for (auto &word : words) {
    if (word.find("aos") != std::string::npos or word.find("array-of-struct") != std::string::npos) {
      dataLayouts.insert(autopas::DataLayoutOption::aos);
    } else if (word.find("soa") != std::string::npos or word.find("-of-array") != std::string::npos) {
      dataLayouts.insert(autopas::DataLayoutOption::soa);
    } else if (word.find("cuda") != std::string::npos) {
      dataLayouts.insert(autopas::DataLayoutOption::cuda);
    } else if (not ignoreUnknownOptions) {
      // hack to initialize the enum out of range as an error value.
      dataLayouts.insert(autopas::DataLayoutOption(-1));
    }
  }
  return dataLayouts;
}

/**
 * Converts a string containing a tuning strategy to an enum. The option is expected to be lower case.
 *
 * Possible options: full-search
 *
 * @param tuningStrategyString String containing the tuning strategy option
 * @return An enum representing the tuningStrategy. If no valid option was found an error value of -1 is returned.
 */
inline autopas::TuningStrategyOption parseTuningStrategyOption(const std::string &tuningStrategyString) {
  // hack to initialize the enum out of range as an error value.
  auto tuningStrategy(autopas::TuningStrategyOption(-1));
  if (tuningStrategyString.find("full") != std::string::npos or tuningStrategyString.find("ex") != std::string::npos) {
    tuningStrategy = autopas::TuningStrategyOption::fullSearch;
  }
  return tuningStrategy;
}

/**
 * Converts a string to a set of doubles.
 * @param doubleString String containing doubles.
 * @param ignoreUnknownOptions If set to false, 'nan' will be inserted in the return set
 * for each not parsable word.
 * @return Set of doubles. If no valid double was found and unknown options are ignored the empty
 * set is returned.
 */
inline std::set<double> parseDoubles(const std::string &doubleString, bool ignoreUnknownOptions = true) {
  auto words = tokenize(doubleString, delimiters);

  std::set<double> doubles;

  for (auto &word : words) {
    try {
      double value = stod(word);
      doubles.insert(value);
    } catch (const std::exception &) {
      if (not ignoreUnknownOptions) {
        doubles.insert(std::nan(""));
      }
    }
  }
  return doubles;
}

/**
 * Converts a string to a NumberSet<double>.
 * @param setString String containing the set.
 * @param ignoreUnknownOptions If set to false, 'nan' will be inserted in the return set
 * for each not parsable word.
 * @return NumberSet<double>. If no valid double was found and unknown options are ignored the empty
 * set is returned.
 */
inline std::unique_ptr<autopas::NumberSet<double>> parseNumberSet(const std::string &setString,
                                                                  bool ignoreUnknownOptions = true) {
  // try to match an interval [x,y]
  std::regex rgx(
      "\\["         // open square bracket
      "([^,]++)"    // any number of non-comma chars (1st Capturing Group)
      ","           // comma
      "([^\\]]++)"  // any number of non-closing-bracket chars (2nd Capturing Group)
      "\\]"         // closing square bracket
  );
  std::smatch matches;
  if (std::regex_match(setString, matches, rgx)) {
    try {
      double min = stod(matches.str(1));
      double max = stod(matches.str(2));
      return std::make_unique<autopas::NumberInterval<double>>(min, max);
    } catch (const std::exception &) {
      // try parseDoubles instead
    }
  }

  std::set<double> values = autopas::utils::StringUtils::parseDoubles(setString, ignoreUnknownOptions);
  return std::make_unique<autopas::NumberSetFinite<double>>(values);
}
}  // namespace StringUtils
}  // namespace utils
}  // namespace autopas
