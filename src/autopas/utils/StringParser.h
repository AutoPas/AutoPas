/**
 * @file StringParser.h
 * @author F. Gratl
 * @date 1/14/19
 */

#pragma once

#include <string>
#include <vector>
#include "autopas/containers/cellPairTraversals/TraversalInterface.h"

namespace autopas {
namespace utils {
/**
 * Some functions to parse enums from (input-) strings.
 */
namespace StringParser {

/**
 * Converts a string of options to a vector of enums. The options are expected to be lower case.
 * @param traversalOptionsString String containing traversal options.
 * @return Vector of TraversalOption enums.
 */
static std::vector<autopas::TraversalOptions> parseTraversalOptions(std::string &traversalOptionsString) {
  std::vector<autopas::TraversalOptions> traversalOptions;
  if (traversalOptionsString.find("c08") != std::string::npos) {
    traversalOptions.push_back(autopas::TraversalOptions::c08);
  }
  if (traversalOptionsString.find("c01") != std::string::npos) {
    traversalOptions.push_back(autopas::TraversalOptions::c01);
  }
  if (traversalOptionsString.find("c18") != std::string::npos) {
    traversalOptions.push_back(autopas::TraversalOptions::c18);
  }
  if (traversalOptionsString.find("sli") != std::string::npos) {
    traversalOptions.push_back(autopas::TraversalOptions::sliced);
  }
  if (traversalOptionsString.find("dir") != std::string::npos) {
    traversalOptions.push_back(autopas::TraversalOptions::directSumTraversal);
  }
  return traversalOptions;
}
}  // namespace StringParser
}  // namespace utils
}  // namespace autopas