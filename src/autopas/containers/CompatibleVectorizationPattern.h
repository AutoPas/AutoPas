/**
 * @file CompatibleVectorizationPattern.h
 *
 * @date 24.12.25
 * @author D. Martin
 */

#pragma once

#include <array>
#include <set>
#include <vector>

#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/VectorizationPatternOption.h"

namespace autopas::compatibleVectorizationPattern {

/**
 * Returns a set of vectorization patterns compatible with the container.
 *
 * @param container
 * @param dataLayout
 * @return compatible load vectorization pattenrs
 */
static std::set<VectorizationPatternOption> allCompatibleVectorizationPattern(const ContainerOption container,
                                                                              const DataLayoutOption dataLayout) {
  switch (dataLayout) {
    case DataLayoutOption::aos:
      // For AoS, Vectorization Patterns are not applicable.
      return std::set<VectorizationPatternOption>{VectorizationPatternOption::NA};
    case DataLayoutOption::soa:
      switch (container) {
      case ContainerOption::verletLists:
      case ContainerOption::verletListsCells:
      case ContainerOption::pairwiseVerletLists:
      case ContainerOption::varVerletListsAsBuild:
      case ContainerOption::directSum: {
        return std::set<VectorizationPatternOption>{VectorizationPatternOption::p1xVec};
      }
      default: {
        return std::set<VectorizationPatternOption>{VectorizationPatternOption::getAllApplicablePatterns()};
      }
      }
    default:
      utils::ExceptionHandler::exception("Unknown data layout {}.", dataLayout.to_string());
      return {};
  }
}

}  // namespace autopas::compatibleVectorizationPattern
