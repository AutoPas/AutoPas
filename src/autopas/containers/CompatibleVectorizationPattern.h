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
#include "autopas/options/VectorizationPatternOption.h"

namespace autopas::compatibleVectorizationPattern {

/**
 * Returns a set of vectorization patterns compatible with the container.
 *
 * @param container
 * @return compatible load vectorization pattenrs
 */
static std::set<autopas::VectorizationPatternOption> allCompatibleVectorizationPattern(
    autopas::ContainerOption container) {
  switch (container) {
    case ContainerOption::verletLists:
    case ContainerOption::verletListsCells:
    case ContainerOption::pairwiseVerletLists:
    case ContainerOption::varVerletListsAsBuild: {
      return std::set<autopas::VectorizationPatternOption>{VectorizationPatternOption::p1xVec};
    }
    default: {
      return std::set<autopas::VectorizationPatternOption>{autopas::VectorizationPatternOption::getAllOptions()};
    }
  }
}

}  // namespace autopas::compatibleVectorizationPattern
