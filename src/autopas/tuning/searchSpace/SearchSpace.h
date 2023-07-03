/**
 * @file SearchSpace.h
 * @author F. Gratl
 * @date 27.06.23
 */

#pragma once

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <vector>

#include "SearchSet.h"
#include "utils/NumberInterval.h"
#include "utils/NumberSetFinite.h"

namespace autopas {

/**
 * Manager class for everything contained in the search space.
 */
class SearchSpace {
 public:
  /**
   * Constructor
   * @param searchSets
   * @param csfsAreContinuous
   */
  explicit SearchSpace(const std::vector<SearchSet> &searchSets);

  /**
   * Merge as many sets as possible to keep the vector size down.
   */
  void simplifyInternalSets();

 private:
  /**
   * The search space is a vector of sets.
   *
   * This allows to create a space of arbitrary configurations (not just a cross-product of options)
   * and associating discrete options with continuous ones.
   */
  std::vector<SearchSet> _searchSets;
  /**
   * Flag to indicate if all cell size factors are continuous intervals or finite sets.
   */
  bool _csfsAreContinuous;
};
}  // namespace autopas
