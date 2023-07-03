/**
 * @file SearchSetIterator.h
 * @author F. Gratl
 * @date 26.06.23
 */

#pragma once

#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/tuning/Configuration.h"

namespace autopas {

/**
 * Forward declaration.
 */
class SearchSet;

/**
 * Iterator for one search set.
 *
 * Iteration order follows Configuration::operator<().
 *
 */
class SearchSetIterator {
 public:
  /**
   * Iterator for the configuration option.
   */
  std::set<Configuration>::iterator configurationsIter;
  /**
   * Current cell size factor.
   */
  double cellSizeFactor;

  /**
   * Constructor
   * @param searchSet
   * @param csfStepSize Step size for sampling a continuous cell size factor.
   */
  explicit SearchSetIterator(SearchSet &searchSet, double csfStepSize = .1);

  /**
   * Constructor with a given position.
   * @param containerOptionIter
   * @param traversalOptionIter
   * @param newton3OptionIter
   * @param dataLayoutOptionIter
   * @param loadEstimatorOptionIter
   * @param cellSizeFactor
   * @param searchSet
   */
  SearchSetIterator(const std::set<Configuration>::iterator &configurationsIter, double cellSizeFactor,
                    SearchSet &searchSet);

  /**
   * Increment operator.
   * @return
   */
  SearchSetIterator &operator++();

  /**
   * Sets the Cell Size Factor to a new value and resets all faster iterating options.
   * @note Useful for manually increasing the iterator in the case that CSF is continuous.
   * @param newCsf
   */
  void setCellSizeFactor(double newCsf);

 private:
  /**
   * Reference to the set that is iterated.
   */
  SearchSet &_searchSet;
  /**
   * Step size for the cell size factor if it is continuous.
   * If it is discrete, the set elements are iterated.
   */
  double _csfStepSize;
};  // namespace autopas
}  // namespace autopas