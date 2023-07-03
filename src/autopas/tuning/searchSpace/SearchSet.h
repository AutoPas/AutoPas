/**
 * @file SearchSet.h
 * @author F. Gratl
 * @date 26.06.23
 */

#pragma once

#include <c++/11/optional>
#include <cwchar>
#include <memory>
#include <optional>
#include <set>
#include <sstream>

#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/NumberSet.h"

namespace autopas {

/**
 * Forward declaration.
 */
class SearchSetIterator;

/**
 * Class representing part of a search space that shares a set of discrete and continuous options.
 * The set that is represented is the cross product of all configurations (=discrete options) and
 * all cellSizeFactors (=continuous options).
 */
class SearchSet {
 public:
  /**
   * Set of discrete options. CSF in these configurations is always -1.
   */
  std::set<Configuration> configurations;
  /**
   * cellSize options.
   */
  std::unique_ptr<NumberSet<double>> cellSizeFactors;

  /**
   * Default constructor, creating an empty set.
   */
  SearchSet() = default;

  /**
   * Constructor.
   * @param configurations Set of
   * @param cellSizeFactors
   */
  SearchSet(const std::set<Configuration> &configurations, std::unique_ptr<NumberSet<double>> cellSizeFactors)
      : configurations(configurations), cellSizeFactors(std::move(cellSizeFactors)) {}

  /**
   * Copy constructor.
   * @param other
   */
  SearchSet(const SearchSet &other)
      : configurations(other.configurations), cellSizeFactors(other.cellSizeFactors->clone()) {}

  /**
   * Move constructor
   * @param other
   */
  SearchSet(SearchSet &&other) noexcept
      : configurations(std::move(other.configurations)), cellSizeFactors(std::move(other.cellSizeFactors)) {}

  /**
   * Swap the members of this and other.
   * @param other
   */
  void swap(SearchSet &other);

  /**
   * Copy assignment operator.
   * @param rhs
   * @return
   */
  SearchSet &operator=(const SearchSet &rhs);

  /**
   * Move assignment operator.
   * @param rhs
   * @return
   */
  SearchSet &operator=(SearchSet &&rhs) noexcept;

  /**
   * Create an iterator at the start of the search set.
   * @param csfStepSize
   * @return
   */
  SearchSetIterator begin(double csfStepSize = .1);

  /**
   * Create an iterator at the end of the search set.
   * @return
   */
  SearchSetIterator end();

  /**
   * Indicates if the set is empty.
   * @return
   */
  bool empty() const;

  /**
   * Removes a config from the set and returns a vector of sets that represent the remaining configuration space.
   *
   * Example:
   *
   * [[c0, c1, c2], [n0, n1, n2]] \ (c1, n1) = {[[c0, c2], [n0, n1, n2]],
   *                                            [[c1], [n0, n2]]}
   * [[c0]        , [n0, n1, n2]] \ (c0, n1) = {[[c0], [n0, n2]]}
   * [[c0, c1, c2], [n0 - nX]]    \ (c1, n1) = {[[c0, c2], [n0 - nX]],
   *                                            [[c1], [n0, n1-1]],
   *                                            [[c1], [n1+1, nX]]}
   * [[c0, c1, c2], [n0 - nX]]    \ (c1, nX) = {[[c0, c2], [n0 - nX]],
   *                                            [[c1], [n0, nX-1]]}
   * [[c0, c1, c2], [n0 - n0]]    \ (c1, n0) = {[[c0, c2], [n0 - nX]]}
   * [[c0]        , [n0]]         \ (c0, n0) = {}
   *
   * @param configuration
   * @return
   */
  [[nodiscard]] std::vector<SearchSet> deleteConfig(const Configuration &configuration) const;

  bool operator==(const SearchSet &other) const {
    return configurations == other.configurations and *cellSizeFactors == *cellSizeFactors;
  }

  /**
   * Creates a string representation of the set.
   * @return
   */
  std::string toString() const;
};
}  // namespace autopas
