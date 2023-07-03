/**
 * @file SearchSetIterator.cpp
 * @author F. Gratl
 * @date 26.06.23
 */

#include "SearchSetIterator.h"

#include "SearchSet.h"

namespace autopas {

SearchSetIterator &SearchSetIterator::operator++() {
  /**
   * Helper function that increases the iterator of one set or resets it to the start if the end was hit.
   * @return True if the iterator was reset / wrapped around.
   */
  auto incAndWrap = [&](auto &iter, const auto &set) {
    ++iter;
    if (iter == set.end()) {
      iter = set.begin();
      return true;
    } else {
      return false;
    }
  };

  if (not incAndWrap(configurationsIter, _searchSet.configurations)) {
    return *this;
  }
  /**
   * Bump the CellSizeFactor and check if it has breached the end of its interval.
   * @return True if now after the update CSF is out of bounds.
   */
  const auto noNewCSFLeft = [&]() {
    if (_searchSet.cellSizeFactors->isFinite()) {
      // Discrete case
      const auto nextCSFIter = std::next(_searchSet.cellSizeFactors->getAll().find(cellSizeFactor));
      if (nextCSFIter == _searchSet.cellSizeFactors->getAll().end()) {
        return true;
      } else {
        cellSizeFactor = *nextCSFIter;
      }
    } else {
      // Continuous case
      cellSizeFactor += _csfStepSize;
      if (cellSizeFactor > _searchSet.cellSizeFactors->getMax()) {
        return true;
      }
    }
    return false;
  };
  if (noNewCSFLeft()) {
    // in this case we have reached the end of the iteration so invalidate every member.
    cellSizeFactor = -1.;
    configurationsIter = _searchSet.configurations.end();
  }

  return *this;
}

void SearchSetIterator::setCellSizeFactor(double newCsf) {
  cellSizeFactor = newCsf;
  configurationsIter = _searchSet.configurations.begin();
}

SearchSetIterator::SearchSetIterator(SearchSet &searchSet, double csfStepSize)
    : configurationsIter(searchSet.configurations.begin()),
      cellSizeFactor(searchSet.cellSizeFactors->getMin()),
      _searchSet(searchSet),
      _csfStepSize(csfStepSize) {}
}  // namespace autopas