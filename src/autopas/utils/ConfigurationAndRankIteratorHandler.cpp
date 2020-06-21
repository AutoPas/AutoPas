/**
 * @file ConfigurationAndRankIteratorHandler.cpp
 * @author W. Thieme
 * @date 21.06.2020
 */

#include "ConfigurationAndRankIteratorHandler.h"

namespace autopas::utils {

inline void ConfigurationAndRankIteratorHandler::advanceConfigIterators() {
  // advance to the next valid config
  ++_newton3It;
  if (_newton3It != _newton3Options.end()) return;
  _newton3It = _newton3Options.begin();
  ++_dataLayoutIt;
  if (_dataLayoutIt != _dataLayoutOptions.end()) return;
  _dataLayoutIt = _dataLayoutOptions.begin();
  ++_traversalIt;
  if (_traversalIt != _allowedAndApplicableTraversalOptions.end()) return;
  _traversalIt = _allowedAndApplicableTraversalOptions.begin();
  ++_cellSizeFactorIt;
  if (_cellSizeFactorIt != _cellSizeFactors.end()) return;
  _cellSizeFactorIt = _cellSizeFactors.begin();
  ++_containerIt;
  selectTraversalsForCurrentContainer();
}

void ConfigurationAndRankIteratorHandler::advanceIterators(const int numConfigs, const int commSize) {
  if (numConfigs >= commSize or _remainingBlockSize == 0) {
    advanceConfigIterators();
  }

  if (commSize >= numConfigs or _remainingBlockSize == 0) {
    // advance to the next rank
    ++_rankIterator;

    // advance offset to the position relative to the first rank with the same configuration
    ++_infiniteCellSizeFactorsOffset;
  }

  // Set information necessary to compute the next block.
  // Block here means either a block of ranks that all have the same configuration or a set of configuration that all
  // have the same ranks.
  if (_remainingBlockSize == 0) {
    if (numConfigs >= commSize) {
      _remainingBlockSize = numConfigs / commSize;
    } else {
      _remainingBlockSize = commSize / numConfigs;

      _infiniteCellSizeFactorsBlockSize = _remainingBlockSize;
      _infiniteCellSizeFactorsOffset = 0;
    }
    if (_remainder > 0) {
      ++_remainingBlockSize;
      --_remainder;
    }
  }

  --_remainingBlockSize;
}

void ConfigurationAndRankIteratorHandler::selectTraversalsForCurrentContainer() {
  // get all traversals of the container and restrict them to the allowed ones
  const std::set<TraversalOption> &allContainerTraversals =
      compatibleTraversals::allCompatibleTraversals(*_containerIt);
  _allowedAndApplicableTraversalOptions.clear();
  std::set_intersection(
      _allowedTraversalOptions.begin(), _allowedTraversalOptions.end(), allContainerTraversals.begin(),
      allContainerTraversals.end(),
      std::inserter(_allowedAndApplicableTraversalOptions, _allowedAndApplicableTraversalOptions.begin()));
  _traversalIt = _allowedAndApplicableTraversalOptions.begin();
}

} // namespace autopas::utils