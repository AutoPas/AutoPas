/**
 * @file ConfigurationAndRankIteratorHandler.cpp
 * @author W. Thieme
 * @date 21.06.2020
 */

#include "ConfigurationAndRankIteratorHandler.h"

#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/containers/CompatibleTraversals.h"

namespace autopas::utils {

inline void ConfigurationAndRankIteratorHandler::advanceConfigIterators() {
  // advance to the next valid config
  ++_newton3It;
  if (_newton3It != _newton3Options.end()) return;
  _newton3It = _newton3Options.begin();
  ++_dataLayoutIt;
  if (_dataLayoutIt != _dataLayoutOptions.end()) return;
  _dataLayoutIt = _dataLayoutOptions.begin();
  ++_loadEstimatorIt;
  if (_loadEstimatorIt != _allowedAndApplicableLoadEstimatorOptions.end()) return;
  _loadEstimatorIt = _allowedAndApplicableLoadEstimatorOptions.begin();
  ++_traversalIt;
  if (_traversalIt != _allowedAndApplicableTraversalOptions.end()) {
    selectLoadEstimatorsForCurrentContainerAndTraversal();
    return;
  } else {
    _traversalIt = _allowedAndApplicableTraversalOptions.begin();
    selectLoadEstimatorsForCurrentContainerAndTraversal();
  }
  ++_cellSizeFactorIt;
  if (_cellSizeFactorIt != _cellSizeFactors.end()) return;
  _cellSizeFactorIt = _cellSizeFactors.begin();
  ++_containerIt;
  if (_containerIt != _containers.end()) {
    selectTraversalsForCurrentContainer();
    selectLoadEstimatorsForCurrentContainerAndTraversal();
  }
}

void ConfigurationAndRankIteratorHandler::advanceIterators(const int numConfigs, const int commSize) {
  if (numConfigs >= commSize or _remainingBlockSize == 0) {
    advanceConfigIterators();
  }

  if (commSize >= numConfigs or _remainingBlockSize == 0) {
    // advance to the next rank
    ++_rankIterator;

    // advance offset to the position relative to the first rank with the same configuration.
    ++_infiniteCellSizeFactorsOffset;
  }

  // Set information necessary to compute the next block.
  // A block is a set of consecutive pairs of configurations and ranks where only the one with more possible values
  // changes. e.g. if there are three configurations for every rank, a block's size is three.
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

void ConfigurationAndRankIteratorHandler::reset(const int numConfigs, const int commSize) {
  _containerIt = _containers.begin();
  _cellSizeFactorIt = _cellSizeFactors.begin();
  _dataLayoutIt = _dataLayoutOptions.begin();
  _newton3It = _newton3Options.begin();
  selectTraversalsForCurrentContainer();
  selectLoadEstimatorsForCurrentContainerAndTraversal();

  _rankIterator = 0;
  _remainingBlockSize = commSize >= numConfigs ? commSize / numConfigs - 1 : numConfigs / commSize - 1;
  _remainder = commSize >= numConfigs ? commSize % numConfigs : numConfigs % commSize;
  _infiniteCellSizeFactorsOffset = 0;
  _infiniteCellSizeFactorsBlockSize = commSize >= numConfigs ? commSize / numConfigs : numConfigs / commSize;
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

void ConfigurationAndRankIteratorHandler::selectLoadEstimatorsForCurrentContainerAndTraversal() {
  // if load estimators are not applicable LoadEstimatorOption::none is returned.
  _allowedAndApplicableLoadEstimatorOptions =
      loadEstimators::getApplicableLoadEstimators(*_containerIt, *_traversalIt, _allowedLoadEstimatorOptions);
  _loadEstimatorIt = _allowedAndApplicableLoadEstimatorOptions.begin();
}

}  // namespace autopas::utils