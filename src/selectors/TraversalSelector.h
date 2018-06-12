/**
 * TraversalSelector.h
 *
 *  Created on: 6/11/18
 *     Aauthor: F. Gratl
 */

#pragma once

#include <array>
#include <vector>
#include <containers/cellPairTraversals/CellPairTraversal.h>
#include <containers/cellPairTraversals/C08Traversal.h>
#include <containers/cellPairTraversals/SlicedTraversal.h>
#include "utils/Logger.h"
namespace autopas {

enum TraversalOptions {
  c08 = 0,
  sliced = 1,
};

/**
 * Provides a way to iterate over the possible choices of TraversalOption.
 */
static std::array<TraversalOptions, 3> allTraversalOptions = {TraversalOptions::c08,
                                                              TraversalOptions::sliced};

template<class ParticleCell>
class TraversalSelector {
 public:
  TraversalSelector(const std::array<unsigned long, 3> &dims,
                    unsigned int retuneInterval,
                    const std::vector<TraversalOptions> &allowedTraversalOptions
  ) : _dims(dims),
      _retuneInterval(retuneInterval),
      _allowedTraversalOptions(allowedTraversalOptions) {
  }

  template<class CellFunctor>
  CellPairTraversal<ParticleCell, CellFunctor> *getOptimalTraversal(CellFunctor &cellFunctor);

 private:
//  template<class CellFunctor>
//  using TraversalType = CellPairTraversal<ParticleCell, CellFunctor>;
//  using TraversalList = std::vector<std::unique_ptr<TraversalSelector::TraversalType>>;

  template<class CellFunctor>
  std::vector<CellPairTraversal<ParticleCell, CellFunctor>> generateTraversals(CellFunctor &cellFunctor);
  template<class CellFunctor>
  CellPairTraversal<ParticleCell, CellFunctor> *chooseOptimalTraversal(std::vector<CellPairTraversal<ParticleCell, CellFunctor>> traversals);

  unsigned int _retuneInterval, _retuneCounter;
//  std::unique_ptr<CellPairTraversal<ParticleCell, CellFunctor>> _optimalTraversal;
  const std::array<unsigned long, 3> &_dims;
  const std::vector<TraversalOptions> &_allowedTraversalOptions;
};

template<class ParticleCell>
template<class CellFunctor>
std::vector<CellPairTraversal<ParticleCell, CellFunctor>> TraversalSelector<ParticleCell>::generateTraversals(CellFunctor &cellFunctor) {

  std::vector<CellPairTraversal<ParticleCell, CellFunctor>*> traversals;

  for (auto &option : _allowedTraversalOptions) {
    switch (option) {
      case c08: {
        traversals.push_back(new C08Traversal<ParticleCell, CellFunctor>(_dims, cellFunctor));
        break;
      }
      case sliced : {
        traversals.push_back(new SlicedTraversal<ParticleCell, CellFunctor>(_dims, cellFunctor));
        break;
      }
      default: {
        AutoPasLogger->warn("Traversal type {} is not a known type!", option);
      }
    }
  }

  assert(traversals.size() > 0);
  return traversals;
}

template<class ParticleCell>
template<class CellFunctor>
CellPairTraversal<ParticleCell, CellFunctor> *TraversalSelector<ParticleCell>::chooseOptimalTraversal(std::vector<CellPairTraversal<ParticleCell, CellFunctor>> traversals) {
  //TODO: Autotuning goes here
//  _optimalTraversal = traversals.front();
  return traversals.front();
}


template<class ParticleCell>
template<class CellFunctor>
CellPairTraversal<ParticleCell, CellFunctor> *TraversalSelector<ParticleCell>::getOptimalTraversal(CellFunctor &cellFunctor) {
//  if (_retuneCounter == 0) {
    return chooseOptimalTraversal<CellFunctor>(generateTraversals<CellFunctor>(cellFunctor));
//    _retuneCounter = _retuneInterval;
//  }
//
//  --_retuneCounter;

//  return _optimalTraversal;
}
}