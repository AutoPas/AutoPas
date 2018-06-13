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
#include "utils/ExceptionHandler.h"

namespace autopas {

enum TraversalOptions {
  c08 = 0,
  sliced = 1,
};

/**
 * Provides a way to iterate over the possible choices of TraversalOption.
 */
static std::vector<TraversalOptions> allTraversalOptions = {TraversalOptions::c08,
                                                            TraversalOptions::sliced};

template<class ParticleCell>
class TraversalSelector {
 public:
  TraversalSelector(const std::array<unsigned long, 3> &dims,
                    unsigned int retuneInterval,
                    const std::vector<TraversalOptions> &allowedTraversalOptions
  ) : _dims(dims),
      _retuneInterval(retuneInterval),
      _retuneCounter(0),
      _allowedTraversalOptions(allowedTraversalOptions) {
  }

  template<class CellFunctor>
  CellPairTraversal<ParticleCell, CellFunctor> *getOptimalTraversal(CellFunctor &cellFunctor);

 private:

  template<class CellFunctor>
  std::vector<CellPairTraversal<ParticleCell, CellFunctor> *> *generateTraversals(CellFunctor &cellFunctor);
  template<class CellFunctor>
  CellPairTraversal<ParticleCell, CellFunctor> *chooseOptimalTraversal(std::vector<CellPairTraversal<ParticleCell,
                                                                                                     CellFunctor> *> &traversals);

  // for each encountered cell processor save the optimal traversal. The cell processor is saved through its hash
  std::unordered_map<size_t, TraversalOptions> _optimalTraversalOptions;
  const std::array<unsigned long, 3> &_dims;
  unsigned int _retuneInterval, _retuneCounter;
  const std::vector<TraversalOptions> &_allowedTraversalOptions;
};

template<class ParticleCell>
template<class CellFunctor>
std::vector<CellPairTraversal<ParticleCell, CellFunctor> *> *TraversalSelector<ParticleCell>::generateTraversals(
    CellFunctor &cellFunctor) {

  auto *traversals = new std::vector<CellPairTraversal<ParticleCell, CellFunctor> *>;

  for (auto &option : _allowedTraversalOptions) {
    switch (option) {
      case c08: {
        traversals->push_back(new C08Traversal<ParticleCell, CellFunctor>(_dims, &cellFunctor));
        break;
      }
      case sliced : {
        traversals->push_back(new SlicedTraversal<ParticleCell, CellFunctor>(_dims, &cellFunctor));
        break;
      }
      default: {
        AutoPasLogger->warn("Traversal type {} is not a known type!", option);
      }
    }
  }

  assert(traversals->size() > 0);
  return traversals;
}

template<class ParticleCell>
template<class CellFunctor>
CellPairTraversal<ParticleCell, CellFunctor> *TraversalSelector<ParticleCell>::chooseOptimalTraversal(std::vector<
    CellPairTraversal<ParticleCell, CellFunctor> *> &traversals) {
  //TODO: Autotuning goes here
  auto optimalTraversal = traversals.front();

  if (dynamic_cast<C08Traversal<ParticleCell, CellFunctor> *>(optimalTraversal))
    _optimalTraversalOptions[typeid(CellFunctor).hash_code()] = TraversalOptions::c08;
  else if (dynamic_cast<SlicedTraversal<ParticleCell, CellFunctor> *>(optimalTraversal))
    _optimalTraversalOptions[typeid(CellFunctor).hash_code()] = TraversalOptions::sliced;
  else {
    utils::ExceptionHandler::exception("Optimal traversal has unknown type!");
  }

  return optimalTraversal;
}

template<class ParticleCell>
template<class CellFunctor>
CellPairTraversal<ParticleCell, CellFunctor> *TraversalSelector<ParticleCell>::getOptimalTraversal(CellFunctor &cellFunctor) {
  CellPairTraversal<ParticleCell, CellFunctor> *traversal;

  if (_retuneCounter == 0 || _optimalTraversalOptions.find(typeid(CellFunctor).hash_code()) == _optimalTraversalOptions.end()) {
    _retuneCounter = _retuneInterval;
    traversal = chooseOptimalTraversal<CellFunctor>(*(generateTraversals<CellFunctor>(cellFunctor)));
  } else {
    switch (_optimalTraversalOptions[typeid(CellFunctor).hash_code()]) {
      case c08 : {
        traversal = new C08Traversal<ParticleCell, CellFunctor>(_dims, &cellFunctor);
        break;
      }
      case sliced : {
        traversal = new SlicedTraversal<ParticleCell, CellFunctor>(_dims, &cellFunctor);
        break;
      }
      default: {
        utils::ExceptionHandler::exception("Optimal traversal option is unknown!");
      }
    }
  }
  --_retuneCounter;
  return traversal;
}
}