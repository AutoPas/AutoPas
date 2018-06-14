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
                    const std::vector<TraversalOptions> &allowedTraversalOptions
  ) : _dims(dims),
      _allowedTraversalOptions(allowedTraversalOptions) {
  }

  template<class CellFunctor>
  std::unique_ptr<CellPairTraversal<ParticleCell, CellFunctor>> getOptimalTraversal(CellFunctor &cellFunctor);
  template<class CellFunctor>
  void tune(CellFunctor &cellFunctor);

 private:

  template<class CellFunctor>
  std::vector<std::unique_ptr<CellPairTraversal<ParticleCell, CellFunctor>>> generateTraversals(CellFunctor &cellFunctor);
  template<class CellFunctor>
  std::unique_ptr<CellPairTraversal<ParticleCell, CellFunctor>>chooseOptimalTraversal(std::vector<std::unique_ptr<CellPairTraversal<ParticleCell,
                                                                                                     CellFunctor>>> &traversals);

  // for each encountered cell processor save the optimal traversal. The cell processor is saved through its hash
  std::unordered_map<size_t, TraversalOptions> _optimalTraversalOptions;
  const std::array<unsigned long, 3> &_dims;
  const std::vector<TraversalOptions> &_allowedTraversalOptions;
};

template<class ParticleCell>
template<class CellFunctor>
std::vector<std::unique_ptr<CellPairTraversal<ParticleCell, CellFunctor>>> TraversalSelector<ParticleCell>::generateTraversals(
    CellFunctor &cellFunctor) {

  std::vector<std::unique_ptr<CellPairTraversal<ParticleCell, CellFunctor>>>traversals;

  for (auto &option : _allowedTraversalOptions) {
    switch (option) {
      case c08: {
        traversals.push_back(std::make_unique<C08Traversal<ParticleCell, CellFunctor>>(_dims, &cellFunctor));
        break;
      }
      case sliced : {
        traversals.push_back(std::make_unique<SlicedTraversal<ParticleCell, CellFunctor>>(_dims, &cellFunctor));
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
std::unique_ptr<CellPairTraversal<ParticleCell, CellFunctor>> TraversalSelector<ParticleCell>::chooseOptimalTraversal(std::vector<std::unique_ptr<
    CellPairTraversal<ParticleCell, CellFunctor>>> &traversals) {
  //TODO: Autotuning goes here
  auto optimalTraversal = std::move(traversals.front());

  if (dynamic_cast<C08Traversal<ParticleCell, CellFunctor> *>(optimalTraversal.get()))
    _optimalTraversalOptions[typeid(CellFunctor).hash_code()] = TraversalOptions::c08;
  else if (dynamic_cast<SlicedTraversal<ParticleCell, CellFunctor> *>(optimalTraversal.get()))
    _optimalTraversalOptions[typeid(CellFunctor).hash_code()] = TraversalOptions::sliced;
  else {
    utils::ExceptionHandler::exception("Optimal traversal has unknown type!");
  }

  return optimalTraversal;
}

template<class ParticleCell>
template<class CellFunctor>
std::unique_ptr<CellPairTraversal<ParticleCell, CellFunctor>>TraversalSelector<ParticleCell>::getOptimalTraversal(CellFunctor &cellFunctor) {
  std::unique_ptr<CellPairTraversal<ParticleCell, CellFunctor>> traversal;

  if (_optimalTraversalOptions.find(typeid(CellFunctor).hash_code()) == _optimalTraversalOptions.end()) {
    auto generatedTraversals = generateTraversals<CellFunctor>(cellFunctor);
    traversal = chooseOptimalTraversal<CellFunctor>(generatedTraversals);
  } else {
    switch (_optimalTraversalOptions[typeid(CellFunctor).hash_code()]) {
      case c08 : {
        traversal = std::make_unique<C08Traversal<ParticleCell, CellFunctor>>(_dims, &cellFunctor);
        break;
      }
      case sliced : {
        traversal = std::make_unique<SlicedTraversal<ParticleCell, CellFunctor>>(_dims, &cellFunctor);
        break;
      }
      default: {
        utils::ExceptionHandler::exception("Optimal traversal option is unknown!");
      }
    }
  }
  return traversal;
}

template<class ParticleCell>
template<class CellFunctor>
void TraversalSelector<ParticleCell>::tune(CellFunctor &cellFunctor) {
  auto generatedTraversals = generateTraversals<CellFunctor>(cellFunctor);
  chooseOptimalTraversal<CellFunctor>(generatedTraversals);
}
}