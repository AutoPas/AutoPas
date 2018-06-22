/**
 * TraversalSelector.h
 *
 *  Created on: 6/11/18
 *     Aauthor: F. Gratl
 */

#pragma once

#include <containers/cellPairTraversals/C08Traversal.h>
#include <containers/cellPairTraversals/CellPairTraversal.h>
#include <containers/cellPairTraversals/SlicedTraversal.h>
#include <array>
#include <vector>
#include "utils/ExceptionHandler.h"
#include "utils/Logger.h"

namespace autopas {

/**
 * Possible choices for the cell pair traversal.
 */
enum TraversalOptions {
  c08 = 0,
  sliced = 1,
};

/**
 * Provides a way to iterate over the possible choices of TraversalOption.
 */
static std::vector<TraversalOptions> allTraversalOptions = {TraversalOptions::c08, TraversalOptions::sliced};

/**
 * Selector for a container traversal.
 * @tparam ParticleCell
 */
template <class ParticleCell>
class TraversalSelector {
 public:
  /**
   * Constructor of the TraversalSelector class.
   * @param dims Array with the dimension lengths of the domain.
   * @param allowedTraversalOptions Vector of traversals the selector can choose from.
   */
  TraversalSelector(const std::array<unsigned long, 3> &dims,
                    const std::vector<TraversalOptions> &allowedTraversalOptions)
      : _dims(dims), _allowedTraversalOptions(allowedTraversalOptions) {}

  /**
   * Gets the optimal traversal for a given cell functor. If no traversal is selected yet a optimum search is started.
   * @tparam CellFunctor
   * @param cellFunctor
   * @return Smartpointer to the optimal traversal.
   */
  template <class CellFunctor>
  std::unique_ptr<CellPairTraversal<ParticleCell, CellFunctor>> getOptimalTraversal(CellFunctor &cellFunctor);

  /**
   * Evaluates to optimal traversal based on a given cell functor.
   * @tparam CellFunctor
   * @param cellFunctor
   */
  template <class CellFunctor>
  void tune(CellFunctor &cellFunctor);

 private:
  template <class CellFunctor>
  std::vector<std::unique_ptr<CellPairTraversal<ParticleCell, CellFunctor>>> generateTraversals(
      CellFunctor &cellFunctor);
  template <class CellFunctor>
  std::unique_ptr<CellPairTraversal<ParticleCell, CellFunctor>> chooseOptimalTraversal(
      std::vector<std::unique_ptr<CellPairTraversal<ParticleCell, CellFunctor>>> &traversals);

  // for each encountered cell processor save the optimal traversal. The cell processor is saved through its hash
  std::unordered_map<size_t, TraversalOptions> _optimalTraversalOptions;
  const std::array<unsigned long, 3> _dims;
  const std::vector<TraversalOptions> _allowedTraversalOptions;
};

template <class ParticleCell>
template <class CellFunctor>
std::vector<std::unique_ptr<CellPairTraversal<ParticleCell, CellFunctor>>>
TraversalSelector<ParticleCell>::generateTraversals(CellFunctor &cellFunctor) {
  std::vector<std::unique_ptr<CellPairTraversal<ParticleCell, CellFunctor>>> traversals;

  for (auto &option : _allowedTraversalOptions) {
    switch (option) {
      case c08: {
        traversals.push_back(std::make_unique<C08Traversal<ParticleCell, CellFunctor>>(_dims, &cellFunctor));
        break;
      }
      case sliced: {
        traversals.push_back(std::make_unique<SlicedTraversal<ParticleCell, CellFunctor>>(_dims, &cellFunctor));
        break;
      }
      default: { AutoPasLogger->warn("Traversal type {} is not a known type!", option); }
    }
  }

  assert(traversals.size() > 0);
  return traversals;
}

template <class ParticleCell>
template <class CellFunctor>
std::unique_ptr<CellPairTraversal<ParticleCell, CellFunctor>> TraversalSelector<ParticleCell>::chooseOptimalTraversal(
    std::vector<std::unique_ptr<CellPairTraversal<ParticleCell, CellFunctor>>> &traversals) {

  // remove traversals which are not applicable to the current domain size
  traversals.erase(std::remove_if(traversals.begin(),
                                  traversals.end(),
                                  [](auto const& t) { return not t->isApplicable(); }),
                   traversals.end());
  assert(traversals.size() > 0);

  // TODO: Autotuning goes here
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

template <class ParticleCell>
template <class CellFunctor>
std::unique_ptr<CellPairTraversal<ParticleCell, CellFunctor>> TraversalSelector<ParticleCell>::getOptimalTraversal(
    CellFunctor &cellFunctor) {
  std::unique_ptr<CellPairTraversal<ParticleCell, CellFunctor>> traversal;

  if (_optimalTraversalOptions.find(typeid(CellFunctor).hash_code()) == _optimalTraversalOptions.end()) {
    auto generatedTraversals = generateTraversals<CellFunctor>(cellFunctor);
    traversal = chooseOptimalTraversal<CellFunctor>(generatedTraversals);
  } else {
    switch (_optimalTraversalOptions[typeid(CellFunctor).hash_code()]) {
      case c08: {
        traversal = std::make_unique<C08Traversal<ParticleCell, CellFunctor>>(_dims, &cellFunctor);
        break;
      }
      case sliced: {
        traversal = std::make_unique<SlicedTraversal<ParticleCell, CellFunctor>>(_dims, &cellFunctor);
        break;
      }
      default: { utils::ExceptionHandler::exception("Invalid saved optimal traversal option for this CellFunctor!"); }
    }
  }
  return traversal;
}

template <class ParticleCell>
template <class CellFunctor>
void TraversalSelector<ParticleCell>::tune(CellFunctor &cellFunctor) {
  auto generatedTraversals = generateTraversals<CellFunctor>(cellFunctor);
  chooseOptimalTraversal<CellFunctor>(generatedTraversals);
}
}  // namespace autopas