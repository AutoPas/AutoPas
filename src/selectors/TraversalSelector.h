/**
 * @file TraversalSelector.h
 * @author F. Gratl
 * @date 11.06.18
 */

#pragma once

#include <containers/cellPairTraversals/C08Traversal.h>
#include <containers/cellPairTraversals/CellPairTraversal.h>
#include <containers/cellPairTraversals/SlicedTraversal.h>
#include <pairwiseFunctors/CellFunctor.h>
#include <array>
#include <vector>
#include "utils/ExceptionHandler.h"
#include "utils/Logger.h"

namespace autopas {

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
   * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
   * @tparam useSoA
   * @tparam useNewton3
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @return Smartpointer to the optimal traversal.
   */
  template <class PairwiseFunctor, bool useSoA, bool useNewton3>
  std::unique_ptr<CellPairTraversal<ParticleCell>> getOptimalTraversal(PairwiseFunctor &pairwiseFunctor);

  /**
   * Evaluates to optimal traversal based on a given cell functor.
   * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
   * @tparam useSoA
   * @tparam useNewton3
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  template <class PairwiseFunctor, bool useSoA, bool useNewton3>
  void tune(PairwiseFunctor &pairwiseFunctor);

 private:
  template <class PairwiseFunctor, bool useSoA, bool useNewton3>
  std::vector<std::unique_ptr<CellPairTraversalInterface>> generateTraversals(PairwiseFunctor &pairwiseFunctor);

  template <class PairwiseFunctor, bool useSoA, bool useNewton3>
  std::unique_ptr<CellPairTraversal<ParticleCell>> chooseOptimalTraversal(
      std::vector<std::unique_ptr<CellPairTraversalInterface>> &traversals);

  // for each encountered cell processor save the optimal traversal. The cell processor is saved through its hash
  std::unordered_map<size_t, TraversalOptions> _optimalTraversalOptions;
  const std::array<unsigned long, 3> _dims;
  const std::vector<TraversalOptions> _allowedTraversalOptions;
};

template <class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3>
std::vector<std::unique_ptr<CellPairTraversalInterface>> TraversalSelector<ParticleCell>::generateTraversals(
    PairwiseFunctor &pairwiseFunctor) {
  std::vector<std::unique_ptr<CellPairTraversalInterface>> traversals;

  for (auto &option : _allowedTraversalOptions) {
    switch (option) {
      case TraversalOptions::c08: {
        traversals.push_back(
            std::make_unique<C08Traversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(_dims, &pairwiseFunctor));
        break;
      }
      case TraversalOptions::sliced: {
        traversals.push_back(std::make_unique<SlicedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(
            _dims, &pairwiseFunctor));
        break;
      }
      default: { AutoPasLogger->warn("Traversal type {} is not a known type!", option); }
    }
  }

  if (traversals.empty())
    utils::ExceptionHandler::exception("TraversalSelector: No traversals were generated.");

  return traversals;
}

template <class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3>
std::unique_ptr<CellPairTraversal<ParticleCell>> TraversalSelector<ParticleCell>::chooseOptimalTraversal(
    std::vector<std::unique_ptr<CellPairTraversalInterface>> &traversals) {
  // remove traversals which are not applicable to the current domain size
  traversals.erase(
      std::remove_if(traversals.begin(), traversals.end(), [](auto const &t) { return not t->isApplicable(); }),
      traversals.end());
  if (traversals.empty())
    utils::ExceptionHandler::exception("TraversalSelector: None of the allowed traversals were applicable.");

  // TODO: Autotuning goes here
  // Tedious downcast
  std::unique_ptr<CellPairTraversal<ParticleCell>> optimalTraversal(
      dynamic_cast<CellPairTraversal<ParticleCell> *>(traversals.front().release()));

  _optimalTraversalOptions[typeid(PairwiseFunctor).hash_code()] = optimalTraversal->getTraversalType();

  return optimalTraversal;
}

template <class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3>
std::unique_ptr<CellPairTraversal<ParticleCell>> TraversalSelector<ParticleCell>::getOptimalTraversal(
    PairwiseFunctor &pairwiseFunctor) {
  std::unique_ptr<CellPairTraversal<ParticleCell>> traversal;

  if (_optimalTraversalOptions.find(typeid(PairwiseFunctor).hash_code()) == _optimalTraversalOptions.end()) {
    auto generatedTraversals = generateTraversals<PairwiseFunctor, useSoA, useNewton3>(pairwiseFunctor);
    traversal = chooseOptimalTraversal<PairwiseFunctor, useSoA, useNewton3>(generatedTraversals);
  } else {
    switch (_optimalTraversalOptions[typeid(PairwiseFunctor).hash_code()]) {
      case TraversalOptions::c08: {
        traversal =
            std::make_unique<C08Traversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(_dims, &pairwiseFunctor);
        break;
      }
      case TraversalOptions::sliced: {
        traversal = std::make_unique<SlicedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(
            _dims, &pairwiseFunctor);
        break;
      }
      default: { utils::ExceptionHandler::exception("Invalid saved optimal traversal option for this CellFunctor!"); }
    }
  }
  return traversal;
}

template <class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3>
void TraversalSelector<ParticleCell>::tune(PairwiseFunctor &pairwiseFunctor) {
  // Workaround for Containers that do not use traversals.
  if (_allowedTraversalOptions.empty()) return;
  auto generatedTraversals = generateTraversals<PairwiseFunctor, useSoA, useNewton3>(pairwiseFunctor);
  chooseOptimalTraversal<PairwiseFunctor, useSoA, useNewton3>(generatedTraversals);
}

}  // namespace autopas