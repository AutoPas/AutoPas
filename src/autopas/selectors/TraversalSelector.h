/**
 * @file TraversalSelector.h
 * @author F. Gratl
 * @date 11.06.18
 */

#pragma once

#include <array>
#include <numeric>
#include <unordered_map>
#include <vector>
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/containers/cellPairTraversals/DummyTraversal.h"
#include "autopas/containers/cellPairTraversals/TraversalInterface.h"
#include "autopas/containers/directSum/DirectSumTraversal.h"
#include "autopas/containers/linkedCells/traversals/C01Traversal.h"
#include "autopas/containers/linkedCells/traversals/C08Traversal.h"
#include "autopas/containers/linkedCells/traversals/C18Traversal.h"
#include "autopas/containers/linkedCells/traversals/SlicedTraversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/C01TraversalVerlet.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/C18TraversalVerlet.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/SlicedTraversalVerlet.h"
#include "autopas/options/SelectorStrategie.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/Logger.h"
#include "autopas/utils/StringUtils.h"
#include "autopas/utils/TrivialHash.h"

namespace autopas {

/**
 * Selector for a container traversal.
 * @tparam ParticleCell
 */
template <class ParticleCell>
class TraversalSelector {
 public:
  /**
   * Dummy constructor such that this class can be used in maps
   */
  TraversalSelector() : _currentTraversal(TraversalOption(-1)), _dims({0, 0, 0}), _allowedTraversalOptions({}) {}
  /**
   * Constructor of the TraversalSelector class.
   * @param dims Array with the dimension lengths of the domain.
   * @param allowedTraversalOptions Vector of traversals the selector can choose from.
   */
  TraversalSelector(const std::array<unsigned long, 3> &dims,
                    const std::vector<TraversalOption> &allowedTraversalOptions)
      : _currentTraversal(TraversalOption(-1)), _dims(dims), _allowedTraversalOptions(allowedTraversalOptions) {}

  /**
   * Gets the optimal traversal for a given cell functor. If no traversal is selected yet a optimum search is started.
   * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
   * @tparam useSoA
   * @tparam useNewton3
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @return Smartpointer to the optimal traversal.
   */
  template <class PairwiseFunctor, bool useSoA, bool useNewton3>
  std::unique_ptr<CellPairTraversal<ParticleCell>> getCurrentTraversal(PairwiseFunctor &pairwiseFunctor);

  /**
   * Generates a given Traversal for the given properties.
   * @tparam PairwiseFunctor
   * @tparam useSoA
   * @tparam useNewton3
   * @param traversalType
   * @param pairwiseFunctor
   * @return Smartpointer to the traversal.
   */
  template <class PairwiseFunctor, bool useSoA, bool useNewton3>
  std::unique_ptr<CellPairTraversal<ParticleCell>> generateTraversal(TraversalOption traversalType,
                                                                     PairwiseFunctor &pairwiseFunctor);

  /**
   * Sets the traversal to the given Option.
   * @param traversalOption
   */
  void selectTraversal(TraversalOption traversalOption);

 private:
  /**
   * The optimal traversal for all functors that are marked relevant.#
   */
  TraversalOption _currentTraversal;
  /**
   * indicating whether or not the optimalTraversalOption is already initialized
   */
  const std::array<unsigned long, 3> _dims;
  const std::vector<TraversalOption> _allowedTraversalOptions;

 public:
  const std::vector<TraversalOption> &getAllowedTraversalOptions() const;

 private:
  struct TimeMeasurement {
    TraversalOption traversal;
    long time;
  };
  // vector of (traversal type, execution time)
  std::vector<TraversalSelector::TimeMeasurement> _traversalTimes;
};

template <class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3>
std::unique_ptr<CellPairTraversal<ParticleCell>> TraversalSelector<ParticleCell>::generateTraversal(
    TraversalOption traversalType, PairwiseFunctor &pairwiseFunctor) {
  std::unique_ptr<CellPairTraversal<ParticleCell>> traversal;
  switch (traversalType) {
    case TraversalOption::directSumTraversal: {
      traversal =
          std::make_unique<DirectSumTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(&pairwiseFunctor);
      break;
    }
    case TraversalOption::c08: {
      traversal =
          std::make_unique<C08Traversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(_dims, &pairwiseFunctor);
      break;
    }
    case TraversalOption::sliced: {
      traversal =
          std::make_unique<SlicedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(_dims, &pairwiseFunctor);
      break;
    }
    case TraversalOption::c18: {
      traversal =
          std::make_unique<C18Traversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(_dims, &pairwiseFunctor);
      break;
    }
    case TraversalOption::c01: {
      traversal =
          std::make_unique<C01Traversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(_dims, &pairwiseFunctor);
      break;
    }
    case TraversalOption::slicedVerlet: {
      traversal = std::make_unique<SlicedTraversalVerlet<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(
          _dims, &pairwiseFunctor);
      break;
    }
    case TraversalOption::c18Verlet: {
      traversal = std::make_unique<C18TraversalVerlet<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(
          _dims, &pairwiseFunctor);
      break;
    }
    case TraversalOption::c01Verlet: {
      traversal = std::make_unique<C01TraversalVerlet<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(
          _dims, &pairwiseFunctor);
      break;
    }
    case TraversalOption::dummyTraversal: {
      traversal = std::make_unique<DummyTraversal<ParticleCell>>(_dims);
      break;
    }
    default: {
      AutoPasLog(warn, "Traversal type {} is not a known type!", utils::StringUtils::to_string(traversalType));
    }
  }
  return traversal;
}

template <class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3>
std::unique_ptr<CellPairTraversal<ParticleCell>> TraversalSelector<ParticleCell>::getCurrentTraversal(
    PairwiseFunctor &pairwiseFunctor) {
  std::unique_ptr<CellPairTraversal<ParticleCell>> traversal;

  if (_currentTraversal == TraversalOption(-1))
    utils::ExceptionHandler::exception(
        "TraversalSelector::getCurrentTraversal(): Traversal selector not initialized (nothing selected yet)!");

  traversal = generateTraversal<PairwiseFunctor, useSoA, useNewton3>(_currentTraversal, pairwiseFunctor);
  return traversal;
}

template <class ParticleCell>
void TraversalSelector<ParticleCell>::selectTraversal(TraversalOption traversalOption) {
  _currentTraversal = traversalOption;
}
template <class ParticleCell>
const std::vector<TraversalOption> &TraversalSelector<ParticleCell>::getAllowedTraversalOptions() const {
  return _allowedTraversalOptions;
}

}  // namespace autopas
