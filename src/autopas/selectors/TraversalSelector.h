/**
 * @file TraversalSelector.h
 * @author F. Gratl
 * @date 11.06.18
 */

#pragma once

#include <autopas/containers/cellPairTraversals/DirectSumTraversal.h>
#include <autopas/containers/cellPairTraversals/DummyTraversal.h>
#include <array>
#include <vector>
#include "autopas/containers/cellPairTraversals/C08Traversal.h"
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/containers/cellPairTraversals/CellPairTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/SlicedTraversal.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/Logger.h"

namespace autopas {

/**
 * Provides a way to iterate over the possible choices of TraversalOption.
 */
static std::vector<TraversalOptions> allTraversalOptions = {TraversalOptions::c08, TraversalOptions::sliced,
                                                            TraversalOptions::directSum};

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
  TraversalSelector() : _dims({0, 0, 0}), _allowedTraversalOptions({}) {}
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
   * Save the runtime of a given traversal if the functor is relevant for tuning.
   * @param pairwiseFunctor
   * @param traversal
   * @param time
   */
  template <class PairwiseFunctor>
  void addTimeMeasurement(PairwiseFunctor &pairwiseFunctor, TraversalOptions traversal, long time) {
    if (pairwiseFunctor.isRelevantForTuning()) {
      struct TimeMeasurement measurement = {traversal, time};
      _traversalTimes.push_back(measurement);
    }
  }

  /**
   * Selects the next allowed and applicable traversal.
   * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
   * @tparam useSoA
   * @tparam useNewton3
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @return Smartpointer to the selected traversal.
   */
  template <class PairwiseFunctor, bool useSoA, bool useNewton3>
  std::unique_ptr<CellPairTraversal<ParticleCell>> selectNextTraversal(PairwiseFunctor &pairwiseFunctor);

  /**
   * Selects the optimal traversal based on saved measurements.
   * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
   * @tparam useSoA
   * @tparam useNewton3
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @return Smartpointer to the selected traversal.
   */
  template <class PairwiseFunctor, bool useSoA, bool useNewton3>
  std::unique_ptr<CellPairTraversal<ParticleCell>> selectOptimalTraversal(PairwiseFunctor &pairwiseFunctor);

 private:
  template <class PairwiseFunctor, bool useSoA, bool useNewton3>
  std::vector<std::unique_ptr<CellPairTraversalInterface>> generateAllAllowedTraversals(
      PairwiseFunctor &pairwiseFunctor);

  template <class PairwiseFunctor, bool useSoA, bool useNewton3>
  std::unique_ptr<CellPairTraversal<ParticleCell>> generateTraversal(TraversalOptions traversalType,
                                                                     PairwiseFunctor &pairwiseFunctor);

  // The optimal traversal for all functors that are marked relevant.
  TraversalOptions _optimalTraversalOptions;
  // indicating whether or not the optimalTraversalOption is already initialized
  bool _isInitialized = false;
  // indicating whether we are currently testing through all options
  bool _isTuning = false;
  const std::array<unsigned long, 3> _dims;
  const std::vector<TraversalOptions> _allowedTraversalOptions;

  struct TimeMeasurement {
    TraversalOptions traversal;
    long time;
  };
  // vector of (traversal type, execution time)
  std::vector<TraversalSelector::TimeMeasurement> _traversalTimes;
};

template <class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3>
std::vector<std::unique_ptr<CellPairTraversalInterface>> TraversalSelector<ParticleCell>::generateAllAllowedTraversals(
    PairwiseFunctor &pairwiseFunctor) {
  std::vector<std::unique_ptr<CellPairTraversalInterface>> traversals;

  for (auto &option : _allowedTraversalOptions) {
    traversals.push_back(generateTraversal<PairwiseFunctor, useSoA, useNewton3>(option, pairwiseFunctor));
  }

  if (traversals.empty()) utils::ExceptionHandler::exception("TraversalSelector: No traversals were generated.");

  return traversals;
}

template <class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3>
std::unique_ptr<CellPairTraversal<ParticleCell>> TraversalSelector<ParticleCell>::generateTraversal(
    TraversalOptions traversalType, PairwiseFunctor &pairwiseFunctor) {
  std::unique_ptr<CellPairTraversal<ParticleCell>> traversal;
  switch (traversalType) {
    case TraversalOptions::directSum: {
      traversal =
          std::make_unique<DirectSumTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(&pairwiseFunctor);
      break;
    }
    case TraversalOptions::c08: {
      traversal =
          std::make_unique<C08Traversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(_dims, &pairwiseFunctor);
      break;
    }
    case TraversalOptions::sliced: {
      traversal =
          std::make_unique<SlicedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(_dims, &pairwiseFunctor);
      break;
    }
    case TraversalOptions::dummyTraversal: {
      traversal = std::make_unique<DummyTraversal<ParticleCell>>(_dims);
      break;
    }
    default: { AutoPasLog(warn, "Traversal type {} is not a known type!", traversalType); }
  }
  return traversal;
}

template <class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3>
std::unique_ptr<CellPairTraversal<ParticleCell>> TraversalSelector<ParticleCell>::selectOptimalTraversal(
    PairwiseFunctor &pairwiseFunctor) {
  // Time measure strategy
  if (_traversalTimes.empty()) {
    utils::ExceptionHandler::exception("TraversalSelector: Trying to determine fastest traversal before measuring!");
  }

  // choose the fastest traversal and reset timings
  // Initialize with something. This will be overridden.
  TraversalOptions optimalTraversalOption = TraversalOptions::c08;
  long optimalTraversalTime = std::numeric_limits<long>::max();
  AutoPasLog(debug, "TraversalSelector: Collected traversal times:");
  for (auto &&t : _traversalTimes) {
    AutoPasLog(debug, "Traversal {} took {} nanoseconds:", t.traversal, t.time);
    if (t.time < optimalTraversalTime) {
      optimalTraversalOption = t.traversal;
      optimalTraversalTime = t.time;
    }
  }
  // measurements are not needed anymore
  _traversalTimes.clear();

  // sanity check
  if (optimalTraversalTime == std::numeric_limits<long>::max()) {
    utils::ExceptionHandler::exception("TraversalSelector: Nothing was faster than max long! o_O");
  }

  // Assumption: the fastest traversal is applicable :O
  auto optimalTraversal =
      generateTraversal<PairwiseFunctor, useSoA, useNewton3>(optimalTraversalOption, pairwiseFunctor);

  _optimalTraversalOptions = optimalTraversal->getTraversalType();
  AutoPasLog(debug, "Selected traversal {}", _optimalTraversalOptions);

  return optimalTraversal;
}

template <class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3>
std::unique_ptr<CellPairTraversal<ParticleCell>> TraversalSelector<ParticleCell>::selectNextTraversal(
    PairwiseFunctor &pairwiseFunctor) {
  std::unique_ptr<CellPairTraversal<ParticleCell>> traversal;
  bool traversalIsApplicable = false;

  // choose new traversals
  while (not traversalIsApplicable) {
    // if no measurements are in yet _optimalTraversalOption is not initialized
    if (not _isTuning) {
      _optimalTraversalOptions = _allowedTraversalOptions.begin().operator*();
      _isTuning = true;
    } else {
      auto selectedTraversalIter =
          std::find(_allowedTraversalOptions.begin(), _allowedTraversalOptions.end(), _optimalTraversalOptions);
      ++selectedTraversalIter;

      // if there is no next return null
      if (selectedTraversalIter >= _allowedTraversalOptions.end()) {
        _isTuning = false;
        return std::unique_ptr<CellPairTraversal<ParticleCell>>(nullptr);
      }
      _optimalTraversalOptions = *selectedTraversalIter;
    }

    traversal = generateTraversal<PairwiseFunctor, useSoA, useNewton3>(_optimalTraversalOptions, pairwiseFunctor);
    traversalIsApplicable = traversal->isApplicable();
  }
  AutoPasLog(debug, "Testing traversal {}", _optimalTraversalOptions);

  _isInitialized = true;
  return traversal;
}

template <class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3>
std::unique_ptr<CellPairTraversal<ParticleCell>> TraversalSelector<ParticleCell>::getOptimalTraversal(
    PairwiseFunctor &pairwiseFunctor) {
  std::unique_ptr<CellPairTraversal<ParticleCell>> traversal;

  if (not _isInitialized)
    utils::ExceptionHandler::exception("TraversalSelector::getOptimalTraversal(): No Traversal selected yet!");

  traversal = generateTraversal<PairwiseFunctor, useSoA, useNewton3>(_optimalTraversalOptions, pairwiseFunctor);
  return traversal;
}
}  // namespace autopas