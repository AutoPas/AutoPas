/**
 * @file TraversalSelector.h
 * @author F. Gratl
 * @date 11.06.18
 */

#pragma once

#include <autopas/containers/cellPairTraversals/DirectSumTraversal.h>
#include <autopas/containers/cellPairTraversals/DummyTraversal.h>
#include <array>
#include <numeric>
#include <unordered_map>
#include <vector>
#include "autopas/containers/cellPairTraversals/C08Traversal.h"
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/containers/cellPairTraversals/CellPairTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/DummyTraversal.h"
#include "autopas/containers/cellPairTraversals/SlicedTraversal.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/Logger.h"
#include "autopas/utils/TrivialHash.h"

namespace autopas {

/**
 * Possible choices for the employed selectors.
 */
enum SelectorStrategy {
  /**
   * Fastest absolute value.
   */
  fastestAbs,
  /**
   * Fastest mean value.
   */
  fastestMean,
  /**
   * Fastest median value
   */
  fastestMedian
};

/**
 * Provides a way to iterate over the possible choices of TraversalOption.
 */
static std::vector<TraversalOptions> allTraversalOptions = {TraversalOptions::c08, TraversalOptions::sliced,
                                                            TraversalOptions::directSumTraversal};

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
   * @param strategy Strategy the selector should employ to choose the best traversal.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @return Smartpointer to the selected traversal.
   */
  template <class PairwiseFunctor, bool useSoA, bool useNewton3>
  std::unique_ptr<CellPairTraversal<ParticleCell>> selectOptimalTraversal(SelectorStrategy strategy,
                                                                          PairwiseFunctor &pairwiseFunctor);

 private:
  void findFastestAbsTraversal();

  void findFastestMeanTraversal();

  void findFastestMedianTraversal();

  template <class PairwiseFunctor, bool useSoA, bool useNewton3>
  std::vector<std::unique_ptr<CellPairTraversalInterface>> generateAllAllowedTraversals(
      PairwiseFunctor &pairwiseFunctor);

  template <class PairwiseFunctor, bool useSoA, bool useNewton3>
  std::unique_ptr<CellPairTraversal<ParticleCell>> generateTraversal(TraversalOptions traversalType,
                                                                     PairwiseFunctor &pairwiseFunctor);

  // The optimal traversal for all functors that are marked relevant.
  TraversalOptions _optimalTraversalOption;
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
    case TraversalOptions::directSumTraversal: {
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
    SelectorStrategy strategy, PairwiseFunctor &pairwiseFunctor) {
  // Time measure strategy
  if (_traversalTimes.empty()) {
    utils::ExceptionHandler::exception("TraversalSelector: Trying to determine fastest traversal before measuring!");
  }

  switch (strategy) {
    case SelectorStrategy::fastestAbs: {
      findFastestAbsTraversal();
      break;
    }
    case SelectorStrategy::fastestMean: {
      findFastestMeanTraversal();
      break;
    }
    case SelectorStrategy::fastestMedian: {
      findFastestMedianTraversal();
      break;
    }
    default:
      utils::ExceptionHandler::exception("TraversalSelector: Unknown selector strategy {}", strategy);
  }

  // measurements are not needed anymore
  _traversalTimes.clear();

  // Assumption: the fastest traversal is applicable :O
  auto traversal = generateTraversal<PairwiseFunctor, useSoA, useNewton3>(_optimalTraversalOption, pairwiseFunctor);

  AutoPasLog(debug, "Selected traversal {}", _optimalTraversalOption);
  return traversal;
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
      _optimalTraversalOption = _allowedTraversalOptions.begin().operator*();
      _isTuning = true;
    } else {
      auto selectedTraversalIter =
          std::find(_allowedTraversalOptions.begin(), _allowedTraversalOptions.end(), _optimalTraversalOption);
      ++selectedTraversalIter;

      // if there is no next return null
      if (selectedTraversalIter >= _allowedTraversalOptions.end()) {
        _isTuning = false;
        return std::unique_ptr<CellPairTraversal<ParticleCell>>(nullptr);
      }
      _optimalTraversalOption = *selectedTraversalIter;
    }

    traversal = generateTraversal<PairwiseFunctor, useSoA, useNewton3>(_optimalTraversalOption, pairwiseFunctor);
    traversalIsApplicable = traversal->isApplicable();
  }
  AutoPasLog(debug, "Testing traversal {}", _optimalTraversalOption);

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

  traversal = generateTraversal<PairwiseFunctor, useSoA, useNewton3>(_optimalTraversalOption, pairwiseFunctor);
  return traversal;
}

template <class ParticleCell>
void TraversalSelector<ParticleCell>::findFastestAbsTraversal() {
  // choose the fastest traversal and reset timings
  // Initialize with something. This will be overridden.
  long optimalTraversalTime = std::numeric_limits<long>::max();
  AutoPasLog(debug, "TraversalSelector: Collected traversal times:");
  for (auto &&t : _traversalTimes) {
    AutoPasLog(debug, "Traversal {} took {} nanoseconds.", t.traversal, t.time);
    if (t.time < optimalTraversalTime) {
      _optimalTraversalOption = t.traversal;
      optimalTraversalTime = t.time;
    }
  }

  // sanity check
  if (optimalTraversalTime == std::numeric_limits<long>::max()) {
    utils::ExceptionHandler::exception("TraversalSelector: Nothing was faster than max long! o_O");
  }
}

template <class ParticleCell>
void TraversalSelector<ParticleCell>::findFastestMeanTraversal() {
  // choose the fastest traversal and reset timings
  // reorder measurements
  std::unordered_map<TraversalOptions, std::vector<long>, TrivialHash> measurementsMap;
  AutoPasLog(debug, "TraversalSelector: Collected traversal times:");
  for (auto &&t : _traversalTimes) {
    AutoPasLog(debug, "Traversal {} took {} nanoseconds", t.traversal, t.time);
    measurementsMap[t.traversal].push_back(t.time);
  }

  long optimalTraversalTime = std::numeric_limits<long>::max();
  // @todo: when verlet list traversals are here apply weights to measurement w/ or w/o vl rebuild
  for (auto &&m : measurementsMap) {
    long meanTime = std::accumulate(m.second.begin(), m.second.end(), 0l) / m.second.size();
    AutoPasLog(debug, "Traversal {} mean: {} nanoseconds", m.first, meanTime);
    if (meanTime < optimalTraversalTime) {
      optimalTraversalTime = meanTime;
      _optimalTraversalOption = m.first;
    }
  }

  // sanity check
  if (optimalTraversalTime == std::numeric_limits<long>::max()) {
    utils::ExceptionHandler::exception("TraversalSelector: Nothing was faster than max long! o_O");
  }
}

template <class ParticleCell>
void TraversalSelector<ParticleCell>::findFastestMedianTraversal() {
  // choose the fastest traversal and reset timings
  // reorder measurements
  std::unordered_map<TraversalOptions, std::vector<long>, TrivialHash> measurementsMap;
  AutoPasLog(debug, "TraversalSelector: Collected traversal times:");
  for (auto &&t : _traversalTimes) {
    AutoPasLog(debug, "Traversal {} took {} nanoseconds", t.traversal, t.time);
    measurementsMap[t.traversal].push_back(t.time);
  }

  long optimalTraversalTime = std::numeric_limits<long>::max();
  for (auto &&m : measurementsMap) {
    std::sort(m.second.begin(), m.second.end());
    long medianTime = m.second[m.second.size() / 2];
    AutoPasLog(debug, "Traversal {} median: {} nanoseconds", m.first, medianTime);
    if (medianTime < optimalTraversalTime) {
      optimalTraversalTime = medianTime;
      _optimalTraversalOption = m.first;
    }
  }

  // sanity check
  if (optimalTraversalTime == std::numeric_limits<long>::max()) {
    utils::ExceptionHandler::exception("TraversalSelector: Nothing was faster than max long! o_O");
  }
}

}  // namespace autopas
