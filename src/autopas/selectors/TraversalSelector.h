/**
 * @file TraversalSelector.h
 * @author F. Gratl
 * @date 11.06.18
 */

#pragma once

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
static std::vector<TraversalOptions> allTraversalOptions = {TraversalOptions::c08, TraversalOptions::sliced};

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
   * Evaluates to optimal traversal based on a given functor.
   * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
   * @tparam useSoA
   * @tparam useNewton3
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @return true if still in tuning phase
   */
  template <class PairwiseFunctor, bool useSoA, bool useNewton3>
  bool tune(PairwiseFunctor &pairwiseFunctor);

  /**
   * Save the runtime of a given traversal with a given functor.
   * @param functor
   * @param traversal
   * @param time
   */
  template <class PairwiseFunctor>
  void addTimeMeasurement(PairwiseFunctor &functor, TraversalOptions traversal, long time) {
    auto functorHash = typeid(PairwiseFunctor).hash_code();
    struct TimeMeasurement measurement = {functorHash, traversal, time};
    _traversalTimes.push_back(measurement);
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

  template <class PairwiseFunctor, bool useSoA, bool useNewton3>
  std::unique_ptr<CellPairTraversal<ParticleCell>> chooseOptimalTraversal(
      std::vector<std::unique_ptr<CellPairTraversalInterface>> &traversals);

  // for each encountered cell processor save the optimal traversal. The cell processor is saved through its hash
  std::unordered_map<size_t, TraversalOptions> _optimalTraversalOptions;
  const std::array<unsigned long, 3> _dims;
  const std::vector<TraversalOptions> _allowedTraversalOptions;

  struct TimeMeasurement {
    size_t functorHash;
    TraversalOptions traversal;
    long time;
  };
  // vector of (functor hash, traversal type, execution time)
  std::vector<TraversalSelector::TimeMeasurement> _traversalTimes;
  bool _currentlyTuning = false;
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
std::unique_ptr<CellPairTraversal<ParticleCell>> TraversalSelector<ParticleCell>::chooseOptimalTraversal(
    std::vector<std::unique_ptr<CellPairTraversalInterface>> &traversals) {
  // remove traversals which are not applicable to the current domain size
  traversals.erase(
      std::remove_if(traversals.begin(), traversals.end(), [](auto const &t) { return not t->isApplicable(); }),
      traversals.end());
  if (traversals.empty())
    utils::ExceptionHandler::exception("TraversalSelector: None of the allowed traversals were applicable.");

  auto functorHash = typeid(PairwiseFunctor).hash_code();
  size_t bestTraversal = 0;
  // Test all options to find the fastest
  // If functor is unknown or no measurements were made by now we start with the first traversal
  if (_traversalTimes.empty() or _optimalTraversalOptions.find(functorHash) == _optimalTraversalOptions.end()) {
    bestTraversal = 0;
  } else if (_currentlyTuning) {
    // if we are in tuning state just select next traversal
    bestTraversal = std::find_if(traversals.begin(), traversals.end(),
                                 [&](const std::unique_ptr<CellPairTraversalInterface> &a) {
                                   return a->getTraversalType() == _optimalTraversalOptions[functorHash];
                                 }) -
                    traversals.begin();
    ++bestTraversal;
    // if the last possible traversal has already been tested choose fastest one and reset timings
    if (bestTraversal >= traversals.size()) {
      _currentlyTuning = false;
      TraversalOptions fastestTraversal;
      long fastestTime = std::numeric_limits<long>::max();
      AutoPasLog(debug, "TraversalSelector: Collected traversals:");
      for (auto t = _traversalTimes.begin(); t != _traversalTimes.end();) {
        // only look at times from the correct functor
        if (t->functorHash != functorHash) {
          ++t;
          continue;
        }
        AutoPasLog(debug, "Traversal {} took {} nanoseconds:", t->traversal, t->time);
        if (t->time < fastestTime) {
          fastestTraversal = t->traversal;
          fastestTime = t->time;
        }
        // when a measurement was used delete it
        t = _traversalTimes.erase(t);
      }
      // sanity check
      if (fastestTime == std::numeric_limits<long>::max()) {
        utils::ExceptionHandler::exception("TraversalSelector: nothing was faster than max long oO");
      }

      // find id of fastest traversal in passed traversal list
      bestTraversal = std::find_if(traversals.begin(), traversals.end(),
                                   [&](const std::unique_ptr<CellPairTraversalInterface> &a) {
                                     return a->getTraversalType() == fastestTraversal;
                                   }) -
                      traversals.begin();
    }
  }

  // Tedious downcast
  std::unique_ptr<CellPairTraversal<ParticleCell>> optimalTraversal(
      dynamic_cast<CellPairTraversal<ParticleCell> *>(traversals[bestTraversal].release()));

  _optimalTraversalOptions[functorHash] = optimalTraversal->getTraversalType();
  AutoPasLog(debug, "{} traversal {}", _currentlyTuning ? "Testing" : "Selected", optimalTraversal->getTraversalType());

  return optimalTraversal;
}

template <class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3>
std::unique_ptr<CellPairTraversal<ParticleCell>> TraversalSelector<ParticleCell>::selectOptimalTraversal(
    PairwiseFunctor &pairwiseFunctor) {
  auto functorHash = typeid(PairwiseFunctor).hash_code();

  // Time measure strategy
  if (_traversalTimes.empty()) {
    utils::ExceptionHandler::exception("TraversalSelector: Trying to determine fastest traversal before measuring!");
  }

  // choose the fastest traversal and reset timings
  TraversalOptions optimalTraversalOption;
  long optimalTraversalTime = std::numeric_limits<long>::max();
  AutoPasLog(debug, "TraversalSelector: Collected traversal times:");
  //  for (auto t = _traversalTimes.begin(); t != _traversalTimes.end(); ++t) {
  for (auto &&t : _traversalTimes) {
    // only look at times from the correct functor
    if (t.functorHash != functorHash) {
      continue;
    }
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

  _optimalTraversalOptions[functorHash] = optimalTraversal->getTraversalType();
  AutoPasLog(debug, "Selected traversal {}", _optimalTraversalOptions[functorHash]);

  return optimalTraversal;
}

template <class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3>
std::unique_ptr<CellPairTraversal<ParticleCell>> TraversalSelector<ParticleCell>::selectNextTraversal(
    PairwiseFunctor &pairwiseFunctor) {
  std::unique_ptr<CellPairTraversal<ParticleCell>> traversal;
  bool traversalIsApplicable = false;

  auto functorHash = typeid(PairwiseFunctor).hash_code();

  // choose new traversals
  while (not traversalIsApplicable) {
    // if the functor is not known yet choose the first allowed traversal
    if (_optimalTraversalOptions.find(functorHash) == _optimalTraversalOptions.end()) {
      _optimalTraversalOptions[functorHash] = _allowedTraversalOptions.begin().operator*();
    } else {
      auto selectedTraversalIter = std::find(_allowedTraversalOptions.begin(), _allowedTraversalOptions.end(),
                                             _optimalTraversalOptions[functorHash]);
      ++selectedTraversalIter;

      // if there is no next return null
      if (selectedTraversalIter >= _allowedTraversalOptions.end()) {
        return std::unique_ptr<CellPairTraversal<ParticleCell>>(nullptr);
      }
      _optimalTraversalOptions[functorHash] = *selectedTraversalIter;
    }

    traversal =
        generateTraversal<PairwiseFunctor, useSoA, useNewton3>(_optimalTraversalOptions[functorHash], pairwiseFunctor);
    traversalIsApplicable = traversal->isApplicable();
  }
  AutoPasLog(debug, "Testing traversal {}", _optimalTraversalOptions[functorHash]);

  return traversal;
}

template <class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3>
std::unique_ptr<CellPairTraversal<ParticleCell>> TraversalSelector<ParticleCell>::getOptimalTraversal(
    PairwiseFunctor &pairwiseFunctor) {
  std::unique_ptr<CellPairTraversal<ParticleCell>> traversal;

  auto functorHash = typeid(PairwiseFunctor).hash_code();

  if (_optimalTraversalOptions.find(functorHash) == _optimalTraversalOptions.end()) {
    utils::ExceptionHandler::exception("TraversalSelector: Traversal not yet selected!");
  } else {
    traversal =
        generateTraversal<PairwiseFunctor, useSoA, useNewton3>(_optimalTraversalOptions[functorHash], pairwiseFunctor);
  }
  return traversal;
}

template <class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3>
bool TraversalSelector<ParticleCell>::tune(PairwiseFunctor &pairwiseFunctor) {
  _currentlyTuning = true;
  // Workaround for Containers that do not use traversals. If there are no traversals there is nothing to tune.
  if (_allowedTraversalOptions.empty()) return false;
  auto generatedTraversals = generateAllAllowedTraversals<PairwiseFunctor, useSoA, useNewton3>(pairwiseFunctor);
  chooseOptimalTraversal<PairwiseFunctor, useSoA, useNewton3>(generatedTraversals);

  return _currentlyTuning;
}
}  // namespace autopas