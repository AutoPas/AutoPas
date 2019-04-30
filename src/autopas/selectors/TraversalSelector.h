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
#include "../options/TraversalOption.h"
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/containers/cellPairTraversals/DummyTraversal.h"
#include "autopas/containers/cellPairTraversals/TraversalInterface.h"
#include "autopas/containers/directSum/DirectSumKokkosTraversal.h"
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
  TraversalSelector() : _dims({0, 0, 0}) {}
  /**
   * Constructor of the TraversalSelector class.
   * @param dims Array with the dimension lengths of the domain.
   */
  TraversalSelector(const std::array<unsigned long, 3> &dims) : _dims(dims) {}

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

 private:
  /**
   * indicating whether or not the optimalTraversalOption is already initialized
   */
  const std::array<unsigned long, 3> _dims;
};

template <class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3>
std::unique_ptr<CellPairTraversal<ParticleCell>> TraversalSelector<ParticleCell>::generateTraversal(
    TraversalOption traversalType, PairwiseFunctor &pairwiseFunctor) {
  switch (traversalType) {
    case TraversalOption::directSumKokkosTraversal: {
      return std::make_unique<DirectSumKokkosTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(
          &pairwiseFunctor);
    }
    case TraversalOption::directSumTraversal: {
      return std::make_unique<DirectSumTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(&pairwiseFunctor);
    }
    case TraversalOption::c08: {
      return std::make_unique<C08Traversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(_dims, &pairwiseFunctor);
    }
    case TraversalOption::sliced: {
      return std::make_unique<SlicedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(_dims,
                                                                                                  &pairwiseFunctor);
    }
    case TraversalOption::c18: {
      return std::make_unique<C18Traversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(_dims, &pairwiseFunctor);
    }
    case TraversalOption::c01: {
      return std::make_unique<C01Traversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(_dims, &pairwiseFunctor);
    }
    case TraversalOption::slicedVerlet: {
      return std::make_unique<SlicedTraversalVerlet<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(
          _dims, &pairwiseFunctor);
    }
    case TraversalOption::c18Verlet: {
      return std::make_unique<C18TraversalVerlet<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(_dims,
                                                                                                     &pairwiseFunctor);
    }
    case TraversalOption::c01Verlet: {
      return std::make_unique<C01TraversalVerlet<ParticleCell, PairwiseFunctor, useSoA, useNewton3>>(_dims,
                                                                                                     &pairwiseFunctor);
    }
    case TraversalOption::dummyTraversal: {
      return std::make_unique<DummyTraversal<ParticleCell>>(_dims);
    }
  }
  autopas::utils::ExceptionHandler::exception("Traversal type {} is not a known type!",
                                              utils::StringUtils::to_string(traversalType));
  return std::unique_ptr<CellPairTraversal<ParticleCell>>(nullptr);
}

}  // namespace autopas
