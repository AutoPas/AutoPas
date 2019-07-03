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
#include "TraversalSelectorInfo.h"
#include "autopas/containers/adaptiveLinkedCells/traversals/C01TraversalAdaptive.h"
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/containers/cellPairTraversals/DummyTraversal.h"
#include "autopas/containers/cellPairTraversals/TraversalInterface.h"
#include "autopas/containers/directSum/DirectSumTraversal.h"
#include "autopas/containers/linkedCells/traversals/C01CudaTraversal.h"
#include "autopas/containers/linkedCells/traversals/C01Traversal.h"
#include "autopas/containers/linkedCells/traversals/C08Traversal.h"
#include "autopas/containers/linkedCells/traversals/C18Traversal.h"
#include "autopas/containers/linkedCells/traversals/SlicedTraversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/TraversalVerlet.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/C01TraversalVerlet.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/C18TraversalVerlet.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/SlicedTraversalVerlet.h"
#include "autopas/options/SelectorStrategyOption.h"
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
   * Generates a given Traversal for the given properties.
   * @tparam PairwiseFunctor
   * @tparam useSoA
   * @tparam useNewton3
   * @param traversalType
   * @param pairwiseFunctor
   * @param info
   * @return Smartpointer to the traversal.
   */
  template <class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
  static std::unique_ptr<CellPairTraversal<ParticleCell, dataLayout, useNewton3>> generateTraversal(
      TraversalOption traversalType, PairwiseFunctor &pairwiseFunctor, const TraversalSelectorInfo<ParticleCell> &info);

  /**
   * Generates a given Traversal for the given properties. Requires less templates but only returns a TraversalInterface
   * smart pointer.
   * @tparam PairwiseFunctor
   * @param traversalType
   * @param pairwiseFunctor
   * @param info
   * @param dataLayout
   * @param useNewton3
   * @return Smartpointer to the traversal.
   */
  template <class PairwiseFunctor>
  static std::unique_ptr<TraversalInterface> generateTraversal(TraversalOption traversalType,
                                                               PairwiseFunctor &pairwiseFunctor,
                                                               const TraversalSelectorInfo<ParticleCell> &info,
                                                               DataLayoutOption dataLayout, Newton3Option useNewton3);
};

template <class ParticleCell>
template <class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
std::unique_ptr<CellPairTraversal<ParticleCell, dataLayout, useNewton3>>
TraversalSelector<ParticleCell>::generateTraversal(TraversalOption traversalType, PairwiseFunctor &pairwiseFunctor,
                                                   const TraversalSelectorInfo<ParticleCell> &info) {
  switch (traversalType) {
    // Direct sum
    case TraversalOption::directSumTraversal: {
      return std::make_unique<DirectSumTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          &pairwiseFunctor);
    }
    // Linked cell
    case TraversalOption::c08: {
      return std::make_unique<C08Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.dims, &pairwiseFunctor, info.cutoff, info.cellLength);
    }
    case TraversalOption::sliced: {
      return std::make_unique<SlicedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.dims, &pairwiseFunctor, info.cutoff, info.cellLength);
    }
    case TraversalOption::c18: {
      return std::make_unique<C18Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.dims, &pairwiseFunctor, info.cutoff, info.cellLength);
    }
    case TraversalOption::c01: {
      return std::make_unique<C01Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.dims, &pairwiseFunctor, info.cutoff, info.cellLength);
    }
    case TraversalOption::c01CombinedSoA: {
      return std::make_unique<C01Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, true>>(
          info.dims, &pairwiseFunctor, info.cutoff, info.cellLength);
    }
    // Adaptive linked cell
    case TraversalOption::c01Adaptive: {
      return std::make_unique<C01TraversalAdaptive<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.dims, &pairwiseFunctor, info.cutoff, info.cellLength);
    }
    // Verlet
    case TraversalOption::slicedVerlet: {
      return std::make_unique<SlicedTraversalVerlet<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.dims, &pairwiseFunctor);
    }
    case TraversalOption::c18Verlet: {
      return std::make_unique<C18TraversalVerlet<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.dims, &pairwiseFunctor);
    }
    case TraversalOption::c01Verlet: {
      return std::make_unique<C01TraversalVerlet<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.dims, &pairwiseFunctor);
    }
    case TraversalOption::c01Cuda: {
      return std::make_unique<C01CudaTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.dims, &pairwiseFunctor);
    }
    case TraversalOption::verletTraversal: {
      return std::make_unique<TraversalVerlet<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(info.dims,
                                                                                                      &pairwiseFunctor);
    }
    case TraversalOption::dummyTraversal: {
      return std::make_unique<DummyTraversal<ParticleCell, dataLayout, useNewton3>>(info.dims);
    }
  }
  autopas::utils::ExceptionHandler::exception("Traversal type {} is not a known type!",
                                              utils::StringUtils::to_string(traversalType));
  return std::unique_ptr<CellPairTraversal<ParticleCell, dataLayout, useNewton3>>(nullptr);
}
template <class ParticleCell>
template <class PairwiseFunctor>
std::unique_ptr<TraversalInterface> TraversalSelector<ParticleCell>::generateTraversal(
    TraversalOption traversalType, PairwiseFunctor &pairwiseFunctor,
    const TraversalSelectorInfo<ParticleCell> &traversalInfo, DataLayoutOption dataLayout, Newton3Option newton3) {
  switch (dataLayout) {
    case DataLayoutOption::aos: {
      if (newton3 == Newton3Option::enabled) {
        return TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayoutOption::aos,
                                                                           true>(traversalType, pairwiseFunctor,
                                                                                 traversalInfo);
      } else {
        return TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayoutOption::aos,
                                                                           false>(traversalType, pairwiseFunctor,
                                                                                  traversalInfo);
      }
    }
    case DataLayoutOption::soa: {
      if (newton3 == Newton3Option::enabled) {
        return TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayoutOption::soa,
                                                                           true>(traversalType, pairwiseFunctor,
                                                                                 traversalInfo);
      } else {
        return TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayoutOption::soa,
                                                                           false>(traversalType, pairwiseFunctor,
                                                                                  traversalInfo);
      }
    }
    case DataLayoutOption::cuda: {
      if (newton3 == Newton3Option::enabled) {
        return TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayoutOption::cuda,
                                                                           true>(traversalType, pairwiseFunctor,
                                                                                 traversalInfo);
      } else {
        return TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayoutOption::cuda,
                                                                           false>(traversalType, pairwiseFunctor,
                                                                                  traversalInfo);
      }
    }
  }

  autopas::utils::ExceptionHandler::exception("Traversal type {} is not a known type!",
                                              utils::StringUtils::to_string(traversalType));
  return std::unique_ptr<TraversalInterface>(nullptr);
}
}  // namespace autopas
