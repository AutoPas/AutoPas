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
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VarVerletTraversalAsBuild.h"
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
   * Generates a given Traversal for the given properties.
   * @tparam PairwiseFunctor
   * @tparam useSoA
   * @tparam useNewton3
   * @param traversalType
   * @param pairwiseFunctor
   * @param info
   * @return Smartpointer to the traversal.
   */
  template <class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
  static std::unique_ptr<CellPairTraversal<ParticleCell>> generateTraversal(
      TraversalOption traversalType, PairwiseFunctor& pairwiseFunctor, const TraversalSelectorInfo<ParticleCell>& info);
};

template <class ParticleCell>
template <class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
std::unique_ptr<CellPairTraversal<ParticleCell>> TraversalSelector<ParticleCell>::generateTraversal(
    TraversalOption traversalType, PairwiseFunctor& pairwiseFunctor, const TraversalSelectorInfo<ParticleCell>& info) {
  switch (traversalType) {
    // Direct sum
    case TraversalOption::directSumTraversal: {
      return std::make_unique<DirectSumTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>>(
          &pairwiseFunctor);
    }
    // Linked cell
    case TraversalOption::c08: {
      return std::make_unique<C08Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>>(
          info.dims, &pairwiseFunctor, info.cutoff, info.cellLength);
    }
    case TraversalOption::sliced: {
      return std::make_unique<SlicedTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>>(
          info.dims, &pairwiseFunctor, info.cutoff, info.cellLength);
    }
    case TraversalOption::c18: {
      return std::make_unique<C18Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>>(
          info.dims, &pairwiseFunctor, info.cutoff, info.cellLength);
    }
    case TraversalOption::c01: {
      return std::make_unique<C01Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>>(
          info.dims, &pairwiseFunctor, info.cutoff, info.cellLength);
    }
    // Verlet
    case TraversalOption::slicedVerlet: {
      return std::make_unique<SlicedTraversalVerlet<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>>(
          info.dims, &pairwiseFunctor);
    }
    case TraversalOption::c18Verlet: {
      return std::make_unique<C18TraversalVerlet<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>>(
          info.dims, &pairwiseFunctor);
    }
    case TraversalOption::c01Verlet: {
      return std::make_unique<C01TraversalVerlet<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>>(
          info.dims, &pairwiseFunctor);
    }
    case TraversalOption::c01Cuda: {
      return std::make_unique<C01CudaTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>>(
          info.dims, &pairwiseFunctor);
    }
    case TraversalOption::verletTraversal: {
      return std::make_unique<TraversalVerlet<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>>(info.dims,
                                                                                                      &pairwiseFunctor);
    }
    case TraversalOption::varVerletTraversalAsBuild: {
      return std::make_unique<
          VarVerletTraversalAsBuild<ParticleCell, typename ParticleCell::ParticleType, PairwiseFunctor, useNewton3>>(
          &pairwiseFunctor);
    }
    case TraversalOption::dummyTraversal: {
      return std::make_unique<DummyTraversal<ParticleCell>>(info.dims);
    }
  }
  autopas::utils::ExceptionHandler::exception("Traversal type {} is not a known type!",
                                              utils::StringUtils::to_string(traversalType));
  return std::unique_ptr<CellPairTraversal<ParticleCell>>(nullptr);
}
}  // namespace autopas
