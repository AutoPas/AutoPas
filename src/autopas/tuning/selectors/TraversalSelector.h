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

#include "autopas/containers/TraversalInterface.h"
#include "autopas/containers/directSum/traversals/DSSequentialTraversal.h"
#include "autopas/containers/directSum/traversals/DSSequentialTraversal3B.h"
#include "autopas/containers/linkedCells/traversals/LCC01Traversal.h"
#include "autopas/containers/linkedCells/traversals/LCC01Traversal3B.h"
#include "autopas/containers/linkedCells/traversals/LCC04CombinedSoATraversal.h"
#include "autopas/containers/linkedCells/traversals/LCC04HCPTraversal.h"
#include "autopas/containers/linkedCells/traversals/LCC04Traversal.h"
#include "autopas/containers/linkedCells/traversals/LCC08Traversal.h"
#include "autopas/containers/linkedCells/traversals/LCC18Traversal.h"
#include "autopas/containers/linkedCells/traversals/LCSlicedBalancedTraversal.h"
#include "autopas/containers/linkedCells/traversals/LCSlicedC02Traversal.h"
#include "autopas/containers/linkedCells/traversals/LCSlicedTraversal.h"
#include "autopas/containers/octree/traversals/OTC01Traversal.h"
#include "autopas/containers/octree/traversals/OTC18Traversal.h"
#include "autopas/containers/verletClusterLists/traversals/VCLC01BalancedTraversal.h"
#include "autopas/containers/verletClusterLists/traversals/VCLC06Traversal.h"
#include "autopas/containers/verletClusterLists/traversals/VCLClusterIterationTraversal.h"
#include "autopas/containers/verletClusterLists/traversals/VCLSlicedBalancedTraversal.h"
#include "autopas/containers/verletClusterLists/traversals/VCLSlicedC02Traversal.h"
#include "autopas/containers/verletClusterLists/traversals/VCLSlicedTraversal.h"
#include "autopas/containers/verletListsCellBased/varVerletLists/traversals/VVLAsBuildTraversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLListIterationTraversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/VLCAllCellsNeighborList.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/VLCCellPairNeighborList.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCC01Traversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCC18Traversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCCellPairC08Traversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCSlicedBalancedTraversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCSlicedC02Traversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCSlicedTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/SelectorStrategyOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/tuning/selectors/TraversalSelectorInfo.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/StringUtils.h"
#include "autopas/utils/TrivialHash.h"
#include "autopas/utils/logging/Logger.h"
#include "autopas/utils/checkFunctorType.h"
#include "autopas/pairwiseFunctors/TriwiseFunctor.h"

namespace autopas {

/**
 * Selector for a container traversal.
 * @tparam ParticleCell
 */
template <class ParticleCell, InteractionTypeOption::Value interactionType>
class TraversalSelector {
 public:
  /**
   * Generates a given Traversal for the given properties.
   * @tparam Functor
   * @tparam useSoA
   * @tparam useNewton3
   * @param traversalType
   * @param functor
   * @param info
   * @return Smartpointer to the traversal.
   */
  template <class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
  static std::unique_ptr<TraversalInterface<interactionType>> generateTraversal(TraversalOption traversalType,
                                                                       Functor &functor,
                                                               const TraversalSelectorInfo &info);

  /**
   * Generates a given Traversal for the given properties.
   * Requires less templates but calls the templated version after a decision tree.
   * @tparam Functor
   * @param traversalType
   * @param functor
   * @param info
   * @param dataLayout
   * @param useNewton3
   * @return Smartpointer to the traversal.
   */
  template <class Functor>
  static std::unique_ptr<TraversalInterface<interactionType>> generateTraversal(TraversalOption traversalType,
                                                                       Functor &functor,
                                                               const TraversalSelectorInfo &info,
                                                               DataLayoutOption dataLayout, Newton3Option useNewton3);

 private:
  template <class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
  static std::unique_ptr<TraversalInterface<InteractionTypeOption::pairwise>> generatePairwiseTraversal(TraversalOption traversalType,
                                                                       PairwiseFunctor &pairwiseFunctor,
                                                                       const TraversalSelectorInfo &info);

  template <class TriwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
  static std::unique_ptr<TraversalInterface<InteractionTypeOption::threeBody>> generateTriwiseTraversal(TraversalOption traversalType,
                                                                                                        TriwiseFunctor &triwiseFunctor,
                                                                               const TraversalSelectorInfo &info);
};

template <class ParticleCell, InteractionTypeOption::Value interactionType>
template <class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
std::unique_ptr<TraversalInterface<InteractionTypeOption::pairwise>> TraversalSelector<ParticleCell, interactionType>::generatePairwiseTraversal(
    TraversalOption traversalType, PairwiseFunctor &pairwiseFunctor, const TraversalSelectorInfo &info) {
  switch (traversalType) {
    // Direct sum
    case TraversalOption::ds_sequential: {
      return std::make_unique<DSSequentialTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          &pairwiseFunctor,
          info.interactionLength /*this is the cutoff, as generated by DirectSum::getTraversalSelectorInfo()!*/);
    }
    // Linked cell
    case TraversalOption::lc_sliced: {
      return std::make_unique<LCSlicedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::lc_sliced_c02: {
      return std::make_unique<LCSlicedC02Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::lc_sliced_balanced: {
      return std::make_unique<LCSlicedBalancedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::lc_c01: {
      return std::make_unique<LCC01Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::lc_c01_combined_SoA: {
      return std::make_unique<LCC01Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, true>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::lc_c04_combined_SoA: {
      return std::make_unique<LCC04CombinedSoATraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::lc_c04: {
      return std::make_unique<LCC04Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::lc_c04_HCP: {
      return std::make_unique<LCC04HCPTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::lc_c08: {
      return std::make_unique<LCC08Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::lc_c18: {
      return std::make_unique<LCC18Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    // Verlet
    case TraversalOption::vl_list_iteration: {
      return std::make_unique<VLListIterationTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          &pairwiseFunctor);
    }
    // Var Verlet Lists
    case TraversalOption::vvl_as_built: {
      return std::make_unique<VVLAsBuildTraversal<ParticleCell, typename ParticleCell::ParticleType, PairwiseFunctor,
                                                  dataLayout, useNewton3>>(&pairwiseFunctor);
    }
    // Verlet List Cells
    case TraversalOption::vlc_sliced: {
      return std::make_unique<VLCSlicedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3,
                                                 VLCAllCellsNeighborList<typename ParticleCell::ParticleType>,
                                                 ContainerOption::verletListsCells>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::vlc_sliced_c02: {
      return std::make_unique<VLCSlicedC02Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3,
                                                    VLCAllCellsNeighborList<typename ParticleCell::ParticleType>,
                                                    ContainerOption::verletListsCells>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::vlc_sliced_balanced: {
      return std::make_unique<VLCSlicedBalancedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3,
                                                         VLCAllCellsNeighborList<typename ParticleCell::ParticleType>,
                                                         ContainerOption::verletListsCells>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::vlc_c01: {
      return std::make_unique<VLCC01Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3,
                                              VLCAllCellsNeighborList<typename ParticleCell::ParticleType>,
                                              ContainerOption::verletListsCells>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::vlc_c18: {
      return std::make_unique<VLCC18Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3,
                                              VLCAllCellsNeighborList<typename ParticleCell::ParticleType>,
                                              ContainerOption::verletListsCells>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    // Verlet Cluster Lists
    case TraversalOption::vcl_cluster_iteration: {
      return std::make_unique<VCLClusterIterationTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          &pairwiseFunctor, info.clusterSize);
    }
    case TraversalOption::vcl_c01_balanced: {
      return std::make_unique<
          VCLC01BalancedTraversal<typename ParticleCell::ParticleType, PairwiseFunctor, dataLayout, useNewton3>>(
          &pairwiseFunctor, info.clusterSize);
    }
    case TraversalOption::vcl_sliced: {
      return std::make_unique<VCLSlicedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength, info.clusterSize);
    }
    case TraversalOption::vcl_sliced_c02: {
      return std::make_unique<VCLSlicedC02Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength, info.clusterSize);
    }
    case TraversalOption::vcl_sliced_balanced: {
      return std::make_unique<VCLSlicedBalancedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength, info.clusterSize);
    }
    case TraversalOption::vcl_c06: {
      return std::make_unique<VCLC06Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(&pairwiseFunctor,
                                                                                                      info.clusterSize);
    }
    // Pairwise Verlet Lists
    case TraversalOption::vlp_sliced: {
      return std::make_unique<VLCSlicedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3,
                                                 VLCCellPairNeighborList<typename ParticleCell::ParticleType>,
                                                 ContainerOption::pairwiseVerletLists>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::vlp_sliced_c02: {
      return std::make_unique<VLCSlicedC02Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3,
                                                    VLCCellPairNeighborList<typename ParticleCell::ParticleType>,
                                                    ContainerOption::pairwiseVerletLists>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::vlp_sliced_balanced: {
      return std::make_unique<VLCSlicedBalancedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3,
                                                         VLCCellPairNeighborList<typename ParticleCell::ParticleType>,
                                                         ContainerOption::pairwiseVerletLists>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::vlp_c01: {
      return std::make_unique<VLCC01Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3,
                                              VLCCellPairNeighborList<typename ParticleCell::ParticleType>,
                                              ContainerOption::pairwiseVerletLists>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::vlp_c18: {
      return std::make_unique<VLCC18Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3,
                                              VLCCellPairNeighborList<typename ParticleCell::ParticleType>,
                                              ContainerOption::pairwiseVerletLists>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::vlp_c08: {
      return std::make_unique<VLCCellPairC08Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    // Octree
    case TraversalOption::ot_c18: {
      using ParticleType = typename ParticleCell::ParticleType;
      return std::make_unique<OTC18Traversal<ParticleType, PairwiseFunctor, dataLayout, useNewton3>>(
          &pairwiseFunctor, info.interactionLength, info.interactionLength);
    }

    case TraversalOption::ot_c01: {
      using ParticleType = typename ParticleCell::ParticleType;
      return std::make_unique<OTC01Traversal<ParticleType, PairwiseFunctor, dataLayout, useNewton3>>(
          &pairwiseFunctor, info.interactionLength, info.interactionLength);
    }
    default: {
      autopas::utils::ExceptionHandler::exception("Traversal type {} is not a known pairwise traversal type!", traversalType.to_string());
      return {nullptr};
    }
  }
}

template <class ParticleCell, InteractionTypeOption::Value interactionType>
template <class TriwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
std::unique_ptr<TraversalInterface<InteractionTypeOption::threeBody>> TraversalSelector<ParticleCell, interactionType>::generateTriwiseTraversal(
    TraversalOption traversalType, TriwiseFunctor &triwiseFunctor, const TraversalSelectorInfo &info) {
  switch (traversalType) {
      // Direct sum
      case TraversalOption::ds_sequential_3b: {
        return std::make_unique<DSSequentialTraversal3B<ParticleCell, TriwiseFunctor, dataLayout, useNewton3>>(
            &triwiseFunctor,
            info.interactionLength /*this is the cutoff, as generated by DirectSum::getTraversalSelectorInfo()!*/);
      }
        // Linked Cells
      case TraversalOption::lc_c01_3b: {
        return std::make_unique<LCC01Traversal3B<ParticleCell, TriwiseFunctor, dataLayout, useNewton3>>(
            info.cellsPerDim, &triwiseFunctor, info.interactionLength, info.cellLength);
      }
      default: {
        autopas::utils::ExceptionHandler::exception("Traversal type {} is not a known 3-body traversal type!", traversalType.to_string());
        return {nullptr};
      }
  }
}

template <class ParticleCell, InteractionTypeOption::Value interactionType>
template <class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
std::unique_ptr<TraversalInterface<interactionType>> TraversalSelector<ParticleCell, interactionType>::generateTraversal(
    TraversalOption traversalType, Functor &functor, const TraversalSelectorInfo &info) {
  if constexpr (utils::isPairwiseFunctor<Functor>()) {
      return generatePairwiseTraversal<Functor, dataLayout, useNewton3>(traversalType, functor, info);
  }
  else if constexpr (utils::isTriwiseFunctor<Functor>()) {
      return generateTriwiseTraversal<Functor, dataLayout, useNewton3>(traversalType, functor, info);
  }
  return {nullptr};
}

template <class ParticleCell, InteractionTypeOption::Value interactionType>
template <class PairwiseFunctor>
std::unique_ptr<TraversalInterface<interactionType>> TraversalSelector<ParticleCell, interactionType>::generateTraversal(
    TraversalOption traversalType, PairwiseFunctor &pairwiseFunctor, const TraversalSelectorInfo &traversalInfo,
    DataLayoutOption dataLayout, Newton3Option newton3) {
  switch (dataLayout) {
    case DataLayoutOption::aos: {
      if (newton3 == Newton3Option::enabled) {
        return TraversalSelector<ParticleCell, interactionType>::template generateTraversal<PairwiseFunctor, DataLayoutOption::aos,
                                                                           true>(traversalType, pairwiseFunctor,
                                                                                 traversalInfo);
      } else {
        return TraversalSelector<ParticleCell, interactionType>::template generateTraversal<PairwiseFunctor, DataLayoutOption::aos,
                                                                           false>(traversalType, pairwiseFunctor,
                                                                                  traversalInfo);
      }
    }
    case DataLayoutOption::soa: {
      if (newton3 == Newton3Option::enabled) {
        return TraversalSelector<ParticleCell, interactionType>::template generateTraversal<PairwiseFunctor, DataLayoutOption::soa,
                                                                           true>(traversalType, pairwiseFunctor,
                                                                                 traversalInfo);
      } else {
        return TraversalSelector<ParticleCell, interactionType>::template generateTraversal<PairwiseFunctor, DataLayoutOption::soa,
                                                                           false>(traversalType, pairwiseFunctor,
                                                                                  traversalInfo);
      }
    }
  }

  autopas::utils::ExceptionHandler::exception("Traversal type {} is not a known type!", traversalType.to_string());
  return {nullptr};
}
}  // namespace autopas
