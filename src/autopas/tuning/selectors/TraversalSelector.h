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
#include "autopas/containers/linkedCells/traversals/LCC01Traversal.h"
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
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCC08Traversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCC18Traversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCCellPairC08Traversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCSlicedBalancedTraversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCSlicedC02Traversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCSlicedTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/SelectorStrategyOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/tuning/selectors/TraversalSelectorInfo.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/StringUtils.h"
#include "autopas/utils/TrivialHash.h"
#include "autopas/utils/checkFunctorType.h"
#include "autopas/utils/logging/Logger.h"

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
   * @tparam Functor
   * @param traversalType
   * @param functor
   * @param traversalInfo
   * @param dataLayout
   * @param useNewton3
   * @return Smartpointer to the traversal.
   */
  template <class Functor>
  static std::optional<std::unique_ptr<TraversalInterface>> generateTraversal(TraversalOption traversalType, Functor &functor,
                                                               const TraversalSelectorInfo &traversalInfo,
                                                               DataLayoutOption dataLayout, bool useNewton3);

  /**
   * Generates a given pairwise Traversal for the given properties.
   * @tparam PairwiseFunctor
   * @param traversalType
   * @param pairwiseFunctor
   * @param traversalInfo
   * @param dataLayout
   * @param useNewton3
   * @return Smartpointer to the traversal.
   */
  template <class PairwiseFunctor>
  static std::optional<std::unique_ptr<TraversalInterface>> generatePairwiseTraversal(TraversalOption traversalType,
                                                                       PairwiseFunctor &pairwiseFunctor,
                                                                       const TraversalSelectorInfo &traversalInfo,
                                                                       DataLayoutOption dataLayout, bool useNewton3);

  /**
   * Generates a given triwise Traversal for the given properties.
   * @tparam TriwiseFunctor
   * @param traversalType
   * @param triwiseFunctor
   * @param traversalInfo
   * @param dataLayout
   * @param useNewton3
   * @return Smartpointer to the traversal.
   */
  template <class TriwiseFunctor>
  static std::optional<std::unique_ptr<TraversalInterface>> generateTriwiseTraversal(TraversalOption traversalType,
                                                                      TriwiseFunctor &triwiseFunctor,
                                                                      const TraversalSelectorInfo &traversalInfo,
                                                                      DataLayoutOption dataLayout, bool useNewton3);
};

template <class ParticleCell>
template <class PairwiseFunctor>
std::optional<std::unique_ptr<TraversalInterface>> TraversalSelector<ParticleCell>::generatePairwiseTraversal(
    TraversalOption traversalType, PairwiseFunctor &pairwiseFunctor, const TraversalSelectorInfo &traversalInfo,
    DataLayoutOption dataLayout, bool useNewton3) {
  std::unique_ptr<TraversalInterface> traversal;
  switch (traversalType) {
    // Direct sum
    case TraversalOption::ds_sequential: {
      traversal = std::make_unique<DSSequentialTraversal<ParticleCell, PairwiseFunctor>>(
          &pairwiseFunctor,
          traversalInfo
              .interactionLength /*this is the cutoff, as generated by DirectSum::getTraversalSelectorInfo()!*/,
          dataLayout, useNewton3);
      break;
    }
    // Linked cell
    case TraversalOption::lc_sliced: {
      traversal = std::make_unique<LCSlicedTraversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_sliced_c02: {
      traversal = std::make_unique<LCSlicedC02Traversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_sliced_balanced: {
      traversal = std::make_unique<LCSlicedBalancedTraversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_c01: {
      traversal = std::make_unique<LCC01Traversal<ParticleCell, PairwiseFunctor, false>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_c01_combined_SoA: {
      traversal = std::make_unique<LCC01Traversal<ParticleCell, PairwiseFunctor, true>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_c04_combined_SoA: {
      traversal = std::make_unique<LCC04CombinedSoATraversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_c04: {
      traversal = std::make_unique<LCC04Traversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_c04_HCP: {
      traversal = std::make_unique<LCC04HCPTraversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_c08: {
      traversal = std::make_unique<LCC08Traversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_c18: {
      traversal = std::make_unique<LCC18Traversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    // Verlet
    case TraversalOption::vl_list_iteration: {
      traversal = std::make_unique<VLListIterationTraversal<ParticleCell, PairwiseFunctor>>(&pairwiseFunctor, dataLayout,
                                                                                       useNewton3);
      break;
    }
    // Var Verlet Lists
    case TraversalOption::vvl_as_built: {
      traversal = std::make_unique<VVLAsBuildTraversal<ParticleCell, typename ParticleCell::ParticleType, PairwiseFunctor>>(
          &pairwiseFunctor, dataLayout, useNewton3);
      break;
    }
    // Verlet List Cells
    case TraversalOption::vlc_sliced: {
      traversal = std::make_unique<VLCSlicedTraversal<ParticleCell, PairwiseFunctor,
                                                 VLCAllCellsNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::verletListsCells);
      break;
    }
    case TraversalOption::vlc_sliced_c02: {
      traversal = std::make_unique<VLCSlicedC02Traversal<ParticleCell, PairwiseFunctor,
                                                    VLCAllCellsNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::verletListsCells);
      break;
    }
    case TraversalOption::vlc_sliced_balanced: {
      traversal = std::make_unique<VLCSlicedBalancedTraversal<ParticleCell, PairwiseFunctor,
                                                         VLCAllCellsNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::verletListsCells);
      break;
    }
    case TraversalOption::vlc_c01: {
      traversal = std::make_unique<
          VLCC01Traversal<ParticleCell, PairwiseFunctor, VLCAllCellsNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::verletListsCells);
      break;
    }
    case TraversalOption::vlc_c18: {
      traversal = std::make_unique<
          VLCC18Traversal<ParticleCell, PairwiseFunctor, VLCAllCellsNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::verletListsCells);
      break;
    }
    case TraversalOption::vlc_c08: {
      traversal = std::make_unique<
          VLCC08Traversal<ParticleCell, PairwiseFunctor, VLCAllCellsNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    // Verlet Cluster Lists
    case TraversalOption::vcl_cluster_iteration: {
      traversal = std::make_unique<VCLClusterIterationTraversal<ParticleCell, PairwiseFunctor>>(
          &pairwiseFunctor, traversalInfo.clusterSize, dataLayout, useNewton3);
      break;
    }
    case TraversalOption::vcl_c01_balanced: {
      traversal = std::make_unique<VCLC01BalancedTraversal<typename ParticleCell::ParticleType, PairwiseFunctor>>(
          &pairwiseFunctor, traversalInfo.clusterSize, dataLayout, useNewton3);
      break;
    }
    case TraversalOption::vcl_sliced: {
      traversal = std::make_unique<VCLSlicedTraversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          traversalInfo.clusterSize, dataLayout, useNewton3);
      break;
    }
    case TraversalOption::vcl_sliced_c02: {
      traversal = std::make_unique<VCLSlicedC02Traversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          traversalInfo.clusterSize, dataLayout, useNewton3);
      break;
    }
    case TraversalOption::vcl_sliced_balanced: {
      traversal = std::make_unique<VCLSlicedBalancedTraversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          traversalInfo.clusterSize, dataLayout, useNewton3);
      break;
    }
    case TraversalOption::vcl_c06: {
      traversal = std::make_unique<VCLC06Traversal<ParticleCell, PairwiseFunctor>>(
          &pairwiseFunctor, traversalInfo.clusterSize, dataLayout, useNewton3);
      break;
    }
    // Pairwise Verlet Lists
    case TraversalOption::vlp_sliced: {
      traversal = std::make_unique<VLCSlicedTraversal<ParticleCell, PairwiseFunctor,
                                                 VLCCellPairNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::pairwiseVerletLists);
      break;
    }
    case TraversalOption::vlp_sliced_c02: {
      traversal = std::make_unique<VLCSlicedC02Traversal<ParticleCell, PairwiseFunctor,
                                                    VLCCellPairNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::pairwiseVerletLists);
      break;
    }
    case TraversalOption::vlp_sliced_balanced: {
      traversal = std::make_unique<VLCSlicedBalancedTraversal<ParticleCell, PairwiseFunctor,
                                                         VLCCellPairNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::pairwiseVerletLists);
      break;
    }
    case TraversalOption::vlp_c01: {
      traversal = std::make_unique<
          VLCC01Traversal<ParticleCell, PairwiseFunctor, VLCCellPairNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::pairwiseVerletLists);
      break;
    }
    case TraversalOption::vlp_c18: {
      traversal = std::make_unique<
          VLCC18Traversal<ParticleCell, PairwiseFunctor, VLCCellPairNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::pairwiseVerletLists);
      break;
    }
    case TraversalOption::vlp_c08: {
      traversal = std::make_unique<VLCCellPairC08Traversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    // Octree
    case TraversalOption::ot_c18: {
      using ParticleType = typename ParticleCell::ParticleType;
      traversal = std::make_unique<OTC18Traversal<ParticleType, PairwiseFunctor>>(
          &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.interactionLength, dataLayout, useNewton3);
      break;
    }

    case TraversalOption::ot_c01: {
      using ParticleType = typename ParticleCell::ParticleType;
      traversal = std::make_unique<OTC01Traversal<ParticleType, PairwiseFunctor>>(
          &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.interactionLength, dataLayout, useNewton3);
      break;
    }
    default: {
      utils::ExceptionHandler::exception("Traversal type {} is not a known pairwise traversal type!",
                                                  traversalType.to_string());
      return std::nullopt;
    }
  }
  // Check if the traversal is applicable.
  if (not traversal->isApplicable()) {
    return std::nullopt;
  }
  // If applicable, return the traversal.
  return std::move(traversal);
}

template <class ParticleCell>
template <class TriwiseFunctor>
std::optional<std::unique_ptr<TraversalInterface>> TraversalSelector<ParticleCell>::generateTriwiseTraversal(
    TraversalOption traversalType, TriwiseFunctor &triwiseFunctor, const TraversalSelectorInfo &traversalInfo,
    DataLayoutOption dataLayout, bool useNewton3) {
  std::unique_ptr<TraversalInterface> traversal;
  switch (traversalType) {
    // Direct sum
    case TraversalOption::ds_sequential: {
      traversal = std::make_unique<DSSequentialTraversal<ParticleCell, TriwiseFunctor>>(
          &triwiseFunctor,
          traversalInfo
              .interactionLength /*this is the cutoff, as generated by DirectSum::getTraversalSelectorInfo()!*/,
          dataLayout, useNewton3);
      break;
    }
      // Linked Cells
    case TraversalOption::lc_c01: {
      traversal = std::make_unique<LCC01Traversal<ParticleCell, TriwiseFunctor, /*combineSoA*/ false>>(
          traversalInfo.cellsPerDim, &triwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    default: {
      utils::ExceptionHandler::exception("Traversal type {} is not a known triwise traversal type!",
                                                  traversalType.to_string());
      return std::nullopt;
    }
  }
  // Check if the traversal is applicable.
  if (not traversal->isApplicable()) {
    return std::nullopt;
  }
  // If applicable, return the traversal.
  return std::move(traversal);
}

template <class ParticleCell>
template <class Functor>
std::optional<std::unique_ptr<TraversalInterface>> TraversalSelector<ParticleCell>::generateTraversal(
    TraversalOption traversalType, Functor &functor, const TraversalSelectorInfo &traversalInfo,
    DataLayoutOption dataLayout, bool useNewton3) {
  if constexpr (utils::isPairwiseFunctor<Functor>()) {
    return generatePairwiseTraversal<Functor>(traversalType, functor, traversalInfo, dataLayout, useNewton3);
  } else if constexpr (utils::isTriwiseFunctor<Functor>()) {
    return generateTriwiseTraversal<Functor>(traversalType, functor, traversalInfo, dataLayout, useNewton3);
  }
  utils::ExceptionHandler::exception(
      "TraversalSelector::generateTraversal(): No Traversals exist for the given Functor: {}", functor.getName());
  return std::nullopt;
}
}  // namespace autopas
