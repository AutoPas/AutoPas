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
#include "autopas/containers/linkedCells/traversals/HGC01Traversal.h"
#include "autopas/containers/linkedCells/traversals/HGColorSoACellToCell.h"
#include "autopas/containers/linkedCells/traversals/HGColorTraversal.h"
#include "autopas/containers/linkedCells/traversals/HGColorTraversalC04.h"
#include "autopas/containers/linkedCells/traversals/HGColorTraversalC04Combined.h"
#include "autopas/containers/linkedCells/traversals/HGColorTraversalC04HCP.h"
#include "autopas/containers/linkedCells/traversals/HGColorTraversalC18.h"
#include "autopas/containers/linkedCells/traversals/HGTestTraversal.h"
#include "autopas/containers/linkedCells/traversals/HGTestTraversal2.h"
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
  static std::unique_ptr<TraversalInterface> generateTraversal(TraversalOption traversalType, Functor &functor,
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
  static std::unique_ptr<TraversalInterface> generatePairwiseTraversal(TraversalOption traversalType,
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
  static std::unique_ptr<TraversalInterface> generateTriwiseTraversal(TraversalOption traversalType,
                                                                      TriwiseFunctor &triwiseFunctor,
                                                                      const TraversalSelectorInfo &traversalInfo,
                                                                      DataLayoutOption dataLayout, bool useNewton3);
};

template <class ParticleCell>
template <class PairwiseFunctor>
std::unique_ptr<TraversalInterface> TraversalSelector<ParticleCell>::generatePairwiseTraversal(
    TraversalOption traversalType, PairwiseFunctor &pairwiseFunctor, const TraversalSelectorInfo &traversalInfo,
    DataLayoutOption dataLayout, bool useNewton3) {
  switch (traversalType) {
    // Direct sum
    case TraversalOption::ds_sequential: {
      return std::make_unique<DSSequentialTraversal<ParticleCell, PairwiseFunctor>>(
          &pairwiseFunctor,
          traversalInfo
              .interactionLength /*this is the cutoff, as generated by DirectSum::getTraversalSelectorInfo()!*/,
          dataLayout, useNewton3);
    }
    // Linked cell
    case TraversalOption::lc_sliced: {
      return std::make_unique<LCSlicedTraversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
    }
    case TraversalOption::lc_sliced_c02: {
      return std::make_unique<LCSlicedC02Traversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
    }
    case TraversalOption::lc_sliced_balanced: {
      return std::make_unique<LCSlicedBalancedTraversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
    }
    case TraversalOption::lc_c01: {
      return std::make_unique<LCC01Traversal<ParticleCell, PairwiseFunctor, false>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
    }
    case TraversalOption::lc_c01_combined_SoA: {
      return std::make_unique<LCC01Traversal<ParticleCell, PairwiseFunctor, true>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
    }
    case TraversalOption::lc_c04_combined_SoA: {
      return std::make_unique<LCC04CombinedSoATraversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
    }
    case TraversalOption::lc_c04: {
      return std::make_unique<LCC04Traversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
    }
    case TraversalOption::lc_c04_HCP: {
      return std::make_unique<LCC04HCPTraversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
    }
    case TraversalOption::lc_c08: {
      return std::make_unique<LCC08Traversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
    }
    case TraversalOption::lc_c18: {
      return std::make_unique<LCC18Traversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
    }
    // Verlet
    case TraversalOption::vl_list_iteration: {
      return std::make_unique<VLListIterationTraversal<ParticleCell, PairwiseFunctor>>(&pairwiseFunctor, dataLayout,
                                                                                       useNewton3);
    }
    // Var Verlet Lists
    case TraversalOption::vvl_as_built: {
      return std::make_unique<VVLAsBuildTraversal<ParticleCell, typename ParticleCell::ParticleType, PairwiseFunctor>>(
          &pairwiseFunctor, dataLayout, useNewton3);
    }
    // Verlet List Cells
    case TraversalOption::vlc_sliced: {
      return std::make_unique<VLCSlicedTraversal<ParticleCell, PairwiseFunctor,
                                                 VLCAllCellsNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::verletListsCells);
    }
    case TraversalOption::vlc_sliced_c02: {
      return std::make_unique<VLCSlicedC02Traversal<ParticleCell, PairwiseFunctor,
                                                    VLCAllCellsNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::verletListsCells);
    }
    case TraversalOption::vlc_sliced_balanced: {
      return std::make_unique<VLCSlicedBalancedTraversal<ParticleCell, PairwiseFunctor,
                                                         VLCAllCellsNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::verletListsCells);
    }
    case TraversalOption::vlc_c01: {
      return std::make_unique<
          VLCC01Traversal<ParticleCell, PairwiseFunctor, VLCAllCellsNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::verletListsCells);
    }
    case TraversalOption::vlc_c18: {
      return std::make_unique<
          VLCC18Traversal<ParticleCell, PairwiseFunctor, VLCAllCellsNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::verletListsCells);
    }
    case TraversalOption::vlc_c08: {
      return std::make_unique<
          VLCC08Traversal<ParticleCell, PairwiseFunctor, VLCAllCellsNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
    }
    // Verlet Cluster Lists
    case TraversalOption::vcl_cluster_iteration: {
      return std::make_unique<VCLClusterIterationTraversal<ParticleCell, PairwiseFunctor>>(
          &pairwiseFunctor, traversalInfo.clusterSize, dataLayout, useNewton3);
    }
    case TraversalOption::vcl_c01_balanced: {
      return std::make_unique<VCLC01BalancedTraversal<typename ParticleCell::ParticleType, PairwiseFunctor>>(
          &pairwiseFunctor, traversalInfo.clusterSize, dataLayout, useNewton3);
    }
    case TraversalOption::vcl_sliced: {
      return std::make_unique<VCLSlicedTraversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          traversalInfo.clusterSize, dataLayout, useNewton3);
    }
    case TraversalOption::vcl_sliced_c02: {
      return std::make_unique<VCLSlicedC02Traversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          traversalInfo.clusterSize, dataLayout, useNewton3);
    }
    case TraversalOption::vcl_sliced_balanced: {
      return std::make_unique<VCLSlicedBalancedTraversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          traversalInfo.clusterSize, dataLayout, useNewton3);
    }
    case TraversalOption::vcl_c06: {
      return std::make_unique<VCLC06Traversal<ParticleCell, PairwiseFunctor>>(
          &pairwiseFunctor, traversalInfo.clusterSize, dataLayout, useNewton3);
    }
    // Pairwise Verlet Lists
    case TraversalOption::vlp_sliced: {
      return std::make_unique<VLCSlicedTraversal<ParticleCell, PairwiseFunctor,
                                                 VLCCellPairNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::pairwiseVerletLists);
    }
    case TraversalOption::vlp_sliced_c02: {
      return std::make_unique<VLCSlicedC02Traversal<ParticleCell, PairwiseFunctor,
                                                    VLCCellPairNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::pairwiseVerletLists);
    }
    case TraversalOption::vlp_sliced_balanced: {
      return std::make_unique<VLCSlicedBalancedTraversal<ParticleCell, PairwiseFunctor,
                                                         VLCCellPairNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::pairwiseVerletLists);
    }
    case TraversalOption::vlp_c01: {
      return std::make_unique<
          VLCC01Traversal<ParticleCell, PairwiseFunctor, VLCCellPairNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::pairwiseVerletLists);
    }
    case TraversalOption::vlp_c18: {
      return std::make_unique<
          VLCC18Traversal<ParticleCell, PairwiseFunctor, VLCCellPairNeighborList<typename ParticleCell::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::pairwiseVerletLists);
    }
    case TraversalOption::vlp_c08: {
      return std::make_unique<VLCCellPairC08Traversal<ParticleCell, PairwiseFunctor>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
    }
    // Octree
    case TraversalOption::ot_c18: {
      using ParticleType = typename ParticleCell::ParticleType;
      return std::make_unique<OTC18Traversal<ParticleType, PairwiseFunctor>>(
          &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.interactionLength, dataLayout, useNewton3);
    }

    case TraversalOption::ot_c01: {
      using ParticleType = typename ParticleCell::ParticleType;
      return std::make_unique<OTC01Traversal<ParticleType, PairwiseFunctor>>(
          &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.interactionLength, dataLayout, useNewton3);
    }
    // Hierarchical Grid
    case TraversalOption::hgrid_c01_iterator: {
      return std::make_unique<HGC01Traversal<ParticleCell, PairwiseFunctor>>(&pairwiseFunctor, dataLayout, useNewton3);
    }
    case TraversalOption::hgrid_color: {
      return std::make_unique<HGColorTraversal<ParticleCell, PairwiseFunctor>>(&pairwiseFunctor, dataLayout,
                                                                               useNewton3);
    }
    case TraversalOption::hgrid_color_soa_cell: {
      return std::make_unique<HGColorSoACellToCell<ParticleCell, PairwiseFunctor>>(&pairwiseFunctor, dataLayout,
                                                                                   useNewton3);
    }
    case TraversalOption::hgrid_test: {
      return std::make_unique<HGTestTraversal<ParticleCell, PairwiseFunctor>>(&pairwiseFunctor, dataLayout, useNewton3);
    }
    case TraversalOption::hgrid_test2: {
      return std::make_unique<HGTestTraversal2<ParticleCell, PairwiseFunctor>>(&pairwiseFunctor, dataLayout, useNewton3);
    }
    case TraversalOption::hgrid_color_c04: {
      return std::make_unique<HGColorTraversalC04<ParticleCell, PairwiseFunctor>>(&pairwiseFunctor, dataLayout,
                                                                                  useNewton3);
    }
    case TraversalOption::hgrid_color_c04_combined: {
      return std::make_unique<HGColorTraversalC04Combined<ParticleCell, PairwiseFunctor>>(&pairwiseFunctor, dataLayout,
                                                                                          useNewton3);
    }
    case TraversalOption::hgrid_color_c04_HCP: {
      return std::make_unique<HGColorTraversalC04HCP<ParticleCell, PairwiseFunctor>>(&pairwiseFunctor, dataLayout,
                                                                                     useNewton3);
    }
    case TraversalOption::hgrid_color_c18: {
      return std::make_unique<HGColorTraversalC18<ParticleCell, PairwiseFunctor>>(&pairwiseFunctor, dataLayout,
                                                                                  useNewton3);
    }
    default: {
      autopas::utils::ExceptionHandler::exception("Traversal type {} is not a known pairwise traversal type!",
                                                  traversalType.to_string());
      return {nullptr};
    }
  }
}

template <class ParticleCell>
template <class TriwiseFunctor>
std::unique_ptr<TraversalInterface> TraversalSelector<ParticleCell>::generateTriwiseTraversal(
    TraversalOption traversalType, TriwiseFunctor &triwiseFunctor, const TraversalSelectorInfo &traversalInfo,
    DataLayoutOption dataLayout, bool useNewton3) {
  switch (traversalType) {
    // Direct sum
    case TraversalOption::ds_sequential: {
      return std::make_unique<DSSequentialTraversal<ParticleCell, TriwiseFunctor>>(
          &triwiseFunctor,
          traversalInfo
              .interactionLength /*this is the cutoff, as generated by DirectSum::getTraversalSelectorInfo()!*/,
          dataLayout, useNewton3);
    }
      // Linked Cells
    case TraversalOption::lc_c01: {
      return std::make_unique<LCC01Traversal<ParticleCell, TriwiseFunctor, /*combineSoA*/ false>>(
          traversalInfo.cellsPerDim, &triwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
    }
    default: {
      autopas::utils::ExceptionHandler::exception("Traversal type {} is not a known triwise traversal type!",
                                                  traversalType.to_string());
      return {nullptr};
    }
  }
}

template <class ParticleCell>
template <class Functor>
std::unique_ptr<TraversalInterface> TraversalSelector<ParticleCell>::generateTraversal(
    TraversalOption traversalType, Functor &functor, const TraversalSelectorInfo &traversalInfo,
    DataLayoutOption dataLayout, bool useNewton3) {
  if constexpr (utils::isPairwiseFunctor<Functor>()) {
    auto pairTraversal =
        generatePairwiseTraversal<Functor>(traversalType, functor, traversalInfo, dataLayout, useNewton3);
    return std::unique_ptr<TraversalInterface>(dynamic_cast<TraversalInterface *>(pairTraversal.release()));
  } else if constexpr (utils::isTriwiseFunctor<Functor>()) {
    auto triTraversal =
        generateTriwiseTraversal<Functor>(traversalType, functor, traversalInfo, dataLayout, useNewton3);
    return std::unique_ptr<TraversalInterface>(dynamic_cast<TraversalInterface *>(triTraversal.release()));
  }
  autopas::utils::ExceptionHandler::exception(
      "TraversalSelector::generateTraversal(): No Traversals exist for the given Functor: {}", functor.getName());
  return {nullptr};
}
}  // namespace autopas
