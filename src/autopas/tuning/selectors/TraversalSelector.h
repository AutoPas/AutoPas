/**
 * @file TraversalSelector.h
 * @author F. Gratl
 * @date 11.06.18
 */

#pragma once

#include <memory>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/cells/ReferenceParticleCell.h"
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
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLListIntersectionTraversalHashing3B.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLListIntersectionTraversalSorted3B.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLListIterationTraversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLPairListIterationTraversal3B.h"
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
#include "autopas/options/TraversalOption.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/selectors/TraversalSelectorInfo.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/checkFunctorType.h"

namespace autopas {

/**
 * Selector for a container traversal.
 */
class TraversalSelector {
 public:
  /**
   * Generates a given Traversal for the given properties.
   * @tparam ParticleCell_T
   * @tparam Functor_T
   * @param traversalType
   * @param functor
   * @param traversalInfo
   * @param dataLayout
   * @param useNewton3
   * @return Smartpointer to the traversal.
   */
  template <class ParticleCell_T, class Functor_T>
  static std::unique_ptr<TraversalInterface> generateTraversal(TraversalOption traversalType, Functor_T &functor,
                                                               const TraversalSelectorInfo &traversalInfo,
                                                               DataLayoutOption dataLayout, bool useNewton3);

  /**
   * Generates a given pairwise Traversal for the given properties.
   * @tparam ParticleCell_T
   * @tparam PairwiseFunctor_T
   * @param traversalType
   * @param pairwiseFunctor
   * @param traversalInfo
   * @param dataLayout
   * @param useNewton3
   * @return Smartpointer to the traversal.
   */
  template <class ParticleCell_T, class PairwiseFunctor_T>
  static std::unique_ptr<TraversalInterface> generatePairwiseTraversal(TraversalOption traversalType,
                                                                       PairwiseFunctor_T &pairwiseFunctor,
                                                                       const TraversalSelectorInfo &traversalInfo,
                                                                       DataLayoutOption dataLayout, bool useNewton3);

  /**
   * Generates a given triwise Traversal for the given properties.
   * @tparam ParticleCell_T
   * @tparam TriwiseFunctor_T
   * @param traversalType
   * @param triwiseFunctor
   * @param traversalInfo
   * @param dataLayout
   * @param useNewton3
   * @return Smartpointer to the traversal.
   */
  template <class ParticleCell_T, class TriwiseFunctor_T>
  static std::unique_ptr<TraversalInterface> generateTriwiseTraversal(TraversalOption traversalType,
                                                                      TriwiseFunctor_T &triwiseFunctor,
                                                                      const TraversalSelectorInfo &traversalInfo,
                                                                      DataLayoutOption dataLayout, bool useNewton3);

  /**
   * Generates a traversal from the given configuration.
   * @tparam Particle_T
   * @tparam Functor_T
   * @param config The configuration to generate the traversal from.
   * @param functor The functor to use in the traversal.
   * @param traversalInfo Additional information for the traversal.
   * @return Smartpointer to the generated traversal, or nullptr if no valid traversal could be generated.
   */
  template <class Particle_T, class Functor_T>
  static std::unique_ptr<TraversalInterface> generateTraversalFromConfig(const Configuration &config,
                                                                         Functor_T &functor,
                                                                         const TraversalSelectorInfo &traversalInfo);
};

template <class ParticleCell_T, class PairwiseFunctor_T>
std::unique_ptr<TraversalInterface> TraversalSelector::generatePairwiseTraversal(
    TraversalOption traversalType, PairwiseFunctor_T &pairwiseFunctor, const TraversalSelectorInfo &traversalInfo,
    DataLayoutOption dataLayout, bool useNewton3) {
  std::unique_ptr<TraversalInterface> traversal;
  switch (traversalType) {
    // Direct sum
    case TraversalOption::ds_sequential: {
      traversal = std::make_unique<DSSequentialTraversal<ParticleCell_T, PairwiseFunctor_T>>(
          &pairwiseFunctor,
          traversalInfo
              .interactionLength /*this is the cutoff, as generated by DirectSum::getTraversalSelectorInfo()!*/,
          dataLayout, useNewton3);
      break;
    }
    // Linked cell
    case TraversalOption::lc_sliced: {
      traversal = std::make_unique<LCSlicedTraversal<ParticleCell_T, PairwiseFunctor_T>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_sliced_c02: {
      traversal = std::make_unique<LCSlicedC02Traversal<ParticleCell_T, PairwiseFunctor_T>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_sliced_balanced: {
      traversal = std::make_unique<LCSlicedBalancedTraversal<ParticleCell_T, PairwiseFunctor_T>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_c01: {
      traversal = std::make_unique<LCC01Traversal<ParticleCell_T, PairwiseFunctor_T, false>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_c01_combined_SoA: {
      traversal = std::make_unique<LCC01Traversal<ParticleCell_T, PairwiseFunctor_T, true>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_c04_combined_SoA: {
      traversal = std::make_unique<LCC04CombinedSoATraversal<ParticleCell_T, PairwiseFunctor_T>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_c04: {
      traversal = std::make_unique<LCC04Traversal<ParticleCell_T, PairwiseFunctor_T>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_c04_HCP: {
      traversal = std::make_unique<LCC04HCPTraversal<ParticleCell_T, PairwiseFunctor_T>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_c08: {
      traversal = std::make_unique<LCC08Traversal<ParticleCell_T, PairwiseFunctor_T>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_c18: {
      traversal = std::make_unique<LCC18Traversal<ParticleCell_T, PairwiseFunctor_T>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    // Verlet
    case TraversalOption::vl_list_iteration: {
      traversal = std::make_unique<VLListIterationTraversal<ParticleCell_T, PairwiseFunctor_T>>(&pairwiseFunctor,
                                                                                                dataLayout, useNewton3);
      break;
    }
    // Var Verlet Lists
    case TraversalOption::vvl_as_built: {
      traversal = std::make_unique<
          VVLAsBuildTraversal<ParticleCell_T, typename ParticleCell_T::ParticleType, PairwiseFunctor_T>>(
          &pairwiseFunctor, dataLayout, useNewton3);
      break;
    }
    // Verlet List Cells
    case TraversalOption::vlc_sliced: {
      traversal = std::make_unique<VLCSlicedTraversal<ParticleCell_T, PairwiseFunctor_T,
                                                      VLCAllCellsNeighborList<typename ParticleCell_T::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::verletListsCells);
      break;
    }
    case TraversalOption::vlc_sliced_c02: {
      traversal =
          std::make_unique<VLCSlicedC02Traversal<ParticleCell_T, PairwiseFunctor_T,
                                                 VLCAllCellsNeighborList<typename ParticleCell_T::ParticleType>>>(
              traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
              dataLayout, useNewton3, ContainerOption::verletListsCells);
      break;
    }
    case TraversalOption::vlc_sliced_balanced: {
      traversal =
          std::make_unique<VLCSlicedBalancedTraversal<ParticleCell_T, PairwiseFunctor_T,
                                                      VLCAllCellsNeighborList<typename ParticleCell_T::ParticleType>>>(
              traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
              dataLayout, useNewton3, ContainerOption::verletListsCells);
      break;
    }
    case TraversalOption::vlc_c01: {
      traversal = std::make_unique<VLCC01Traversal<ParticleCell_T, PairwiseFunctor_T,
                                                   VLCAllCellsNeighborList<typename ParticleCell_T::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::verletListsCells);
      break;
    }
    case TraversalOption::vlc_c18: {
      traversal = std::make_unique<VLCC18Traversal<ParticleCell_T, PairwiseFunctor_T,
                                                   VLCAllCellsNeighborList<typename ParticleCell_T::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::verletListsCells);
      break;
    }
    case TraversalOption::vlc_c08: {
      traversal = std::make_unique<VLCC08Traversal<ParticleCell_T, PairwiseFunctor_T,
                                                   VLCAllCellsNeighborList<typename ParticleCell_T::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    // Verlet Cluster Lists
    case TraversalOption::vcl_cluster_iteration: {
      traversal = std::make_unique<VCLClusterIterationTraversal<ParticleCell_T, PairwiseFunctor_T>>(
          &pairwiseFunctor, traversalInfo.clusterSize, dataLayout, useNewton3);
      break;
    }
    case TraversalOption::vcl_c01_balanced: {
      traversal = std::make_unique<VCLC01BalancedTraversal<typename ParticleCell_T::ParticleType, PairwiseFunctor_T>>(
          &pairwiseFunctor, traversalInfo.clusterSize, dataLayout, useNewton3);
      break;
    }
    case TraversalOption::vcl_sliced: {
      traversal = std::make_unique<VCLSlicedTraversal<ParticleCell_T, PairwiseFunctor_T>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          traversalInfo.clusterSize, dataLayout, useNewton3);
      break;
    }
    case TraversalOption::vcl_sliced_c02: {
      traversal = std::make_unique<VCLSlicedC02Traversal<ParticleCell_T, PairwiseFunctor_T>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          traversalInfo.clusterSize, dataLayout, useNewton3);
      break;
    }
    case TraversalOption::vcl_sliced_balanced: {
      traversal = std::make_unique<VCLSlicedBalancedTraversal<ParticleCell_T, PairwiseFunctor_T>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          traversalInfo.clusterSize, dataLayout, useNewton3);
      break;
    }
    case TraversalOption::vcl_c06: {
      traversal = std::make_unique<VCLC06Traversal<ParticleCell_T, PairwiseFunctor_T>>(
          &pairwiseFunctor, traversalInfo.clusterSize, dataLayout, useNewton3);
      break;
    }
    // Pairwise Verlet Lists
    case TraversalOption::vlp_sliced: {
      traversal = std::make_unique<VLCSlicedTraversal<ParticleCell_T, PairwiseFunctor_T,
                                                      VLCCellPairNeighborList<typename ParticleCell_T::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::pairwiseVerletLists);
      break;
    }
    case TraversalOption::vlp_sliced_c02: {
      traversal =
          std::make_unique<VLCSlicedC02Traversal<ParticleCell_T, PairwiseFunctor_T,
                                                 VLCCellPairNeighborList<typename ParticleCell_T::ParticleType>>>(
              traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
              dataLayout, useNewton3, ContainerOption::pairwiseVerletLists);
      break;
    }
    case TraversalOption::vlp_sliced_balanced: {
      traversal =
          std::make_unique<VLCSlicedBalancedTraversal<ParticleCell_T, PairwiseFunctor_T,
                                                      VLCCellPairNeighborList<typename ParticleCell_T::ParticleType>>>(
              traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
              dataLayout, useNewton3, ContainerOption::pairwiseVerletLists);
      break;
    }
    case TraversalOption::vlp_c01: {
      traversal = std::make_unique<VLCC01Traversal<ParticleCell_T, PairwiseFunctor_T,
                                                   VLCCellPairNeighborList<typename ParticleCell_T::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::pairwiseVerletLists);
      break;
    }
    case TraversalOption::vlp_c18: {
      traversal = std::make_unique<VLCC18Traversal<ParticleCell_T, PairwiseFunctor_T,
                                                   VLCCellPairNeighborList<typename ParticleCell_T::ParticleType>>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3, ContainerOption::pairwiseVerletLists);
      break;
    }
    case TraversalOption::vlp_c08: {
      traversal = std::make_unique<VLCCellPairC08Traversal<ParticleCell_T, PairwiseFunctor_T>>(
          traversalInfo.cellsPerDim, &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    // Octree
    case TraversalOption::ot_c18: {
      using ParticleType = typename ParticleCell_T::ParticleType;
      traversal = std::make_unique<OTC18Traversal<ParticleType, PairwiseFunctor_T>>(
          &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.interactionLength, dataLayout, useNewton3);
      break;
    }

    case TraversalOption::ot_c01: {
      using ParticleType = typename ParticleCell_T::ParticleType;
      traversal = std::make_unique<OTC01Traversal<ParticleType, PairwiseFunctor_T>>(
          &pairwiseFunctor, traversalInfo.interactionLength, traversalInfo.interactionLength, dataLayout, useNewton3);
      break;
    }
    default: {
      utils::ExceptionHandler::exception("Traversal type {} is not a known pairwise traversal type!",
                                         traversalType.to_string());
      return nullptr;
    }
  }
  // Check if the traversal is applicable.
  if (not traversal->isApplicable()) {
    return nullptr;
  }
  // If applicable, return the traversal.
  return std::move(traversal);
}

template <class ParticleCell_T, class TriwiseFunctor_T>
std::unique_ptr<TraversalInterface> TraversalSelector::generateTriwiseTraversal(
    TraversalOption traversalType, TriwiseFunctor_T &triwiseFunctor, const TraversalSelectorInfo &traversalInfo,
    DataLayoutOption dataLayout, bool useNewton3) {
  std::unique_ptr<TraversalInterface> traversal;
  switch (traversalType) {
    // Direct sum
    case TraversalOption::ds_sequential: {
      traversal = std::make_unique<DSSequentialTraversal<ParticleCell_T, TriwiseFunctor_T>>(
          &triwiseFunctor,
          traversalInfo
              .interactionLength /*this is the cutoff, as generated by DirectSum::getTraversalSelectorInfo()!*/,
          dataLayout, useNewton3);
      break;
    }
      // Linked Cells
    case TraversalOption::lc_c01: {
      traversal = std::make_unique<LCC01Traversal<ParticleCell_T, TriwiseFunctor_T, /*combineSoA*/ false>>(
          traversalInfo.cellsPerDim, &triwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_c04: {
      traversal = std::make_unique<LCC04Traversal<ParticleCell_T, TriwiseFunctor_T>>(
          traversalInfo.cellsPerDim, &triwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_c08: {
      traversal = std::make_unique<LCC08Traversal<ParticleCell_T, TriwiseFunctor_T>>(
          traversalInfo.cellsPerDim, &triwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_sliced: {
      traversal = std::make_unique<LCSlicedTraversal<ParticleCell_T, TriwiseFunctor_T>>(
          traversalInfo.cellsPerDim, &triwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    case TraversalOption::lc_sliced_c02: {
      traversal = std::make_unique<LCSlicedC02Traversal<ParticleCell_T, TriwiseFunctor_T>>(
          traversalInfo.cellsPerDim, &triwiseFunctor, traversalInfo.interactionLength, traversalInfo.cellLength,
          dataLayout, useNewton3);
      break;
    }
    // VerletLists
    case TraversalOption::vl_list_iteration: {
      traversal = std::make_unique<VLListIterationTraversal<ParticleCell_T, TriwiseFunctor_T>>(&triwiseFunctor, dataLayout,
                                                                                          useNewton3);
      break;
    }
    case TraversalOption::vl_list_intersection_sorted_3b: {
      traversal = std::make_unique<VLListIntersectionTraversalSorted3B<ParticleCell_T, TriwiseFunctor_T>>(
          &triwiseFunctor, dataLayout, useNewton3);
      break;
    }
    case TraversalOption::vl_list_intersection_hashing_3b: {
      traversal = std::make_unique<VLListIntersectionTraversalHashing3B<ParticleCell_T, TriwiseFunctor_T>>(
          &triwiseFunctor, dataLayout, useNewton3);
      break;
    }
    case TraversalOption::vl_pair_list_iteration_3b: {
      traversal = std::make_unique<VLPairListIterationTraversal3B<ParticleCell_T, TriwiseFunctor_T>>(&triwiseFunctor,
                                                                                                dataLayout, useNewton3);
      break;
    }
    default: {
      utils::ExceptionHandler::exception("Traversal type {} is not a known triwise traversal type!",
                                         traversalType.to_string());
      return nullptr;
    }
  }
  // Check if the traversal is applicable.
  if (not traversal->isApplicable()) {
    return nullptr;
  }
  // If applicable, return the traversal.
  return std::move(traversal);
}

template <class ParticleCell_T, class Functor_T>
std::unique_ptr<TraversalInterface> TraversalSelector::generateTraversal(TraversalOption traversalType,
                                                                         Functor_T &functor,
                                                                         const TraversalSelectorInfo &traversalInfo,
                                                                         DataLayoutOption dataLayout, bool useNewton3) {
  if constexpr (utils::isPairwiseFunctor<Functor_T>()) {
    return generatePairwiseTraversal<ParticleCell_T, Functor_T>(traversalType, functor, traversalInfo, dataLayout,
                                                                useNewton3);
  } else if constexpr (utils::isTriwiseFunctor<Functor_T>()) {
    return generateTriwiseTraversal<ParticleCell_T, Functor_T>(traversalType, functor, traversalInfo, dataLayout,
                                                               useNewton3);
  }
  utils::ExceptionHandler::exception(
      "TraversalSelector::generateTraversal(): No Traversals exist for the given Functor: {}", functor.getName());
  return nullptr;
}

template <class Particle_T, class Functor_T>
std::unique_ptr<TraversalInterface> TraversalSelector::generateTraversalFromConfig(
    const Configuration &config, Functor_T &functor, const TraversalSelectorInfo &traversalInfo) {
  switch (config.container) {
    case ContainerOption::Value::linkedCellsReferences:
      return TraversalSelector::generateTraversal<ReferenceParticleCell<Particle_T>, Functor_T>(
          config.traversal, functor, traversalInfo, config.dataLayout, config.newton3);
    default:
      return TraversalSelector::generateTraversal<FullParticleCell<Particle_T>, Functor_T>(
          config.traversal, functor, traversalInfo, config.dataLayout, config.newton3);
  }
}
}  // namespace autopas
