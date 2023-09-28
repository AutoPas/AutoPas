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
   * @tparam useSorting For traversals that use the CellFunctor: if the CellFunctor should apply sorting of particles
   * @param traversalType
   * @param pairwiseFunctor
   * @param info
   * @return Smartpointer to the traversal.
   */
  template <class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3, bool useSorting = true>
  static std::unique_ptr<TraversalInterface> generateTraversal(TraversalOption traversalType,
                                                               PairwiseFunctor &pairwiseFunctor,
                                                               const TraversalSelectorInfo &info);

  /**
   * Generates a given Traversal for the given properties.
   * Requires less templates but calls the templated version after a decision tree.
   * @tparam PairwiseFunctor
   * @param traversalType
   * @param pairwiseFunctor
   * @param info
   * @param dataLayout
   * @param useNewton3
   * @param useSorting For traversals that use the CellFunctor: if the CellFunctor should apply sorting of particles
   * @return Smartpointer to the traversal.
   */
  template <class PairwiseFunctor>
  static std::unique_ptr<TraversalInterface> generateTraversal(TraversalOption traversalType,
                                                               PairwiseFunctor &pairwiseFunctor,
                                                               const TraversalSelectorInfo &info,
                                                               DataLayoutOption dataLayout, Newton3Option useNewton3,
                                                               bool useSorting = true);
};

template <class ParticleCell>
template <class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3, bool useSorting>
std::unique_ptr<TraversalInterface> TraversalSelector<ParticleCell>::generateTraversal(
    TraversalOption traversalType, PairwiseFunctor &pairwiseFunctor, const TraversalSelectorInfo &info) {
  switch (traversalType) {
    // Direct sum
    case TraversalOption::ds_sequential: {
      return std::make_unique<DSSequentialTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, useSorting>>(
          &pairwiseFunctor,
          info.interactionLength /*this is the cutoff, as generated by DirectSum::getTraversalSelectorInfo()!*/);
    }
    // Linked cell
    case TraversalOption::lc_sliced: {
      return std::make_unique<LCSlicedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, useSorting>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::lc_sliced_c02: {
      return std::make_unique<LCSlicedC02Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, useSorting>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::lc_sliced_balanced: {
      return std::make_unique<
          LCSlicedBalancedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, useSorting>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::lc_c01: {
      return std::make_unique<LCC01Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, false, useSorting>>(
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
      return std::make_unique<LCC04Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, useSorting>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::lc_c04_HCP: {
      return std::make_unique<LCC04HCPTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, useSorting>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::lc_c08: {
      return std::make_unique<LCC08Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, useSorting>>(
          info.cellsPerDim, &pairwiseFunctor, info.interactionLength, info.cellLength);
    }
    case TraversalOption::lc_c18: {
      return std::make_unique<LCC18Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, useSorting>>(
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
      return std::make_unique<OTC18Traversal<ParticleType, PairwiseFunctor, dataLayout, useNewton3, useSorting>>(
          &pairwiseFunctor, info.interactionLength, info.interactionLength);
    }

    case TraversalOption::ot_c01: {
      using ParticleType = typename ParticleCell::ParticleType;
      return std::make_unique<OTC01Traversal<ParticleType, PairwiseFunctor, dataLayout, useNewton3, useSorting>>(
          &pairwiseFunctor, info.interactionLength, info.interactionLength);
    }
  }
  autopas::utils::ExceptionHandler::exception("Traversal type {} is not a known type!", traversalType.to_string());
  return {nullptr};
}
template <class ParticleCell>
template <class PairwiseFunctor>
std::unique_ptr<TraversalInterface> TraversalSelector<ParticleCell>::generateTraversal(
    TraversalOption traversalType, PairwiseFunctor &pairwiseFunctor, const TraversalSelectorInfo &traversalInfo,
    DataLayoutOption dataLayout, Newton3Option newton3, bool useSorting) {
  switch (dataLayout) {
    case DataLayoutOption::aos: {
      if (newton3 == Newton3Option::enabled) {
        if (useSorting) {
          return TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayoutOption::aos,
                                                                             true, true>(traversalType, pairwiseFunctor,
                                                                                         traversalInfo);
        } else {
          return TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayoutOption::aos,
                                                                             true, false>(
              traversalType, pairwiseFunctor, traversalInfo);
        }

      } else {
        if (useSorting) {
          return TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayoutOption::aos,
                                                                             false, true>(
              traversalType, pairwiseFunctor, traversalInfo);
        } else {
          return TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayoutOption::aos,
                                                                             false, false>(
              traversalType, pairwiseFunctor, traversalInfo);
        }
      }
    }
    case DataLayoutOption::soa: {
      if (newton3 == Newton3Option::enabled) {
        if (useSorting) {
          return TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayoutOption::soa,
                                                                             true, true>(traversalType, pairwiseFunctor,
                                                                                         traversalInfo);
        } else {
          return TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayoutOption::soa,
                                                                             true, false>(
              traversalType, pairwiseFunctor, traversalInfo);
        }

      } else {
        if (useSorting) {
          return TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayoutOption::soa,
                                                                             false, true>(
              traversalType, pairwiseFunctor, traversalInfo);

        } else {
          return TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayoutOption::soa,
                                                                             false, false>(
              traversalType, pairwiseFunctor, traversalInfo);
        }
      }
    }
  }

  autopas::utils::ExceptionHandler::exception("Traversal type {} is not a known type!", traversalType.to_string());
  return {nullptr};
}
}  // namespace autopas
