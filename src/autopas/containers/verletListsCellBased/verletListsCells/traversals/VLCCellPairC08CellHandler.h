/**
 * @file VLCCellPairC08CellHandler.h
 * @author tirgendetwas
 * @date 11.01.21
 */
#pragma once

#include "autopas/containers/linkedCells/traversals/LCC08CellHandlerUtility.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/VLCCellPairNeighborList.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCCellPairTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {
/**
 * This class provides the base logic for the c08 traversal for VLCCellPairNeighborList.
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 */
template <class ParticleCell, class PairwiseFunctor>
class VLCCellPairC08CellHandler {
 public:
  /**
   * Constructor of the VLCCellPairC08CellHandler. The offsets for the base step are computed here using the
   * implementation from the superclass.
   * @param dims The number of cells per dimension.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   */
  VLCCellPairC08CellHandler(const std::array<unsigned long, 3> &dims, double interactionLength,
                            const std::array<double, 3> &cellLength)
      : _cellPairOffsets{LCC08CellHandlerUtility::computePairwiseCellOffsetsC08<
            LCC08CellHandlerUtility::C08OffsetMode::c08CellPairs>(dims, cellLength, interactionLength)} {}

  /**
   * Executes a c08 base step for the cell at cellIndex.
   * The offsets computed in the constructor are used to identify the pairs of interacting cells.
   * @param neighborList pairs of main particle and partner particles for each pair of cells
   * @param cellIndex index of the current base cell
   * @param pairwiseFunctor functor that defines the interaction of two particles
   * @param layout data layout (AoS or SoA)
   * @param soa Structure of Arrays where the particles are loaded if the SoA data layout is used
   * @param useNewton3
   */
  void processCellListsC08(VLCCellPairNeighborList<typename ParticleCell::ParticleType> &neighborList,
                           unsigned long cellIndex, PairwiseFunctor *pairwiseFunctor, DataLayoutOption layout,
                           SoA<typename ParticleCell::ParticleType::SoAArraysType> *soa, bool useNewton3) {
    const auto &aosNeighborList = neighborList.getAoSNeighborList();
    const auto &soaNeighborList = neighborList.getSoANeighborList();
    const auto &globalToLocalIndex = neighborList.getGlobalToLocalMap();

    // for all interaction pairs defined via the c08 base step
    for (const auto &[offsetA, offsetB] : this->_cellPairOffsets) {
      // the lists are built with a c18 traversal
      // the interaction will always be saved in the smaller cell's neighbor list
      // std::minmax(a, b) returns references. Hence, we can't use temporaries as arguments.
      const auto offsetCellA = cellIndex + offsetA;
      const auto offsetCellB = cellIndex + offsetB;
      const auto [offsetCell1, offsetCell2] = std::minmax(offsetCellA, offsetCellB);

      const auto cell2Local = globalToLocalIndex[offsetCell1].find(offsetCell2);

      // check if cell2 exists (in cell1's neighbor list)
      if (cell2Local != globalToLocalIndex[offsetCell1].end()) {
        // if aos, send every pair of interacting particles to functor
        if (layout == DataLayoutOption::aos) {
          // vector of pairs {particle, list}
          const auto &currentList = aosNeighborList[offsetCell1][cell2Local->second];
          for (auto &[particleBasePtr, particleList] : currentList) {
            for (auto *particlePartnerPtr : particleList) {
              pairwiseFunctor->AoSFunctor(*particleBasePtr, *particlePartnerPtr, useNewton3);
            }
          }
        }

        // if soa, send particle and corresponding neighbor list to the functor
        else if (layout == DataLayoutOption::soa) {
          // vector of pairs {particle, list}
          const auto &currentList = soaNeighborList[offsetCell1][cell2Local->second];
          for (const auto &[particleIndex, particleList] : currentList) {
            if (not particleList.empty()) {
              pairwiseFunctor->SoAFunctorVerlet(*soa, particleIndex, particleList, useNewton3);
            }
          }
        }
      }

      // if newton3 is off, find cell1 in cell2's neighbor list ("switch" the pair from above) and repeat the
      // interaction from above
      if (not useNewton3) {
        const auto cell2LocalNoN3 = globalToLocalIndex[offsetCell2].find(offsetCell1);
        // exclude interaction within same cell - already handled in
        if (cell2LocalNoN3 != globalToLocalIndex[offsetCell2].end() and offsetCell1 != offsetCell2) {
          // if aos, send every pair of interacting particles to functor
          if (layout == DataLayoutOption::aos) {
            // vector of pairs {particle, list}
            const auto &currentList = aosNeighborList[offsetCell2][cell2LocalNoN3->second];
            for (auto &[particleBasePtr, particleList] : currentList) {
              for (auto *particlePartnerPtr : particleList) {
                pairwiseFunctor->AoSFunctor(*particleBasePtr, *particlePartnerPtr, useNewton3);
              }
            }
          }

          // if soa, send particle and corresponding neighbor list to the functor
          else if (layout == DataLayoutOption::soa) {
            // vector of pairs {particle, list}
            const auto &currentList = soaNeighborList[offsetCell2][cell2LocalNoN3->second];
            for (const auto &[particleIndex, particleList] : currentList) {
              if (not particleList.empty()) {
                pairwiseFunctor->SoAFunctorVerlet(*soa, particleIndex, particleList, useNewton3);
              }
            }
          }
        }
      }
    }
  }

 private:
  /**
   * Member containng the cell pair offsets for processCellListsC08
   */
  std::vector<LCC08CellHandlerUtility::OffsetPair> _cellPairOffsets;
};
}  // namespace autopas
