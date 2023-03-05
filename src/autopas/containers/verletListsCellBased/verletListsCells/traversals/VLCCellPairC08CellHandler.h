/**
 * @file VLCCellPairC08CellHandler.h
 * @author tirgendetwas
 * @date 11.01.21
 */
#pragma once
#include "autopas/containers/linkedCells/traversals/LCC08CellHandler.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/VLCCellPairNeighborList.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCCellPairTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {
/**
 * This class provides the base logic for the c08 traversal for VLCCellPairNeighborList.
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class VLCCellPairC08CellHandler : public LCC08CellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3> {
 public:
  /**
   * Constructor of the VLCCellPairC08CellHandler. The offsets for the base step are computed here using the
   * implementation from the superclass.
   * @param dims The number of cells per dimension.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param overlap number of overlapping cells in each direction as result from cutoff and cellLength.
   */
  explicit VLCCellPairC08CellHandler(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                     double interactionLength, const std::array<double, 3> &cellLength,
                                     const std::array<unsigned long, 3> &overlap)
      : LCC08CellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(
            pairwiseFunctor, dims, interactionLength, cellLength, overlap) {
    this->computeOffsets(dims);
  }

  /**
   * Executes a c08 base step for the cell at cellIndex.
   * The offsets computed in the constructor are used to identify the pairs of interacting cells.
   * @param neighborList pairs of main particle and partner particles for each pair of cells
   * @param cellIndex index of the current base cell
   * @param pairwiseFunctor functor that defines the interaction of two particles
   * @param layout data layout (AoS or SoA)
   * @param _soa Structure of Arrays where the particles are loaded if the SoA data layout is used
   * @param dims number of cells per dimension
   */
  void processCellListsC08(VLCCellPairNeighborList<typename ParticleCell::ParticleType> &neighborList,
                           unsigned long cellIndex, PairwiseFunctor *pairwiseFunctor, DataLayoutOption::Value layout,
                           SoA<typename ParticleCell::ParticleType::SoAArraysType> *_soa,
                           std::array<unsigned long, 3> dims) {
    const auto &aosNeighborList = neighborList.getAoSNeighborList();
    const auto &soaNeighborList = neighborList.getSoANeighborList();
    const auto &globalToLocalIndex = neighborList.getGlobalToLocalMap();

    // for all interaction pairs defined via the c08 base step
    for (const auto &offsets : this->_cellPairOffsets) {
      // the lists are built with a c18 traversal
      // the interaction will always be saved in the smaller cell's neighbor list
      const auto [offsetCell1, offsetCell2] =
          std::minmax(cellIndex + std::get<0>(offsets), cellIndex + std::get<1>(offsets));

      const auto cell2Local = globalToLocalIndex[offsetCell1].find(offsetCell2);

      // check if cell2 exists (in cell1's neighbor list)
      if (cell2Local != globalToLocalIndex[offsetCell1].end()) {
        // if aos, send every pair of interacting particles to functor
        if (dataLayout == DataLayoutOption::aos) {
          // vector of pairs {particle, list}
          const auto &currentList = aosNeighborList[offsetCell1][cell2Local->second];
          for (auto &[particleBasePtr, particleList] : currentList) {
            for (auto *particlePartnerPtr : particleList) {
              pairwiseFunctor->AoSFunctor(*particleBasePtr, *particlePartnerPtr, useNewton3);
            }
          }
        }

        // if soa, send particle and corresponding neighbor list to the functor
        else if (dataLayout == DataLayoutOption::soa) {
          // vector of pairs {particle, list}
          const auto &currentList = soaNeighborList[offsetCell1][cell2Local->second];
          for (const auto &[particleIndex, particleList] : currentList) {
            if (not particleList.empty()) {
              pairwiseFunctor->SoAFunctorVerlet(*_soa, particleIndex, particleList, useNewton3);
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
          if (dataLayout == DataLayoutOption::aos) {
            // vector of pairs {particle, list}
            const auto &currentList = aosNeighborList[offsetCell2][cell2LocalNoN3->second];
            for (auto &[particleBasePtr, particleList] : currentList) {
              for (auto *particlePartnerPtr : particleList) {
                pairwiseFunctor->AoSFunctor(*particleBasePtr, *particlePartnerPtr, useNewton3);
              }
            }
          }

          // if soa, send particle and corresponding neighbor list to the functor
          else if (dataLayout == DataLayoutOption::soa) {
            // vector of pairs {particle, list}
            const auto &currentList = soaNeighborList[offsetCell2][cell2LocalNoN3->second];
            for (const auto &[particleIndex, particleList] : currentList) {
              if (not particleList.empty()) {
                pairwiseFunctor->SoAFunctorVerlet(*_soa, particleIndex, particleList, useNewton3);
              }
            }
          }
        }
      }
    }
  }
};
}  // namespace autopas
