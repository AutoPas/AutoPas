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
   * Computes a base c08 step for the base cell with index cellIndex. The offsets computed in the constructor are used
   * to define the pairs of interacting cells.
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
    auto &aosNeighborList = neighborList.getAoSNeighborList();
    auto &soaNeighborList = neighborList.getSoANeighborList();
    auto &globalToLocalIndex = neighborList.getGlobalToLocalMap();

    for (auto &offsets : this->_cellPairOffsets) {
      auto offsetCell1 = cellIndex + std::get<0>(offsets);
      auto offsetCell2 = cellIndex + std::get<1>(offsets);

      // the lists are built with a c18 traversal and the interaction will always be saved in the smaller cell's
      // neighbor list
      if (offsetCell1 > offsetCell2) {
        auto temp = offsetCell1;
        offsetCell1 = offsetCell2;
        offsetCell2 = temp;
      }

      auto cell2Local = globalToLocalIndex[offsetCell1].find(offsetCell2);

      // check if cell2 exists (in cell1's neighbor list)
      if (cell2Local != globalToLocalIndex[offsetCell1].end()) {
        // if aos, send every pair of interacting particles to functor
        if (dataLayout == DataLayoutOption::aos) {
          // vector of pairs {particle, list}
          auto &currentList = aosNeighborList[offsetCell1][cell2Local->second];
          for (auto &[particleBasePtr, particleList] : currentList) {
            auto &particleBase = *particleBasePtr;
            for (auto particlePartnerPtr : particleList) {
              auto &particlePartner = *particlePartnerPtr;
              pairwiseFunctor->AoSFunctor(particleBase, particlePartner, useNewton3);
            }
          }
        }

        // if soa, send particle and corresponding neighbor list to the functor
        else if (dataLayout == DataLayoutOption::soa) {
          // vector of pairs {particle, list}
          auto &currentList = soaNeighborList[offsetCell1][cell2Local->second];
          for (auto &[particle, particleList] : currentList) {
            if (!particleList.empty()) {
              pairwiseFunctor->SoAFunctorVerlet(*_soa, particle, particleList, useNewton3);
            }
          }
        }
      }

      // if newton3 is off, find cell1 in cell2's neighbor list ("switch" the pair from above) and repeat the
      // interaction from above
      if (useNewton3 == false) {
        cell2Local = globalToLocalIndex[offsetCell2].find(offsetCell1);
        // exclude interaction within same cell - already handled in
        if (cell2Local != globalToLocalIndex[offsetCell2].end() && offsetCell1 != offsetCell2) {
          // if aos, send every pair of interacting particles to functor
          if (dataLayout == DataLayoutOption::aos) {
            // vector of pairs {particle, list}
            auto &currentList = aosNeighborList[offsetCell2][cell2Local->second];
            for (auto &[particleBasePtr, particleList] : currentList) {
              auto &particleBase = *particleBasePtr;
              for (auto particlePartnerPtr : particleList) {
                auto &particlePartner = *particlePartnerPtr;
                pairwiseFunctor->AoSFunctor(particleBase, particlePartner, useNewton3);
              }
            }
          }

          // if soa, send particle and corresponding neighbor list to the functor
          else if (dataLayout == DataLayoutOption::soa) {
            // vector of pairs {particle, list}
            auto &currentList = soaNeighborList[offsetCell2][cell2Local->second];
            for (auto &[particle, particleList] : currentList) {
              if (!particleList.empty()) {
                pairwiseFunctor->SoAFunctorVerlet(*_soa, particle, particleList, useNewton3);
              }
            }
          }
        }
      }
    }
  }
};
}  // namespace autopas
