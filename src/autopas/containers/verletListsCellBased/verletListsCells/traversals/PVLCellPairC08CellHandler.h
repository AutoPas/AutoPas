/**
 * @file VLCCellPairC08CellHandler.h
 * @author tirgendetwas
 * @date 11.01.21
 */
#pragma once
#include "autopas/cells/SortedCellView.h"
#include "autopas/containers/linkedCells/traversals/LCC08CellHandler.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/PVLCellPairNeighborList.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/PVLCellPairTraversalInterface.h"
#include "autopas/iterators/SingleCellIterator.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/particles/Particle.h"

namespace autopas {
/**
 * This class provides the base logic for the c08 traversal for VLCCellPairNeighborList.
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class Particle, class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class PVLCellPairC08CellHandler : public LCC08CellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3> {
 public:
  /**
   * Constructor of the PVLCellPairC08CellHandler. The offsets for the base step are computed here using the
   * implementation from the superclass.
   * @param dims The number of cells per dimension.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param overlap number of overlapping cells in each direction as result from cutoff and cellLength.
   */
  explicit PVLCellPairC08CellHandler(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
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
  void processCellListsC08(PVLCellPairNeighborList<typename ParticleCell::ParticleType> &neighborList,
                           std::vector<ParticleCell> &cells,
                           unsigned long cellIndex, PairwiseFunctor *pairwiseFunctor, DataLayoutOption::Value layout,
                           SoA<typename ParticleCell::ParticleType::SoAArraysType> *_soa,
                           std::array<unsigned long, 3> dims) {
    auto &aosNeighborList = neighborList.getAoSNeighborList();
    auto &soaNeighborList = neighborList.getSoANeighborList();
    auto &globalToLocalIndex = neighborList.getGlobalToLocalMap();


    for (auto &offsets : this->_cellPairOffsets) {
      // the lists are built with a c18 traversal
      // the interaction will always be saved in the smaller cell's neighbor list
      const auto [offsetCell1, offsetCell2] =
          std::minmax({cellIndex + std::get<0>(offsets), cellIndex + std::get<1>(offsets)});

      auto cell2Local = globalToLocalIndex[offsetCell1].find(offsetCell2);

      // check if cell2 exists (in cell1's neighbor list)
      if (cell2Local != globalToLocalIndex[offsetCell1].end()) {
        std::array<double ,3> sortingDirection = std::get<2>(offsets);
        //create the view of cell2 needed for the 3D distance calculations
        SortedCellView<Particle, ParticleCell> cell2Sorted(cells[offsetCell2], sortingDirection);
        // if aos, send every pair of interacting particles to functor
        if (dataLayout == DataLayoutOption::aos) {
          // vector of pairs {particle, particle}
          auto &currentList = aosNeighborList[offsetCell1][cell2Local->second];
          for (auto &[particleBasePtr, particleList] : currentList) {
            auto &particleBase = *particleBasePtr;
            //the list contains only 1 particle anyway
            for (auto entryPVL : particleList) {
              //loop through the view of cell2 and perform the 3D distance calculation for every particle up until the PVL entry
              for (auto &cell2ParticlePtr : cell2Sorted._particles){
                Particle &cell2Particle = *cell2ParticlePtr.second;
                pairwiseFunctor->AoSFunctor(particleBase, cell2Particle, useNewton3);
                //if the last particle taken from cell2 is the PVL entry, stop the iteration
                if (entryPVL = cell2ParticlePtr.second){
                  break;
                }
              }
            }
          }
        }

        // if soa, send particle and corresponding neighbor list to the functor
        /*else if (dataLayout == DataLayoutOption::soa) {
          // vector of pairs {particle, list}
          auto &currentList = soaNeighborList[offsetCell1][cell2Local->second];
          for (auto &[particle, particleList] : currentList) {
            if (!particleList.empty()) {
              pairwiseFunctor->SoAFunctorVerlet(*_soa, particle, particleList, useNewton3);
            }
          }
        }*/
      }

      // if newton3 is off, find cell1 in cell2's neighbor list ("switch" the pair from above) and repeat the
      // interaction from above
      if (not useNewton3) {
        cell2Local = globalToLocalIndex[offsetCell2].find(offsetCell1);
        // exclude interaction within same cell - already handled in
        if (cell2Local != globalToLocalIndex[offsetCell2].end() and offsetCell1 != offsetCell2) {
          //create the view of cell2Local (the original cell1) needed for the 3D distance calculations
          std::array<double ,3> sortingDirection = {-std::get<2>(offsets)[0], -std::get<2>(offsets)[1],-std::get<2>(offsets)[2]};
          SortedCellView<Particle, ParticleCell> cell2LocalSorted(cells[offsetCell1], sortingDirection);
          // if aos, send every pair of interacting particles to functor
          if (dataLayout == DataLayoutOption::aos) {
            // vector of pairs {particle, particle}
            auto &currentList = aosNeighborList[offsetCell2][cell2Local->second];
            for (auto &[particleBasePtr, particleList] : currentList) {
              auto &particleBase = *particleBasePtr;
              //the list contains only 1 particle anyway
              for (auto entryPVL : particleList) {
                //loop through the view of cell2Local and perform the 3D distance calculation for every particle up until the PVL entry
                for (auto &cell2LocalParticlePtr : cell2LocalSorted._particles) {
                  Particle &cell2LocalParticle = *cell2LocalParticlePtr.second;
                  pairwiseFunctor->AoSFunctor(particleBase, cell2LocalParticle, useNewton3);
                  // if the last particle taken from cell2 is the PVL entry, stop the iteration
                  if (entryPVL = cell2LocalParticlePtr.second) {
                    break;
                  }
                }
              }
            }
          }

          // if soa, send particle and corresponding neighbor list to the functor
          /*else if (dataLayout == DataLayoutOption::soa) {
            // vector of pairs {particle, list}
            auto &currentList = soaNeighborList[offsetCell2][cell2Local->second];
            for (auto &[particle, particleList] : currentList) {
              if (not particleList.empty()) {
                pairwiseFunctor->SoAFunctorVerlet(*_soa, particle, particleList, useNewton3);
              }
            }
          }*/
        }
      }
    }
  }
};
}  // namespace autopas
