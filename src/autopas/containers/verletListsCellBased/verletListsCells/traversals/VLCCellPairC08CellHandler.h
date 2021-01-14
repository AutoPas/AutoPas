//
// Created by TinaVl on 1/11/2021.
//
#pragma once
#include "autopas/options/DataLayoutOption.h"
#include "autopas/containers/linkedCells/traversals/LCC08CellHandler.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/VLCCellPairNeighborList.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCCellPairTraversalInterface.h"

namespace autopas
{
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class VLCCellPairC08CellHandler : public LCC08CellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>
{
 public:
  explicit VLCCellPairC08CellHandler(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                   double interactionLength, const std::array<double, 3> &cellLength,
                                     const std::array<unsigned long, 3> &overlap)
      : LCC08CellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(pairwiseFunctor, dims, interactionLength, cellLength, overlap)
  {
            this->computeOffsets(dims);
  }

  void processCellListsC08(VLCCellPairNeighborList<typename ParticleCell::ParticleType> &neighborList, unsigned long cellIndex,
      PairwiseFunctor *pairwiseFunctor, DataLayoutOption::Value layout, SoA<typename ParticleCell::ParticleType::SoAArraysType>* _soa,
                           std::array<unsigned long, 3> dims)
  {
    auto &aosNeighborList = neighborList.getAoSNeighborList();
    auto &soaNeighborList = neighborList.getSoANeighborList();
    auto &globalToLocalIndex = neighborList.getGlobalToLocalMap();

    for(auto &offsets : this->_cellPairOffsets){
      auto offsetCell1 = cellIndex + std::get<0>(offsets);
      auto offsetCell2 = cellIndex + std::get<1>(offsets);

      //the lists are built with a c18 traversal and the interaction will always be saved in the smaller cell's neighbor list
      /*if(offsetCell1 > offsetCell2)
      {
        auto temp = offsetCell1;
        offsetCell1 = offsetCell2;
        offsetCell2 = temp;
      }*/

      auto cell2Local = globalToLocalIndex[offsetCell1].find(offsetCell2);

      if(cell2Local != globalToLocalIndex[offsetCell1].end() && offsetCell1 != offsetCell2)
      {
        if(dataLayout == DataLayoutOption::aos)
        {
          //vector of pairs {particle, list}
          auto &currentList = aosNeighborList[offsetCell1][cell2Local->second];
          for(auto &[particleBasePtr, particleList] : currentList) {
            auto &particleBase = *particleBasePtr;
            for (auto particlePartnerPtr : particleList) {
              auto &particlePartner = *particlePartnerPtr;
              pairwiseFunctor->AoSFunctor(particleBase, particlePartner, useNewton3);
            }
          }
        }

        else if(dataLayout == DataLayoutOption::soa)
        {
          //vector of pairs {particle, list}
          auto &currentList = soaNeighborList[offsetCell1][cell2Local->second];
          for(auto &[particle, particleList] : currentList) {
            pairwiseFunctor->SoAFunctorVerlet(*_soa, particle, particleList, useNewton3);
          }
        }
      }

      cell2Local = globalToLocalIndex[offsetCell2].find(offsetCell1);
      if(cell2Local != globalToLocalIndex[offsetCell2].end())
      {
        if(dataLayout == DataLayoutOption::aos)
        {
          //vector of pairs {particle, list}
          auto &currentList = aosNeighborList[offsetCell2][cell2Local->second];
          for(auto &[particleBasePtr, particleList] : currentList) {
            auto &particleBase = *particleBasePtr;
            for (auto particlePartnerPtr : particleList) {
              auto &particlePartner = *particlePartnerPtr;
              pairwiseFunctor->AoSFunctor(particleBase, particlePartner, useNewton3);
              //newton3 off: do the cellOffsets have the pairs twice or I have to do it in both directions
              // or just tell the AoSFunctor and it will deal with it
            }
          }
        }

        else if(dataLayout == DataLayoutOption::soa)
        {
          //vector of pairs {particle, list}
          auto &currentList = soaNeighborList[offsetCell2][cell2Local->second];
          for(auto &[particle, particleList] : currentList) {
            pairwiseFunctor->SoAFunctorVerlet(*_soa, particle, particleList, useNewton3);
          }
        }
      }
    }
  }

};
}
