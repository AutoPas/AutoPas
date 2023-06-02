/**
* @file VLCPartialCellPairNeighborList.h
* @author Luis Gall
* @date 31.05.23
*/

#pragma once
#include "VLCCellPairNeighborList.h"

namespace autopas {


template <class Particle>
class VLCPartialCellPairNeighborList : public VLCCellPairNeighborList<Particle> {
public:

 [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::partialPairwiseVerletLists; }


 void buildAoSNeighborList(LinkedCells<Particle> &linkedCells, bool useNewton3, double cutoff, double skin,
                           double interactionLength, const TraversalOption buildTraversalOption,
                           typename VerletListsCellsHelpers<Particle>::VLCBuildType::Value buildType, bool partialRebuilding) override {
   this->_internalLinkedCells = &linkedCells;
   auto &cells = linkedCells.getCells();
   const auto cellsSize = cells.size();
   const auto cellLength = linkedCells.getCellBlock().getCellLength();
   const auto interactionLengthSquare = linkedCells.getInteractionLength() * linkedCells.getInteractionLength();

   if (!partialRebuilding) {
     this->_aosNeighborList.clear();
     this->_globalToLocalIndex.clear();
     this->_particleToCellMap.clear();

     this->_aosNeighborList.resize(cellsSize);
     this->_globalToLocalIndex.resize(cellsSize);
   }

   std::array<long, 3> overlap{};
   for (unsigned int d = 0; d < 3; d++) {
     overlap[d] = std::ceil(linkedCells.getInteractionLength() / cellLength[d]);
   }

   // count number of neighbor cells
   size_t neighborCells = 0;
   for (int x = -static_cast<int>(overlap[0]); x < overlap[0] + 1; x++) {
     std::array<double, 3> pos{};
     pos[0] = std::max(0l, (std::abs(x) - 1l)) * cellLength[0];
     for (int y = -static_cast<int>(overlap[1]); y < overlap[1] + 1; y++) {
       pos[1] = std::max(0l, (std::abs(y) - 1l)) * cellLength[1];
       for (int z = -static_cast<int>(overlap[2]); z < overlap[2] + 1; z++) {
         pos[2] = std::max(0l, (std::abs(z) - 1l)) * cellLength[2];
         const double distSquare = utils::ArrayMath::dot(pos, pos);
         if (distSquare <= interactionLengthSquare) {
           neighborCells++;
         }
       }
     }
   }

   // when N3 is used we only need half of the cells (rounded up)
   if (useNewton3) {
     neighborCells /= 2;
     if (neighborCells % 2 != 0) {
       neighborCells++;
     }
   }

   // initialize lists for every particle-cell pair
   for (size_t firstCellIndex = 0; firstCellIndex < cellsSize; ++firstCellIndex) {

     // either not partial rebuilding or the current cell is dirty
     if (!partialRebuilding || cells[firstCellIndex].getDirty() || cells[firstCellIndex].getInflowDirty() ||cells[firstCellIndex].getOutflowDirty()) {
       size_t numParticlesFirstCell = cells[firstCellIndex].numParticles();
       this->_aosNeighborList[firstCellIndex].clear();
       this->_aosNeighborList[firstCellIndex].resize(neighborCells);
       for (auto &cellPair : this->_aosNeighborList[firstCellIndex]) {
         // reserve vector of neighbor lists for every particle in cell1
         cellPair.reserve(numParticlesFirstCell);
         size_t particleIndexCurrentCell = 0;
         for (auto &particle : cells[firstCellIndex]) {
           // for each particle in cell1 make a pair of particle and neighbor list
           cellPair.emplace_back(std::make_pair(&particle, std::vector<Particle *>()));

           // add a pair of cell's index and particle's index in the cell
           this->_particleToCellMap[&particle] = std::make_pair(firstCellIndex, particleIndexCurrentCell);
           particleIndexCurrentCell++;
         }
       }
     }
     // partial rebuilding and cell is not dirty
     else {
       // have a look at this cells neighboring cells
       for (auto [globalIndex, localIndex] : this->_globalToLocalIndex.at(firstCellIndex)) {
         // neighboring cell is dirty
         if (cells.at(globalIndex).getDirty() || cells.at(globalIndex).getOutflowDirty() || cells.at(globalIndex).getInflowDirty()) {
           // clear all neighbor lists for this cell pair
           this->_aosNeighborList[firstCellIndex][localIndex].clear();
           size_t numParticlesFirstCell = cells[firstCellIndex].numParticles();
           this->_aosNeighborList[firstCellIndex][localIndex].reserve(numParticlesFirstCell);
           size_t particleIndexCurrentCell = 0;
           for (auto &particle : cells[firstCellIndex]) {
             this->_aosNeighborList[firstCellIndex][localIndex].emplace_back(std::make_pair(&particle, std::vector<Particle *>()));
             this->_particleToCellMap[&particle] = std::make_pair(firstCellIndex, particleIndexCurrentCell);
             particleIndexCurrentCell++;
           }
         }
       }
     }
   }

   // fill the lists
   linkedCells.setOnlyDirtyCells(partialRebuilding);
   this->applyBuildFunctor(linkedCells, useNewton3, cutoff, skin, interactionLength, buildTraversalOption, buildType, partialRebuilding);
 }
};
}  // namespace autopas
