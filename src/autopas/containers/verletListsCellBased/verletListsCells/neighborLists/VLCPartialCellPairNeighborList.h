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
     this->_particleToCellMap.reserve(linkedCells.getNumberOfParticles());
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
   // TODO : think about chunk size
#pragma omp parallel for schedule(dynamic, 15)
   for (size_t firstCellIndex = 0; firstCellIndex < cellsSize; ++firstCellIndex) {

     // either not partial rebuilding or the current cell is dirty
     if (!partialRebuilding || cells[firstCellIndex].getDirty() || cells[firstCellIndex].getInflowDirty() /*|| cells[firstCellIndex].getOutflowDirty()*/) {
       size_t numParticlesFirstCell = cells[firstCellIndex].numParticles();
       this->_aosNeighborList[firstCellIndex].clear();
       this->_aosNeighborList[firstCellIndex].resize(neighborCells);
       // this->_globalToLocalIndex[firstCellIndex].clear();
       // this->_globalToLocalIndex[firstCellIndex].reserve(neighborCells);
       for (auto &cellPair : this->_aosNeighborList[firstCellIndex]) {
         // reserve vector of neighbor lists for every particle in cell1
         cellPair.reserve(numParticlesFirstCell);
         size_t particleIndexCurrentCell = 0;
         for (auto &particle : cells[firstCellIndex]) {
           if (particle.isDummy()) {
             continue;
           }

           // for each particle in cell1 make a pair of particle and neighbor list
           cellPair.emplace_back(std::make_pair(&particle, std::vector<Particle *>{}));
           // add a pair of cell's index and particle's index in the cell
           // critical if a new entry is necessary
           if (this->_particleToCellMap.find(particle.getID()) == this->_particleToCellMap.end()) {
#pragma omp critical
             this->_particleToCellMap[particle.getID()] = std::make_pair(firstCellIndex, particleIndexCurrentCell);
           }
           else {
             this->_particleToCellMap[particle.getID()] = std::make_pair(firstCellIndex, particleIndexCurrentCell);
           }
           particleIndexCurrentCell++;
         }
       }
     }
     // partial rebuilding and cell is inflow dirty
     /*
     else if (cells.at(firstCellIndex).getInflowDirty()) {
       for (auto [globalIndex, localIndex] : this->_globalToLocalIndex.at(firstCellIndex)) {
         // only delete neighbour lists which are pointing to same cell
         if (globalIndex == firstCellIndex) {
           // clear only the neighbor list which relates the particles in the same cell
           this->_aosNeighborList[firstCellIndex][localIndex].clear();
           size_t numParticlesFirstCell = cells[firstCellIndex].numParticles();
           this->_aosNeighborList[firstCellIndex][localIndex].reserve(numParticlesFirstCell);
           size_t particleIndexCurrentCell = 0;
           for (auto &particle : cells[firstCellIndex]) {
             this->_aosNeighborList[firstCellIndex][localIndex].emplace_back(std::make_pair(&particle, std::vector<Particle *>{}));
             this->_particleToCellMap[&particle] = std::make_pair(firstCellIndex, particleIndexCurrentCell);
             particleIndexCurrentCell++;
           }
         }
       }
     }

     // partial rebuilding and base cell has outflow particles
     else if (cells[firstCellIndex].getOutflowDirty()) {
       // delete all dummy particles from the neighbor lists
       for (auto [globalIndex, localIndex] : this->_globalToLocalIndex.at(firstCellIndex)) {
         // neighbour list pointing to same cell
         if (firstCellIndex == globalIndex) {
           for (auto& [particle, neighborList] : this->_aosNeighborList[firstCellIndex][localIndex]) {
             // delete dummy particles' neighbor lists and replace them with real dummy
              if (particle->isDummy()) {
                auto& [cellIndex, particleIndex] = this->_particleToCellMap[particle];
                this->_particleToCellMap.erase(particle);
                this->_aosNeighborList[firstCellIndex][localIndex][particleIndex].second.clear();
                Particle dummyParticle {};
                dummyParticle.setOwnershipState(OwnershipState::dummy);
                this->_aosNeighborList[firstCellIndex][localIndex][particleIndex].first = &dummyParticle;
              }
              // delete every dummy particle in a neighbor list
              for (auto neighborIter = neighborList.begin(); neighborIter != neighborList.end() ;) {
                if ((*neighborIter)->isDummy()) {
                  neighborIter = neighborList.erase(neighborIter);
                }
                else {
                  ++neighborIter;
                }
              }
           }
         }
       }
     }
     */
     // partial rebuilding and base cell is not dirty
     else {
       // have a look at this cells neighboring cells
       for (auto [globalIndex, localIndex] : this->_globalToLocalIndex.at(firstCellIndex)) {
         // neighboring cell is dirty
         if (cells.at(globalIndex).getDirty() || cells.at(globalIndex).getInflowDirty() /*|| cells.at(globalIndex).getOutflowDirty()*/) {
           // clear all neighbor lists for this cell pair
           this->_aosNeighborList[firstCellIndex][localIndex].clear();
           size_t numParticlesFirstCell = cells[firstCellIndex].numParticles();
           this->_aosNeighborList[firstCellIndex][localIndex].reserve(numParticlesFirstCell);
           size_t particleIndexCurrentCell = 0;
           for (auto &particle : cells[firstCellIndex]) {
             if (particle.isDummy()) {
               continue ;
             }

             this->_aosNeighborList[firstCellIndex][localIndex].emplace_back(std::make_pair(&particle, std::vector<Particle *>{}));
             // critical if a new entry is necessary
             if (this->_particleToCellMap.find(particle.getID()) == this->_particleToCellMap.end()) {
#pragma omp critical
               this->_particleToCellMap[particle.getID()] = std::make_pair(firstCellIndex, particleIndexCurrentCell);
             }
             else {
               this->_particleToCellMap[particle.getID()] = std::make_pair(firstCellIndex, particleIndexCurrentCell);
             }
             particleIndexCurrentCell++;
           }
         }

         /*
         // if neighbour has outflowing particles -> delete them from the neighbour lists
         else if (cells.at(globalIndex).getOutflowDirty()) {
           for (auto particleNeighboursIter = this->_aosNeighborList[firstCellIndex][localIndex].begin();
                particleNeighboursIter != this->_aosNeighborList[firstCellIndex][localIndex].end(); ++particleNeighboursIter) {

             for (auto neighborIter = particleNeighboursIter->second.begin(); neighborIter != particleNeighboursIter->second.end();) {
                 if ((*neighborIter)->isDummy()) {
                   neighborIter = particleNeighboursIter->second.erase(neighborIter);
                 }
                 else {
                   ++neighborIter;
                 }
             }
           }
         }
         */
       }
     }
   }
   // std::cout << "Size of particle to cell map: " << this->_particleToCellMap.size() << std::endl;
   // fill the lists
   linkedCells.setOnlyDirtyCells(partialRebuilding);
   this->applyBuildFunctor(linkedCells, useNewton3, cutoff, skin, interactionLength, TraversalOption::lc_c01_b08, buildType, partialRebuilding);
 }
};
}  // namespace autopas
