#pragma once

#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsHelpers.h"

namespace autopas
{
template<class Particle>
class VerletListsCellsNeighborList
{
 public:
  VerletListsCellsNeighborList() : _aosNeighborList{}, _particleToCellMap{} {}
  void buildAoSNeighborList(LinkedCells<typename VerletListsCellsHelpers<Particle>::VLCCellType> &linkedCells, bool useNewton3)
  {
    // Initialize a neighbor list for each cell.
    _aosNeighborList.clear();
    auto &cells = linkedCells.getCells();
    size_t cellsSize = cells.size();
    _aosNeighborList.resize(cellsSize);
    for (size_t cellIndex = 0; cellIndex < cellsSize; ++cellIndex) {
      _aosNeighborList[cellIndex].reserve(cells[cellIndex].numParticles());
      size_t particleIndexWithinCell = 0;
      for (auto iter = cells[cellIndex].begin(); iter.isValid(); ++iter, ++particleIndexWithinCell) {
        Particle *particle = &*iter;
        _aosNeighborList[cellIndex].emplace_back(particle, std::vector<Particle *>());
        // In a cell with N particles, reserve space for 5N neighbors.
        // 5 is an empirically determined magic number that provides good speed.
        _aosNeighborList[cellIndex].back().second.reserve(cells[cellIndex].numParticles() * 5);
        _particleToCellMap[particle] = std::make_pair(cellIndex, particleIndexWithinCell);
      }
    }
  }

  auto &getParticleToCellMap() {return _particleToCellMap;}
  auto &getParticleToCellMapConst() const {return _particleToCellMap;}

  typename VerletListsCellsHelpers<Particle>::NeighborListsType &getAoSNeighborList() {return _aosNeighborList;}
  //typename VerletListsCellsHelpers<Particle>::NeighborListsType &getAoSConst() const {return _aosNeighborList;}
  //const typename VerletListsCellsHelpers<Particle>::NeighborListsType &getAoSConstRef() {return _aosNeighborList;}
  const typename VerletListsCellsHelpers<Particle>::NeighborListsType &getAoSConstAll() const {return _aosNeighborList;}
  //typename VerletListsCellsHelpers<Particle>::NeighborListsType getAoS() {return _aosNeighborList;}
  //typename VerletListsCellsHelpers<Particle>::NeighborListsType getAoSBlandConst() const {return _aosNeighborList;}

 public:
  typename VerletListsCellsHelpers<Particle>::NeighborListsType _aosNeighborList;
  std::unordered_map<Particle *, std::pair<size_t, size_t>> _particleToCellMap;

};
}
