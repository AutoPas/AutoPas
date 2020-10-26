#pragma once

#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsHelpers.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/selectors/TraversalSelector.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StaticSelectorMacros.h"

namespace autopas
{
template<class Particle>
class VerletListsCellsNeighborList
{

  using LinkedParticleCell = typename VerletListsCellsHelpers<Particle>::VLCCellType;

 public:
  VerletListsCellsNeighborList() : _aosNeighborList{}, _particleToCellMap{} {}
  void buildAoSNeighborList(LinkedCells<typename VerletListsCellsHelpers<Particle>::VLCCellType> &linkedCells, bool useNewton3,
                            double cutoff, double skin, double interactionLength, const TraversalOption buildTraversalOption)
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

    applyBuildFunctor(linkedCells, useNewton3, cutoff, skin, interactionLength, buildTraversalOption);
  }

  auto &getParticleToCellMapConst() const {return _particleToCellMap;}

  typename VerletListsCellsHelpers<Particle>::NeighborListsType &getAoSNeighborList() {return _aosNeighborList;}
  const typename VerletListsCellsHelpers<Particle>::NeighborListsType &getAoSConstAll() const {return _aosNeighborList;}

 private:

  void applyBuildFunctor(LinkedCells<typename VerletListsCellsHelpers<Particle>::VLCCellType> &linkedCells, bool useNewton3,
                         double cutoff, double skin, double interactionLength, const TraversalOption buildTraversalOption)
  {
    typename VerletListsCellsHelpers<Particle>::VerletListGeneratorFunctor f(_aosNeighborList, _particleToCellMap,
                                                                             cutoff + skin);

    // Generate the build traversal with the traversal selector and apply the build functor with it.
    TraversalSelector<LinkedParticleCell> traversalSelector;
    // Argument "cluster size" does not matter here.
    TraversalSelectorInfo traversalSelectorInfo(linkedCells.getCellBlock().getCellsPerDimensionWithHalo(),
                                                interactionLength,
                                                linkedCells.getCellBlock().getCellLength(), 0);
    autopas::utils::withStaticBool(useNewton3, [&](auto n3) {
      auto buildTraversal = traversalSelector.template generateTraversal<decltype(f), DataLayoutOption::aos, n3>(
          buildTraversalOption, f, traversalSelectorInfo);
      linkedCells.iteratePairwise(buildTraversal.get());
    });
  }

  typename VerletListsCellsHelpers<Particle>::NeighborListsType _aosNeighborList;
  std::unordered_map<Particle *, std::pair<size_t, size_t>> _particleToCellMap;

};
}
