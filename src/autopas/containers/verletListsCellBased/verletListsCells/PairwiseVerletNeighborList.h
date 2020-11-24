/**
 * @file PairwiseVerletNeighborList.h
 * @author tirgendetwas
 * @date 03.11.20 od. 07.11.20
 */

#pragma once
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsNeighborListInterface.h"
#include "autopas/selectors/TraversalSelector.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCTraversalInterface.h"

namespace autopas
{
template <class Particle>
class PairwiseVerletNeighborList : public VerletListsCellsNeighborListInterface<Particle>
{
  using SoAPairOfParticleAndList = std::pair<size_t, std::vector<size_t, autopas::AlignedAllocator<size_t>>>;

 public:
  PairwiseVerletNeighborList(): _aosNeighborList{}, _particleToCellMap{}, _globalToLocalIndex{}, _soaNeighborList{} {}

  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::pairwiseVerletLists; }

  const size_t &getNumberOfPartners(const Particle *particle) const override {
    size_t localSize = 0;
    const auto [firstCellIndex, particleInCellIndex] = _particleToCellMap.at(const_cast<Particle *>(particle));
    for(size_t secondCellIndex = 0; secondCellIndex < 27; secondCellIndex++)
    {
      localSize+=_aosNeighborList[firstCellIndex][secondCellIndex][particleInCellIndex].second.size();
    }
    //auto toReturnSize = std::make_unique<size_t>(localSize);
    //auto ret = *(toReturnSize);
    return localSize;
  }

  typename VerletListsCellsHelpers<Particle>::PairwiseNeighborListsType &getAoSNeighborList() { return _aosNeighborList; }

  auto doCast(TraversalInterface *traversal)
  {
    return dynamic_cast<autopas::VLCTraversalInterface<Particle, typename VerletListsCellsHelpers<Particle>::PairwiseNeighborListsType> *>(traversal);
  }

  void buildAoSNeighborList(LinkedCells<Particle> &linkedCells,
                            bool useNewton3, double cutoff, double skin, double interactionLength,
                            const TraversalOption buildTraversalOption) override
  {
    _aosNeighborList.clear();
    auto &cells = linkedCells.getCells();
    auto cellsSize = cells.size();
    _aosNeighborList.resize(cellsSize);
    _globalToLocalIndex.resize(cellsSize);
    for(size_t firstCellIndex = 0; firstCellIndex < cellsSize; ++firstCellIndex)
    {
      _aosNeighborList[firstCellIndex].resize(27); //every cell has max 26 neighbors
      for(size_t secondCellIndex = 0; secondCellIndex < 27; ++secondCellIndex)
      {
        _aosNeighborList[firstCellIndex][secondCellIndex].reserve(cells[firstCellIndex].numParticles());
        size_t particleIndexCurrentCell = 0;
        for(auto particleIter = cells[firstCellIndex].begin(); particleIter.isValid(); ++particleIter)
        {
          Particle* currentParticle = &*particleIter;
          _aosNeighborList[firstCellIndex][secondCellIndex].emplace_back(std::make_pair(currentParticle, std::vector<Particle *>()));
          //using that magic number 5 from VLCneighborlist
          _aosNeighborList[firstCellIndex][secondCellIndex].back().second.reserve(cells[firstCellIndex].numParticles() * 5);
          _particleToCellMap[currentParticle] = std::make_pair(firstCellIndex, particleIndexCurrentCell);
          particleIndexCurrentCell++;
        }
      }
    }

    applyBuildFunctor(linkedCells, useNewton3, cutoff, skin, interactionLength, buildTraversalOption);
    generateSoAFromAoS(linkedCells);
  }

  void generateSoAFromAoS(LinkedCells<Particle> &linkedCells)
  {
    _soaNeighborList.clear();

    std::unordered_map<Particle*, size_t> _particleToIndex;
    _particleToIndex.reserve(linkedCells.getNumParticles());
    size_t i = 0;
    for (auto iter = linkedCells.begin(IteratorBehavior::haloOwnedAndDummy); iter.isValid(); ++iter, ++i) {
      _particleToIndex[&(*iter)] = i;
    }

    auto &cells = linkedCells.getCells();
    size_t cellsSize = cells.size();
    _soaNeighborList.resize(cellsSize);

    for(size_t firstCellIndex = 0; firstCellIndex < cellsSize; ++firstCellIndex) {
      _soaNeighborList[firstCellIndex].resize(27);  // every cell has max 26 neighbors

      for (size_t secondCellIndex = 0; secondCellIndex < 27; ++secondCellIndex) {
        _soaNeighborList[firstCellIndex][secondCellIndex].reserve(cells[firstCellIndex].numParticles());
        size_t particleIndexCurrentCell = 0;

        for (auto particleIter = cells[firstCellIndex].begin(); particleIter.isValid(); ++particleIter) {
          Particle *currentParticle = &*particleIter;
          size_t currentIndex = _particleToIndex.at(currentParticle);

          _soaNeighborList[firstCellIndex][secondCellIndex].emplace_back(std::make_pair(currentIndex, std::vector<size_t, autopas::AlignedAllocator<size_t>>()));

          for(size_t neighbor = 0; neighbor < _aosNeighborList[firstCellIndex][secondCellIndex][currentIndex].second.size(); neighbor++)
          {
            _soaNeighborList[firstCellIndex][secondCellIndex][currentIndex].second.emplace_back(_particleToIndex.at(_aosNeighborList[firstCellIndex][secondCellIndex][currentIndex].second[neighbor]));
          }
        }
      }
    }
  }


 private:

  void applyBuildFunctor(LinkedCells<Particle> &linkedCells, bool useNewton3, double cutoff, double skin,
                         double interactionLength, const TraversalOption buildTraversalOption) {
    typename VerletListsCellsHelpers<Particle>::PairwiseVerletListGeneratorFunctor f(_aosNeighborList, _particleToCellMap,
                                                                             _globalToLocalIndex,
                                                                             cutoff + skin);

    // Generate the build traversal with the traversal selector and apply the build functor with it.
    TraversalSelector<FullParticleCell<Particle>> traversalSelector;
    // Argument "cluster size" does not matter here.
    TraversalSelectorInfo traversalSelectorInfo(linkedCells.getCellBlock().getCellsPerDimensionWithHalo(),
                                                interactionLength, linkedCells.getCellBlock().getCellLength(), 0);
    autopas::utils::withStaticBool(useNewton3, [&](auto n3) {
      auto buildTraversal = traversalSelector.template generateTraversal<decltype(f), DataLayoutOption::soa, n3>( ///beware soa
          TraversalOption::lc_c08, f, traversalSelectorInfo); ///is it the correct traversal and is this how it works?
      linkedCells.iteratePairwise(buildTraversal.get());
    });
  }
  /**
   * Internal neighbor list structure in AoS format - Verlet lists for each particle for each cell.
   */
  typename VerletListsCellsHelpers<Particle>::PairwiseNeighborListsType _aosNeighborList;

  /**
   * Mapping of each particle to its corresponding cell and id within this cell.
   */
  std::unordered_map<Particle *, std::pair<size_t, size_t>> _particleToCellMap;
  std::vector<std::unordered_map<size_t, size_t>> _globalToLocalIndex;

  std::vector<std::vector<std::vector<SoAPairOfParticleAndList>>> _soaNeighborList;
  SoA<typename Particle::SoAArraysType> _soa;

  //LinkedCells<Particle> _internalLinkedCells;
};
}
