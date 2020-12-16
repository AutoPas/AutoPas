/**
 * @file PairwiseVerletNeighborList.h
 * @author tirgendetwas
 * @date 03.11.20 od. 07.11.20
 */

#pragma once
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsNeighborListInterface.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCTraversalInterface.h"
#include "autopas/utils/StaticBoolSelector.h"

namespace autopas {
template <class ParticleCell>
class TraversalSelector;
/**
 * Neighbor list to be used with VerletListsCells container.
 * Pairwise verlet lists iterates through each pair of neighboring cells
 * and generates a neighbor list for each particle from cell1, which consists of its (potential) partners from cell2.
 * @tparam Particle Type of particle to be used for this neighbor list.
 * */
template <class Particle>
class PairwiseVerletNeighborList : public VerletListsCellsNeighborListInterface<Particle> {
 public:
  /**
   * Type of the data structure used to save the neighbor lists.
   * */
  using listType = typename VerletListsCellsHelpers<Particle>::PairwiseNeighborListsType;

  //TODO
  using SoAPairOfParticleAndList = std::pair<size_t, std::vector<size_t, autopas::AlignedAllocator<size_t>>>;

  //TODO
  using soaListType = typename std::vector<std::vector<std::vector<SoAPairOfParticleAndList>>>;

  /**
   * Constructor for PairwiseVerletNeighborList. Initializes private attributes.
   * */
  PairwiseVerletNeighborList() : _aosNeighborList{}, _particleToCellMap{}, _globalToLocalIndex{}, _soaNeighborList{}, _soa{} {}

  /**
   * @copydoc VerletListsCellsNeighborListInterface::getContainerType()
   * */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::pairwiseVerletLists; }

  /**
   * @copydoc VerletListsCellsNeighborListInterface::getNumberOfPartners()
   * */
  const size_t getNumberOfPartners(const Particle *particle) const override {
    size_t localSize = 0;
    const auto [firstCellIndex, particleInCellIndex] = _particleToCellMap.at(const_cast<Particle *>(particle));
    size_t numberOfCellsToInteract = 27;
    for (size_t secondCellIndex = 0; secondCellIndex < numberOfCellsToInteract; secondCellIndex++) {
      localSize += _aosNeighborList[firstCellIndex][secondCellIndex][particleInCellIndex].second.size();
    }
    return localSize;
  }

  /**
   * Returns the neighbor list in AoS layout.
   * @return Neighbor list in AoS layout.
   * */
  typename VerletListsCellsHelpers<Particle>::PairwiseNeighborListsType &getAoSNeighborList() {
    return _aosNeighborList;
  }

  auto &getSoANeighborList() {return _soaNeighborList;}

  /**
   * @copydoc VerletListsCellsNeighborListInterface::buildAoSNeighborList()
   * */
  void buildAoSNeighborList(LinkedCells<Particle> &linkedCells, bool useNewton3, double cutoff, double skin,
                            double interactionLength, const TraversalOption buildTraversalOption) override {
    _internalLinkedCells = &linkedCells;
    _aosNeighborList.clear();
    auto &cells = linkedCells.getCells();
    auto cellsSize = cells.size();
    _aosNeighborList.resize(cellsSize);
    _globalToLocalIndex.resize(cellsSize);

    for (size_t firstCellIndex = 0; firstCellIndex < cellsSize; ++firstCellIndex) {
      size_t numCellsInteracting = 27;  // every cell has max 26 neighbors + interaction with itself
      _aosNeighborList[firstCellIndex].resize(numCellsInteracting);
      for (size_t secondCellIndex = 0; secondCellIndex < numCellsInteracting; ++secondCellIndex) {
        // reserve vector of neighbor lists for every particle in cell1
        _aosNeighborList[firstCellIndex][secondCellIndex].reserve(cells[firstCellIndex].numParticles());
        size_t particleIndexCurrentCell = 0;
        for (auto particleIter = cells[firstCellIndex].begin(); particleIter.isValid(); ++particleIter) {
          // for each particle in cell1 make pair of particle and neighbor list
          Particle *currentParticle = &*particleIter;
          _aosNeighborList[firstCellIndex][secondCellIndex].emplace_back(
              std::make_pair(currentParticle, std::vector<Particle *>()));

          // magic number 5 doesn't make sense here anymore
          // how much should we actually reserve?
          _aosNeighborList[firstCellIndex][secondCellIndex].back().second.reserve(cells[firstCellIndex].numParticles() *
                                                                                  5 / 27);
          // add pair of cell's index and particle's index in the cell
          _particleToCellMap[currentParticle] = std::make_pair(firstCellIndex, particleIndexCurrentCell);
          particleIndexCurrentCell++;
        }
      }
    }

    applyBuildFunctorSoA(linkedCells, useNewton3, cutoff, skin, interactionLength, buildTraversalOption);
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

          for(size_t neighbor = 0; neighbor < _aosNeighborList[firstCellIndex][secondCellIndex][particleIndexCurrentCell].second.size(); neighbor++)
          {
            _soaNeighborList[firstCellIndex][secondCellIndex][particleIndexCurrentCell].second.emplace_back
                (_particleToIndex.at(_aosNeighborList[firstCellIndex][secondCellIndex][particleIndexCurrentCell].second[neighbor]));
          }
          particleIndexCurrentCell++;
        }
      }
    }
  }

  template <class TFunctor>
  auto *loadSoA(TFunctor *f) {
    _soa.clear();
    size_t offset = 0;
    size_t index = 0;
    for (auto &cell : _internalLinkedCells->getCells()) {
      f->SoALoader(cell, _soa, offset);
      offset += cell.numParticles();
      index++;
    }
    return &_soa;
  }

  template <class TFunctor>
  void extractSoA(TFunctor *f) {
    size_t offset = 0;
    for (auto &cell : _internalLinkedCells->getCells()) {
      f->SoAExtractor(cell, _soa, offset);
      offset += cell.numParticles();
    }
  }

 private:
  /**
   * Creates and applies generator functor for the building of the neighbor list.
   * @param linkedCells Linked Cells object used to build the neighbor list.
   * @param useNewton3 Whether Newton 3 should be used for the neighbor list.
   * @param cutoff Cutoff radius.
   * @param skin Skin of the verlet list.
   * @param interactionLength Interaction length of the underlying linked cells object.
   * @param buildTraversalOption Traversal option necessary for generator functor.
   * */
  void applyBuildFunctor(LinkedCells<Particle> &linkedCells, bool useNewton3, double cutoff, double skin,
                         double interactionLength, const TraversalOption buildTraversalOption) {
    typename VerletListsCellsHelpers<Particle>::PairwiseVerletListGeneratorFunctor f(
        _aosNeighborList, _particleToCellMap, _globalToLocalIndex, cutoff + skin);

    // Generate the build traversal with the traversal selector and apply the build functor with it.
    TraversalSelector<FullParticleCell<Particle>> traversalSelector;
    // Argument "cluster size" does not matter here.
    TraversalSelectorInfo traversalSelectorInfo(linkedCells.getCellBlock().getCellsPerDimensionWithHalo(),
                                                interactionLength, linkedCells.getCellBlock().getCellLength(), 0);
    autopas::utils::withStaticBool(useNewton3, [&](auto n3) {
      auto buildTraversal = traversalSelector.template generateTraversal<decltype(f), DataLayoutOption::aos, n3>( ///beware soa
          buildTraversalOption, f, traversalSelectorInfo);
      linkedCells.iteratePairwise(buildTraversal.get());
    });
  }

  void applyBuildFunctorSoA(LinkedCells<Particle> &linkedCells, bool useNewton3, double cutoff, double skin,
                         double interactionLength, const TraversalOption buildTraversalOption) {

    typename VerletListsCellsHelpers<Particle>::PairwiseVerletListGeneratorFunctor f(
        _aosNeighborList, _particleToCellMap, _globalToLocalIndex, cutoff + skin);

    //this->loadSoA(f);

    // Generate the build traversal with the traversal selector and apply the build functor with it.
    TraversalSelector<FullParticleCell<Particle>> traversalSelector;
    // Argument "cluster size" does not matter here.
    TraversalSelectorInfo traversalSelectorInfo(linkedCells.getCellBlock().getCellsPerDimensionWithHalo(),
                                                interactionLength, linkedCells.getCellBlock().getCellLength(), 0);
    autopas::utils::withStaticBool(useNewton3, [&](auto n3) {
      auto buildTraversal = traversalSelector.template generateTraversal<decltype(f), DataLayoutOption::soa, n3>( ///beware soa
          buildTraversalOption, f, traversalSelectorInfo);
      linkedCells.iteratePairwise(buildTraversal.get());
    });

    //this->extractSoA(f);
  }

  /**
   * Internal neighbor list structure in AoS format - Verlet lists for each particle for each cell.
   */
  typename VerletListsCellsHelpers<Particle>::PairwiseNeighborListsType _aosNeighborList;

  /**
   * Mapping of each particle to its corresponding cell and id within this cell.
   */
  std::unordered_map<Particle *, std::pair<size_t, size_t>> _particleToCellMap;

  /**
   * For each cell1: a mapping of the "absolute" index of cell2 (in the base linked cells structure) to its "relative"
   * index (0 to 26) in cell1's neighbors.
   */
  std::vector<std::unordered_map<size_t, size_t>> _globalToLocalIndex;

  std::vector<std::vector<std::vector<SoAPairOfParticleAndList>>> _soaNeighborList;
  SoA<typename Particle::SoAArraysType> _soa;

  LinkedCells<Particle>* _internalLinkedCells;
};
}  // namespace autopas
