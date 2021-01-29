/**
 * @file VLCCellPairNeighborList.h
 * @author tirgendetwas
 * @date 03.11.20 od. 07.11.20
 */

#pragma once
#include "VLCCellPairGeneratorFunctor.h"
#include "VLCNeighborListInterface.h"
#include "autopas/utils/StaticBoolSelector.h"

namespace autopas {

/**
 * TraversalSelector is used for the construction of the list in the applyBuildFunctor method.
 * Forward declaration necessary to avoid circle of includes:
 * TraversalSelector includes all VLC traversals include VLCTraversalInterface includes VLCCellPairNeighborList
 */
template <class ParticleCell>
class TraversalSelector;
/**
 * Neighbor list to be used with VerletListsCells container.
 * Pairwise verlet lists iterates through each pair of neighboring cells
 * and generates a neighbor list for each particle from cell1, which consists of its (potential) partners from cell2.
 * @tparam Particle Type of particle to be used for this neighbor list.
 */
template <class Particle>
class VLCCellPairNeighborList : public VLCNeighborListInterface<Particle> {
 public:
  /**
   * Type of the data structure used to save the neighbor lists.
   */
  using listType = typename VerletListsCellsHelpers<Particle>::PairwiseNeighborListsType;

  /**
   * @copydoc VLCNeighborListInterface::getContainerType()
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::pairwiseVerletLists; }

  /**
   * @copydoc VLCNeighborListInterface::getNumberOfPartners()
   */
  const size_t getNumberOfPartners(const Particle *particle) const override {
    size_t listSize = 0;
    const auto &[firstCellIndex, particleInCellIndex] = _particleToCellMap.at(const_cast<Particle *>(particle));
    for (auto &cellPair : _aosNeighborList[firstCellIndex]) {
      listSize += cellPair[particleInCellIndex].second.size();
    }
    return listSize;
  }

  /**
   * Returns the neighbor list in AoS layout.
   * @return Neighbor list in AoS layout.
   */
  typename VerletListsCellsHelpers<Particle>::PairwiseNeighborListsType &getAoSNeighborList() {
    return _aosNeighborList;
  }

  /**
   * @copydoc VLCNeighborListInterface::buildAoSNeighborList()
   */
  void buildAoSNeighborList(LinkedCells<Particle> &linkedCells, bool useNewton3, double cutoff, double skin,
                            double interactionLength, const TraversalOption buildTraversalOption) override {
    _aosNeighborList.clear();
    auto &cells = linkedCells.getCells();
    auto cellsSize = cells.size();
    _aosNeighborList.resize(cellsSize);
    _globalToLocalIndex.resize(cellsSize);

    size_t result = 0;
    auto cellLength = linkedCells.getCellBlock().getCellLength();
    auto cellsPerDimension = linkedCells.getCellBlock().getCellsPerDimensionWithHalo();
    auto interactionLengthSquare = linkedCells.getInteractionLength() * linkedCells.getInteractionLength();
    std::array<long, 3> overlap;
    size_t neighborCells = 0;

    for (unsigned int d = 0; d < 3; d++) {
      overlap[d] = std::ceil(linkedCells.getInteractionLength() / cellLength[d]);
    }

    for (int x = -overlap[0]; x < overlap[0] + 1; x++) {
      for (int y = -overlap[1]; y < overlap[1] + 1; y++) {
        for (int z = -overlap[2]; z < overlap[2] + 1; z++) {
          std::array<double, 3> pos = {};
          pos[0] = std::max(0l, (std::abs(x) - 1l)) * cellLength[0];
          pos[1] = std::max(0l, (std::abs(y) - 1l)) * cellLength[1];
          pos[2] = std::max(0l, (std::abs(z) - 1l)) * cellLength[2];
          const double distSquare = utils::ArrayMath::dot(pos, pos);
          if (distSquare <= interactionLengthSquare) {
            neighborCells++;
          }
        }
      }
    }

    for (size_t firstCellIndex = 0; firstCellIndex < cellsSize; ++firstCellIndex) {
      _aosNeighborList[firstCellIndex].resize(neighborCells);
      size_t numParticlesFirstCell = cells[firstCellIndex].numParticles();
      for (auto &cellPair : _aosNeighborList[firstCellIndex]) {
        // reserve vector of neighbor lists for every particle in cell1
        cellPair.reserve(numParticlesFirstCell);
        size_t particleIndexCurrentCell = 0;
        for (auto particleIter = cells[firstCellIndex].begin(); particleIter.isValid(); ++particleIter) {
          // for each particle in cell1 make pair of particle and neighbor list
          Particle *currentParticle = &*particleIter;
          cellPair.emplace_back(std::make_pair(currentParticle, std::vector<Particle *>()));
          // TODO reserve experiment - Tina

          // add pair of cell's index and particle's index in the cell
          _particleToCellMap[currentParticle] = std::make_pair(firstCellIndex, particleIndexCurrentCell);
          particleIndexCurrentCell++;
        }
      }
    }

    applyBuildFunctor(linkedCells, useNewton3, cutoff, skin, interactionLength, buildTraversalOption);
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
   */
  void applyBuildFunctor(LinkedCells<Particle> &linkedCells, bool useNewton3, double cutoff, double skin,
                         double interactionLength, const TraversalOption buildTraversalOption) {
    VLCCellPairGeneratorFunctor<Particle> f(_aosNeighborList, _particleToCellMap, _globalToLocalIndex, cutoff + skin);

    // Generate the build traversal with the traversal selector and apply the build functor with it.
    TraversalSelector<FullParticleCell<Particle>> traversalSelector;
    // Argument "cluster size" does not matter here.
    TraversalSelectorInfo traversalSelectorInfo(linkedCells.getCellBlock().getCellsPerDimensionWithHalo(),
                                                interactionLength, linkedCells.getCellBlock().getCellLength(), 0);
    autopas::utils::withStaticBool(useNewton3, [&](auto n3) {
      auto buildTraversal = traversalSelector.template generateTraversal<std::remove_reference_t<decltype(f)>, DataLayoutOption::aos, n3>(
          buildTraversalOption, f, traversalSelectorInfo);
      linkedCells.iteratePairwise(buildTraversal.get());
    });
  }

  /**
   * Internal neighbor list structure in AoS format - Verlet lists for each particle for each pair of cells.
   */
  typename VerletListsCellsHelpers<Particle>::PairwiseNeighborListsType _aosNeighborList =
      std::vector<std::vector<std::vector<std::pair<Particle *, std::vector<Particle *>>>>>();

  /**
   * Mapping of each particle to its corresponding cell and id within this cell.
   */
  std::unordered_map<Particle *, std::pair<size_t, size_t>> _particleToCellMap =
      std::unordered_map<Particle *, std::pair<size_t, size_t>>();

  /**
   * For each cell1: a mapping of the "absolute" index of cell2 (in the base linked cells structure) to its "relative"
   * index in cell1's neighbors.
   */
  std::vector<std::unordered_map<size_t, size_t>> _globalToLocalIndex =
      std::vector<std::unordered_map<size_t, size_t>>();
};
}  // namespace autopas
