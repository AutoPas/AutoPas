/**
 * @file VLCCellPairNeighborList.h
 * @author tirgendetwas
 * @date 07.11.20
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

template <class Particle>
class VLCCellPairTraversalInterface;
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
   * Helper type definition. Pair of particle and neighbor list for SoA layout.
   */
  using SoAPairOfParticleAndList = std::pair<size_t, std::vector<size_t, autopas::AlignedAllocator<size_t>>>;

  /**
   * Helper type definition. Vector of cells, for each cell a vector of neighbors.
   * For each pair of cells, a vector of mappings from particle to its neighbor list.
   */
  using SoAListType = typename std::vector<std::vector<std::vector<SoAPairOfParticleAndList>>>;

  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::pairwiseVerletLists; }

  size_t getNumberOfPartners(const Particle *particle) const override {
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
   * Returns a map of each cell's global index to its local index in another cell's neighbor list. More specifically,
   * for each cell1: a mapping of the "absolute" index of cell2 (in the base linked cells structure) to its "relative"
   * index in cell1's neighbors.
   * @return a map of each cell's global index to its local index in another cell's neighbor list
   */
  auto &getGlobalToLocalMap() { return _globalToLocalIndex; }

  /**
   * Returns the neighbor list in SoA layout.
   * @return Neighbor list in SoA layout.
   */
  auto &getSoANeighborList() { return _soaNeighborList; }

  void buildAoSNeighborList(LinkedCells<Particle> &linkedCells, bool useNewton3, double cutoff, double skin,
                            double interactionLength, const TraversalOption buildTraversalOption,
                            typename VerletListsCellsHelpers<Particle>::VLCBuildType::Value buildType) override {
    this->_internalLinkedCells = &linkedCells;
    _aosNeighborList.clear();
    _globalToLocalIndex.clear();
    _particleToCellMap.clear();
    auto &cells = linkedCells.getCells();
    const auto cellsSize = cells.size();
    _aosNeighborList.resize(cellsSize);
    _globalToLocalIndex.resize(cellsSize);

    const auto cellLength = linkedCells.getCellBlock().getCellLength();
    const auto interactionLengthSquare = linkedCells.getInteractionLength() * linkedCells.getInteractionLength();

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

    // initialize empty lists for every particle-cell pair
    for (size_t firstCellIndex = 0; firstCellIndex < cellsSize; ++firstCellIndex) {
      _aosNeighborList[firstCellIndex].resize(neighborCells);
      size_t numParticlesFirstCell = cells[firstCellIndex].numParticles();
      for (auto &cellPair : _aosNeighborList[firstCellIndex]) {
        // reserve vector of neighbor lists for every particle in cell1
        cellPair.reserve(numParticlesFirstCell);
        size_t particleIndexCurrentCell = 0;
        for (auto &particle : cells[firstCellIndex]) {
          // for each particle in cell1 make a pair of particle and neighbor list
          cellPair.emplace_back(std::make_pair(&particle, std::vector<Particle *>()));

          // add a pair of cell's index and particle's index in the cell
          _particleToCellMap[&particle] = std::make_pair(firstCellIndex, particleIndexCurrentCell);
          particleIndexCurrentCell++;
        }
      }
    }

    // fill the lists
    applyBuildFunctor(linkedCells, useNewton3, cutoff, skin, interactionLength, buildTraversalOption, buildType);
  }

  void generateSoAFromAoS(LinkedCells<Particle> &linkedCells) override {
    _soaNeighborList.clear();

    // particle pointer to global index of particle
    std::unordered_map<Particle *, size_t> particleToIndex;
    particleToIndex.reserve(linkedCells.getNumberOfParticles());
    size_t i = 0;
    for (auto iter = linkedCells.begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter, ++i) {
      particleToIndex[&(*iter)] = i;
    }

    _soaNeighborList.resize(linkedCells.getCells().size());

    // iterate over cells
    for (size_t firstCellIndex = 0; firstCellIndex < _aosNeighborList.size(); ++firstCellIndex) {
      const auto &aosCurrentCell = _aosNeighborList[firstCellIndex];
      auto &soaCurrentCell = _soaNeighborList[firstCellIndex];
      soaCurrentCell.resize(aosCurrentCell.size());

      // iterate over each cell's neighboring cells
      for (size_t secondCellIndex = 0; secondCellIndex < aosCurrentCell.size(); ++secondCellIndex) {
        const auto &aosCellPair = aosCurrentCell[secondCellIndex];
        auto &soaCellPair = soaCurrentCell[secondCellIndex];
        soaCellPair.reserve(aosCellPair.capacity());

        // iterate over pairs of particle and neighbor list
        for (const auto &aosParticleAndParticleList : aosCellPair) {
          Particle *currentParticle = aosParticleAndParticleList.first;
          // global index of current particle
          size_t currentParticleGlobalIndex = particleToIndex.at(currentParticle);

          // create SoA neighbor list for current particle
          std::vector<size_t, autopas::AlignedAllocator<size_t>> currentSoANeighborList{};

          // fill the SoA neighbor list with the indices of the particles from the corresponding AoS neighbor list
          for (const auto &neighborOfCurrentParticle : aosParticleAndParticleList.second) {
            currentSoANeighborList.emplace_back(particleToIndex.at(neighborOfCurrentParticle));
          }

          // add newly constructed pair of particle index and SoA particle neighbor list to cell pair
          soaCellPair.emplace_back(std::make_pair(currentParticleGlobalIndex, currentSoANeighborList));
        }
      }
    }
  }

  void setUpTraversal(TraversalInterface *traversal) override {
    auto vTraversal = dynamic_cast<VLCCellPairTraversalInterface<Particle> *>(traversal);

    if (vTraversal) {
      vTraversal->setVerletList(*this);
    } else {
      auto traversal2 = dynamic_cast<VLCTraversalInterface<Particle, VLCCellPairNeighborList<Particle>> *>(traversal);
      if (traversal2) {
        traversal2->setVerletList(*this);
      } else {
        autopas::utils::ExceptionHandler::exception(
            "Trying to use a traversal of wrong type in VerletListCells.h. TraversalID: {}",
            traversal->getTraversalType());
      }
    }
  }

 private:
  void applyBuildFunctor(LinkedCells<Particle> &linkedCells, bool useNewton3, double cutoff, double skin,
                         double interactionLength, const TraversalOption buildTraversalOption,
                         typename VerletListsCellsHelpers<Particle>::VLCBuildType::Value buildType) override {
    VLCCellPairGeneratorFunctor<Particle> f(_aosNeighborList, _particleToCellMap, _globalToLocalIndex, cutoff + skin);

    // Generate the build traversal with the traversal selector and apply the build functor with it.
    TraversalSelector<FullParticleCell<Particle>> traversalSelector;
    // Argument "cluster size" does not matter here.
    TraversalSelectorInfo traversalSelectorInfo(linkedCells.getCellBlock().getCellsPerDimensionWithHalo(),
                                                interactionLength, linkedCells.getCellBlock().getCellLength(), 0);

    // Build the AoS list using the AoS or SoA functor depending on buildType
    if (buildType == VerletListsCellsHelpers<Particle>::VLCBuildType::Value::aosBuild) {
      autopas::utils::withStaticBool(useNewton3, [&](auto n3) {
        auto buildTraversal =
            traversalSelector
                .template generateTraversal<std::remove_reference_t<decltype(f)>, DataLayoutOption::aos, n3>(
                    buildTraversalOption, f, traversalSelectorInfo);
        linkedCells.iteratePairwise(buildTraversal.get());
      });
    }

    else if (buildType == VerletListsCellsHelpers<Particle>::VLCBuildType::Value::soaBuild) {
      autopas::utils::withStaticBool(useNewton3, [&](auto n3) {
        auto buildTraversal = traversalSelector.template generateTraversal<decltype(f), DataLayoutOption::soa, n3>(
            buildTraversalOption, f, traversalSelectorInfo);
        linkedCells.iteratePairwise(buildTraversal.get());
      });
    }
  }

  /**
   * Internal neighbor list structure in AoS format - Verlet lists for each particle for each cell pair.
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

  /**
   * Internal neighbor list structure in SoA format - Verlet lists for each particle for each cell pair.
   * Contrary to aosNeighborList it saves global particle indices instead of particle pointers.
   */
  SoAListType _soaNeighborList = std::vector<
      std::vector<std::vector<std::pair<size_t, std::vector<size_t, autopas::AlignedAllocator<size_t>>>>>>();
};
}  // namespace autopas
