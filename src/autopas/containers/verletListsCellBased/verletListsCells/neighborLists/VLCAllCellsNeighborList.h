/**
 * @file VLCAllCellsNeighborList.h
 * @author tirgendetwas
 * @date 25.10.20
 *
 * originally from
 * @file VerletListsCells.h
 * @author nguyen
 * @date 30.08.18
 */

#pragma once

#include "VLCAllCellsGeneratorFunctor.h"
#include "VLCNeighborListInterface.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StaticBoolSelector.h"

namespace autopas {

/**
 * TraversalSelector is used for the construction of the list in the applyBuildFunctor method.
 * Forward declaration necessary to avoid circle of includes:
 * TraversalSelector includes all VLC traversals include VLCTraversalInterface includes VLCAllCellsNeighborList
 */
template <class ParticleCell>
class TraversalSelector;

template <class Particle, class NeighborList>
class VLCTraversalInterface;

/**
 * Neighbor list to be used with VerletListsCells container. Classic implementation of verlet lists based on linked
 * cells.
 * @tparam Particle Type of particle to be used for this neighbor list.
 */
template <class Particle>
class VLCAllCellsNeighborList : public VLCNeighborListInterface<Particle> {
 public:
  /**
   * Type of the data structure used to save the neighbor lists.
   */
  using listType = typename VerletListsCellsHelpers<Particle>::NeighborListsType;

  /**
   * Helper type definition. Pair of particle and neighbor list for SoA layout.
   */
  using SoAPairOfParticleAndList = std::pair<size_t, std::vector<size_t, autopas::AlignedAllocator<size_t>>>;

  /**
   * Helper type definition. Vector of cells, for each cell a vector of neighbors.
   * For each pair of cells, a vector of mappings from particle to its neighbor list.
   */
  using SoAListType = typename std::vector<std::vector<SoAPairOfParticleAndList>>;

  /**
   * @copydoc VLCNeighborListInterface::getContainerType()
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::verletListsCells; }

  /**
   * @copydoc VLCNeighborListInterface::buildAoSNeighborList()
   */
  void buildAoSNeighborList(LinkedCells<Particle> &linkedCells, bool useNewton3, double cutoff, double skin,
                            double interactionLength, const TraversalOption buildTraversalOption,
                            typename VerletListsCellsHelpers<Particle>::VLCBuildType::Value buildType) override {
    // Initialize a neighbor list for each cell.
    _aosNeighborList.clear();
    this->_internalLinkedCells = &linkedCells;
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

    applyBuildFunctor(linkedCells, useNewton3, cutoff, skin, interactionLength, buildTraversalOption, buildType);
  }

  /**
   * @copydoc VLCNeighborListInterface::getNumberOfPartners()
   */
  size_t getNumberOfPartners(const Particle *particle) const override {
    const auto &[cellIndex, particleIndexInCell] = _particleToCellMap.at(const_cast<Particle *>(particle));
    return _aosNeighborList.at(cellIndex).at(particleIndexInCell).second.size();
  }

  /**
   * Returns the neighbor list in AoS layout.
   * @return Neighbor list in AoS layout.
   */
  typename VerletListsCellsHelpers<Particle>::NeighborListsType &getAoSNeighborList() { return _aosNeighborList; }

  /**
   * Returns the neighbor list in SoA layout.
   * @return Neighbor list in SoA layout.
   */
  auto &getSoANeighborList() { return _soaNeighborList; }

  /**
   * @copydoc VLCNeighborListInterface::generateSoAFromAoS()
   */
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
      soaCurrentCell.reserve(aosCurrentCell.capacity());

      // iterate over pairs of particle and neighbor list
      for (const auto &aosParticleAndParticleList : aosCurrentCell) {
        Particle *currentParticle = aosParticleAndParticleList.first;
        // global index of current particle
        size_t currentParticleGlobalIndex = particleToIndex.at(currentParticle);

        // create SoA neighbor list for current particle
        std::vector<size_t, autopas::AlignedAllocator<size_t>> currentSoANeighborList{};

        // fill the SoA neighbor list with the indices of the particles from the corresponding AoS neighbor list
        for (const auto &neighborOfCurrentParticle : aosParticleAndParticleList.second) {
          currentSoANeighborList.emplace_back(particleToIndex.at(neighborOfCurrentParticle));
        }

        // add newly constructed pair of particle index and SoA neighbor list to cell
        soaCurrentCell.emplace_back(std::make_pair(currentParticleGlobalIndex, currentSoANeighborList));
      }
    }
  }

  void setUpTraversal(TraversalInterface *traversal) override {
    auto vTraversal = dynamic_cast<VLCTraversalInterface<Particle, VLCAllCellsNeighborList<Particle>> *>(traversal);

    if (vTraversal) {
      vTraversal->setVerletList(*this);
    } else {
      autopas::utils::ExceptionHandler::exception(
          "Trying to use a traversal of wrong type in VerletListCells.h. TraversalID: {}",
          traversal->getTraversalType());
    }
  }

 private:
  /**
   * @copydoc VLCNeighborListInterface::applyBuildFunctor()
   */
  void applyBuildFunctor(LinkedCells<Particle> &linkedCells, bool useNewton3, double cutoff, double skin,
                         double interactionLength, const TraversalOption buildTraversalOption,
                         typename VerletListsCellsHelpers<Particle>::VLCBuildType::Value buildType) override {
    VLCAllCellsGeneratorFunctor<Particle> f(_aosNeighborList, _particleToCellMap, cutoff + skin);

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
   * Internal neighbor list structure in AoS format - Verlet lists for each particle for each cell.
   */
  typename VerletListsCellsHelpers<Particle>::NeighborListsType _aosNeighborList{};

  /**
   * Mapping of each particle to its corresponding cell and id within this cell.
   */
  std::unordered_map<Particle *, std::pair<size_t, size_t>> _particleToCellMap{};

  /**
   * Internal neighbor list structure in SoA format - Verlet lists for each particle for each cell.
   * Contrary to aosNeighborList it saves global particle indices instead of particle pointers.
   */
  SoAListType _soaNeighborList{};
};
}  // namespace autopas
