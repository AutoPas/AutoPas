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
  using listType = typename VerletListsCellsHelpers::AllCellsNeighborListsType<Particle>;

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
                            double interactionLength, const TraversalOption vlcTraversalOpt,
                            typename VerletListsCellsHelpers::VLCBuildType buildType) override {
    using namespace utils::ArrayMath::literals;
    this->_internalLinkedCells = &linkedCells;
    auto &cells = linkedCells.getCells();
    const auto cellsPerDim = this->_internalLinkedCells->getCellBlock().getCellsPerDimensionWithHalo();

    const auto boxSizeWithHalo = this->_internalLinkedCells->getBoxMax() - this->_internalLinkedCells->getBoxMin() +
                                 std::array<double, 3>{interactionLength, interactionLength, interactionLength} * 2.;
    const auto listLengthEstimate = VerletListsCellsHelpers::estimateListLength(
        this->_internalLinkedCells->getNumberOfParticles(IteratorBehavior::ownedOrHalo), boxSizeWithHalo,
        interactionLength, 1.3);

    const auto offsetsC08 = VerletListsCellsHelpers::buildBaseStep(
        utils::ArrayUtils::static_cast_copy_array<int>(
            this->_internalLinkedCells->getCellBlock().getCellsPerDimensionWithHalo()),
        TraversalOption::vlc_c08);

    // Helper function to estimate the number of neighbor lists for one base step
    // TODO: This is a generous and rough estimate and can probably be improved!
    const auto estimateNumLists = [&](size_t baseCellIndex) {
      // If the cell is near the end of any dimension we can reduce the estimate to the particles in the base cell,
      // because all other cells' interactions would be from halos.
      const auto cellIndex3D = utils::ThreeDimensionalMapping::oneToThreeD(baseCellIndex, cellsPerDim);
      const auto numTouchedBoundaries = [&]() {
        int acc = 0;
        for (size_t i = 0; i < 3; ++i) {
          // false == 0 ; true == 1
          acc += static_cast<int>(cellIndex3D[i] == (cellsPerDim[i] - 1));
        }
        return acc;
      }();
      if (numTouchedBoundaries > 0) {
        if (useNewton3) {
          return cells[baseCellIndex].size();
        } else {
          // In this case we have to accommodate the lists for the reverse interactions from all non-halo neighbors.
          // 1x for lists from the cell itself and once more per untouched border.
          return cells[baseCellIndex].size() * (4 - numTouchedBoundaries);
        }
      }
      size_t estimate = 0;
      std::map<int, double> offsetFactors{};
      std::map<int, double> offsetFactorsNoN3{};
      for (const auto [offsetA, offsetB, factor] : offsetsC08) {
        offsetFactors[offsetA] = std::max(offsetFactors[offsetA], factor);
        offsetFactorsNoN3[offsetB] = std::max(offsetFactors[offsetB], factor);
      }
      for (const auto &[offset, factor] : offsetFactors) {
        estimate += cells[baseCellIndex + offset].size() * factor;
      }
      if (not useNewton3) {
        for (const auto &[offset, factor] : offsetFactorsNoN3) {
          estimate += cells[baseCellIndex + offset].size() * factor;
        }
      }
      return estimate;
    };

    // Initialize a map of neighbor lists for each cell.
    _aosNeighborList.clear();
    const size_t numCells = cells.size();
    _aosNeighborList.resize(numCells);
    for (size_t cellIndex = 0; cellIndex < numCells; ++cellIndex) {
      const auto estimateForNumberOfLists = [&]() {
        // Usually each cell holds one list per particle, except vlc_c08, which holds all lists of a base step,
        // which involves also lists from other cells' particles.
        if (vlcTraversalOpt == TraversalOption::vlc_c08) {
          return estimateNumLists(cellIndex);
        } else {
          return cells[cellIndex].size();
        }
      }();
      auto &cell = cells[cellIndex];
      _aosNeighborList[cellIndex].reserve(estimateForNumberOfLists);
      size_t particleIndexWithinCell = 0;
      for (auto iter = cell.begin(); iter != cell.end(); ++iter, ++particleIndexWithinCell) {
        Particle *particle = &*iter;
        _aosNeighborList[cellIndex].emplace_back(particle, std::vector<Particle *>());
        _aosNeighborList[cellIndex].back().second.reserve(listLengthEstimate);
        _particleToCellMap[particle] = std::make_pair(cellIndex, particleIndexWithinCell);
      }
    }

    applyBuildFunctor(linkedCells, useNewton3, cutoff, skin, interactionLength, vlcTraversalOpt, buildType);
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
  typename VerletListsCellsHelpers::AllCellsNeighborListsType<Particle> &getAoSNeighborList() {
    return _aosNeighborList;
  }

  /**
   * Returns a Mapping of particles to its corresponding cell and index within this cell.
   * @return Mapping of particles to its corresponding cell and index within this cell.
   */
  auto &getParticleToCellMap() { return _particleToCellMap; }

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
    std::unordered_map<Particle *, size_t> particlePtrToIndex;
    particlePtrToIndex.reserve(linkedCells.size());
    size_t i = 0;
    for (auto iter = linkedCells.begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter, ++i) {
      particlePtrToIndex[&(*iter)] = i;
    }

    _soaNeighborList.resize(linkedCells.getCells().size());

    // iterate over cells and for each create the soa lists from the aos lists
    for (size_t firstCellIndex = 0; firstCellIndex < _aosNeighborList.size(); ++firstCellIndex) {
      const auto &aosLists = _aosNeighborList[firstCellIndex];
      auto &soaLists = _soaNeighborList[firstCellIndex];
      soaLists.reserve(aosLists.size());

      // iterate over pairs of particle and neighbor list
      for (const auto &[particlePtr, neighbors] : aosLists) {
        // global index of current particle
        size_t currentParticleGlobalIndex = particlePtrToIndex.at(particlePtr);

        // create SoA neighbor list for current particle
        std::vector<size_t, autopas::AlignedAllocator<size_t>> currentSoANeighborList{};
        currentSoANeighborList.reserve(neighbors.size());

        // fill the SoA neighbor list with the indices of the particles from the corresponding AoS neighbor list
        for (const auto &neighborOfCurrentParticle : neighbors) {
          currentSoANeighborList.emplace_back(particlePtrToIndex.at(neighborOfCurrentParticle));
        }

        // add the newly constructed pair of particle index and SoA neighbor list to cell
        soaLists.emplace_back(currentParticleGlobalIndex, currentSoANeighborList);
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
                         double interactionLength, const TraversalOption &vlcTraversalOpt,
                         typename VerletListsCellsHelpers::VLCBuildType buildType) override {
    // Generate the build traversal with the traversal selector and apply the build functor with it.
    TraversalSelector<FullParticleCell<Particle>> traversalSelector;
    // Argument "cluster size" does not matter here.
    TraversalSelectorInfo traversalSelectorInfo(linkedCells.getCellBlock().getCellsPerDimensionWithHalo(),
                                                interactionLength, linkedCells.getCellBlock().getCellLength(), 0);

    const auto dataLayout =
        buildType == VerletListsCellsHelpers::VLCBuildType::aosBuild ? DataLayoutOption::aos : DataLayoutOption::soa;

    // LC traversal that will be used to build the List
    const TraversalOption lcBuildTraversalOpt =
        vlcTraversalOpt == TraversalOption::vlc_c08 ? TraversalOption::lc_c08 : TraversalOption::lc_c18;

    using namespace utils::ArrayMath::literals;
    const auto boxSizeWithHalo = linkedCells.getBoxMax() - linkedCells.getBoxMin() +
                                 std::array<double, 3>{interactionLength, interactionLength, interactionLength} * 2.;
    const auto listLengthEstimate = VerletListsCellsHelpers::estimateListLength(
        linkedCells.getNumberOfParticles(IteratorBehavior::ownedOrHalo), boxSizeWithHalo, interactionLength, 1.3);

    if (vlcTraversalOpt == TraversalOption::vlc_c08) {
      VLCAllCellsGeneratorFunctor<Particle, TraversalOption::vlc_c08> f(
          _aosNeighborList, _particleToCellMap, cutoff + skin, listLengthEstimate,
          linkedCells.getCellBlock().getCellsPerDimensionWithHalo());
      f.setCells(&this->_internalLinkedCells->getCells());
      // Build the AoS list using the AoS or SoA functor depending on buildType
      auto buildTraversal = traversalSelector.template generateTraversal<std::remove_reference_t<decltype(f)>>(
          lcBuildTraversalOpt, f, traversalSelectorInfo, dataLayout, useNewton3);
      linkedCells.iteratePairwise(buildTraversal.get());
    } else {
      VLCAllCellsGeneratorFunctor<Particle, TraversalOption::vlc_c18> f(
          _aosNeighborList, _particleToCellMap, cutoff + skin, listLengthEstimate,
          linkedCells.getCellBlock().getCellsPerDimensionWithHalo());
      // Build the AoS list using the AoS or SoA functor depending on buildType
      auto buildTraversal = traversalSelector.template generateTraversal<std::remove_reference_t<decltype(f)>>(
          lcBuildTraversalOpt, f, traversalSelectorInfo, dataLayout, useNewton3);
      linkedCells.iteratePairwise(buildTraversal.get());
    }
  }

  /**
   * Internal neighbor list structure in AoS format - Verlet lists for each particle for each cell.
   */
  typename VerletListsCellsHelpers::AllCellsNeighborListsType<Particle> _aosNeighborList{};

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
