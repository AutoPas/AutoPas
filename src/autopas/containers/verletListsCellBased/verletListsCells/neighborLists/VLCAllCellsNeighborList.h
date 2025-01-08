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
   * Special case of building the neighbor lists in c08 and c18 style where all lists that belong to one base step are
   * stored together.
   * @param traversal The TraversalOption with which the function is called.
   * @param linkedCells A reference to the linked cells container.
   * @param verletBuiltNewton3 Boolean to specify if newton3 is enabled or not.
   */
  void rebuildNeighborListsVLC(TraversalOption traversal, LinkedCells<Particle> &linkedCells, bool verletBuiltNewton3) {
    using namespace utils::ArrayMath::literals;
    // Sanity check.
    if (linkedCells.getCellBlock().getCellsPerInteractionLength() > 1) {
      utils::ExceptionHandler::exception(
          "VerletListsCells::rebuildNeighborListsVLC() was called with a CSF < 1 but it only supports CSF>=1.");
    }
    // Define some aliases
    auto &neighborLists = getAoSNeighborList();
    auto &cells = linkedCells.getCells();
    const auto interactionLength = linkedCells.getInteractionLength();
    const auto interactionLengthSquared = interactionLength * interactionLength;
    const auto boxSizeWithHalo = linkedCells.getBoxMax() - linkedCells.getBoxMin() +
                                 std::array<double, 3>{interactionLength, interactionLength, interactionLength} * 2.;
    // Create an estimate for the average length of a neighbor list.
    // This assumes homogeneous distribution and some overestimation.
    const auto listLengthEstimate = VerletListsCellsHelpers::estimateListLength(
        linkedCells.getNumberOfParticles(IteratorBehavior::ownedOrHalo), boxSizeWithHalo, interactionLength, 1.3);

    // Reset lists. Don't free any memory, only mark as unused.
    this->setLinkedCellsPointer(&linkedCells);
    for (auto &cellLists : neighborLists) {
      for (auto &[particlePtr, neighbors] : cellLists) {
        particlePtr = nullptr;
        neighbors.clear();
      }
    }
    neighborLists.resize(cells.size());

    /* This must not be a doc comment (with two **) to not confuse doxygen.
     * Helper function to insert a pointer into a list of the base cell.
     * It considers the cases that neither particle is in the base cell
     * and in that case finds or creates the appropriate list.
     *
     * @param p1 Reference to source particle.
     * @param p1Index Index of p1 in its cell.
     * @param p2 Reference to target particle.
     * @param cellIndex1 Index of cell where p1 resides.
     * @param cellIndexBase Index of the base cell of the c08 step.
     * @param neighborList Reference to the list where the particle pair should be stored.
     */
    auto insert = [&](auto &p1, auto p1Index, auto &p2, auto cellIndex1, auto cellIndexBase, auto &neighborList) {
      // Easy case: cell1 is the base cell
      if (cellIndexBase == cellIndex1) {
        neighborList[p1Index].second.push_back(&p2);
      } else {
        // Otherwise, check if the base cell already has a list for p1
        auto iter = std::find_if(neighborList.begin(), neighborList.end(), [&](const auto &pair) {
          const auto &[particlePtr, list] = pair;
          return particlePtr == &p1;
        });
        // If yes, append p2 to it.
        if (iter != neighborList.end()) {
          iter->second.push_back(&p2);
        } else {
          // If no, create one (or reuse an empty pair), reserve space for the list and emplace p2
          if (auto insertLoc = std::find_if(neighborList.begin(), neighborList.end(),
                                            [&](const auto &pair) {
                                              const auto &[particlePtr, list] = pair;
                                              return particlePtr == nullptr;
                                            });
              insertLoc != neighborList.end()) {
            auto &[particlePtr, neighbors] = *insertLoc;
            particlePtr = &p1;
            neighbors.reserve(listLengthEstimate);
            neighbors.push_back(&p2);
          } else {
            neighborList.emplace_back(&p1, std::vector<Particle *>{});
            neighborList.back().second.reserve(listLengthEstimate);
            neighborList.back().second.push_back(&p2);
          }
        }
      }
    };

    const auto &cellsPerDim =
        utils::ArrayUtils::static_cast_copy_array<int>(linkedCells.getCellBlock().getCellsPerDimensionWithHalo());
    // Vector of offsets from the base cell for the c08 base step
    // and respective factors for the fraction of particles per cell that need neighbor lists in the base cell.
    const auto offsets = VerletListsCellsHelpers::buildBaseStep(cellsPerDim, traversal);

    int xEnd{};
    int yEnd{};
    int zEnd{};

    switch (traversal) {
      case TraversalOption::vlc_c08:
        // Go over all cells except the very last layer and create lists per base step.
        xEnd = cellsPerDim[0] - 1;
        yEnd = cellsPerDim[1] - 1;
        zEnd = cellsPerDim[2] - 1;
        break;
      default:
        xEnd = cellsPerDim[0];
        yEnd = cellsPerDim[1];
        zEnd = cellsPerDim[2];
        break;
    }

    // Since there are no loop dependencies merge all for loops and create 10 chunks per thread.
    AUTOPAS_OPENMP(parallel for collapse(3) schedule(dynamic, std::max(cells.size() / (autopas::autopas_get_max_threads() * 10), 1ul)))
    for (int z = 0; z < zEnd; ++z) {
      for (int y = 0; y < yEnd; ++y) {
        for (int x = 0; x < xEnd; ++x) {
          // aliases
          const auto cellIndexBase = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, cellsPerDim);
          auto &baseCell = cells[cellIndexBase];
          auto &baseCellsLists = neighborLists[cellIndexBase];

          // Allocate memory for ptr-list pairs for this cell.
          baseCellsLists.resize(VerletListsCellsHelpers::estimateNumLists(
              cellIndexBase, verletBuiltNewton3, cells, offsets,
              utils::ArrayUtils::static_cast_copy_array<size_t>(cellsPerDim)));
          // Re-initialize a neighbor list for all particles in the cell but at most for as many as there are lists
          const size_t minCellSizeVsNumLists = std::min(baseCell.size(), baseCellsLists.size());
          for (size_t i = 0; i < minCellSizeVsNumLists; ++i) {
            auto &[particlePtr, neighbors] = baseCellsLists[i];
            particlePtr = &baseCell[i];
            neighbors.reserve(listLengthEstimate);
          }
          // For any remaining particles create a new list.
          // This case can only happen if estimateNumLists can return values smaller than baseCell.size()
          for (size_t i = minCellSizeVsNumLists; i < baseCell.size(); ++i) {
            baseCellsLists.emplace_back(&baseCell[i], std::vector<Particle *>{});
            baseCellsLists.back().second.reserve(listLengthEstimate);
          }

          // Build c08 lists for this base step according to predefined cell pairs
          for (const auto &[offset1, offset2, _] : offsets) {
            const auto cellIndex1 = cellIndexBase + offset1;
            const auto cellIndex2 = cellIndexBase + offset2;

            // For all traversals ensures the partner cell is not outside the boundary
            const auto cell2Coords = utils::ThreeDimensionalMapping::oneToThreeD(cellIndex2, cellsPerDim);
            if (cell2Coords[0] >= cellsPerDim[0] or cell2Coords[0] < 0 or cell2Coords[1] >= cellsPerDim[1] or
                cell2Coords[1] < 0 or cell2Coords[2] >= cellsPerDim[2] or cell2Coords[2] < 0) {
              continue;
            }

            // Skip if both cells only contain halos
            if (not(cells[cellIndex1].getPossibleParticleOwnerships() == OwnershipState::owned) and
                not(cells[cellIndex2].getPossibleParticleOwnerships() == OwnershipState::owned)) {
              continue;
            }

            // Go over all particle pairs in the two cells and insert close pairs into their respective lists
            for (size_t particleIndexCell1 = 0; particleIndexCell1 < cells[cellIndex1].size(); ++particleIndexCell1) {
              auto &p1 = cells[cellIndex1][particleIndexCell1];

              // Determine starting index for the second particle in the second cell
              // If both cells are the same start after the current particle in the first cell.
              // For vlc_c01 we have to iterate over all pairs like p1<->p2 and p2<->p1 because we do not make use of
              // newton3
              size_t startIndexCell2 = 0;
              if (cellIndex1 == cellIndex2 and traversal != TraversalOption::vlc_c01) {
                startIndexCell2 = particleIndexCell1 + 1;
              }

              for (size_t particleIndexCell2 = startIndexCell2; particleIndexCell2 < cells[cellIndex2].size();
                   ++particleIndexCell2) {
                auto &p2 = cells[cellIndex2][particleIndexCell2];
                // Ignore dummies and self interaction
                if (&p1 == &p2 or p1.isDummy() or p2.isDummy()) {
                  continue;
                }

                // If the distance is less than interaction length add the pair to the list
                const auto distVec = p2.getR() - p1.getR();
                const auto distSquared = utils::ArrayMath::dot(distVec, distVec);
                if (distSquared < interactionLengthSquared) {
                  insert(p1, particleIndexCell1, p2, cellIndex1, cellIndexBase, baseCellsLists);
                  // If the traversal does not use Newton3 the inverse interaction also needs to be stored in p2's list
                  if (not verletBuiltNewton3 and not(traversal == TraversalOption::vlc_c01)) {
                    insert(p2, particleIndexCell2, p1, cellIndex2, cellIndexBase, baseCellsLists);
                  }
                }
              }
            }
          }
        }
      }
    }

    // Cleanup: Remove any unused ptr-list pairs to avoid accessing nullptr
    for (auto &cellLists : neighborLists) {
      cellLists.erase(std::remove_if(cellLists.begin(), cellLists.end(),
                                     [](const auto &pair) {
                                       const auto &[particlePtr, neighbors] = pair;
                                       return particlePtr == nullptr;
                                     }),
                      cellLists.end());
    }
  }

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
    for (const auto &cellsLists : _aosNeighborList) {
      for (const auto &particlesLists : cellsLists) {
        if (particlesLists.first == particle) {
          return particlesLists.second.size();
        }
      }
    }
    return 0lu;
  }

  /**
   * Returns the neighbor list in AoS layout.
   * @return Neighbor list in AoS layout.
   */
  typename VerletListsCellsHelpers::AllCellsNeighborListsType<Particle> &getAoSNeighborList() {
    return _aosNeighborList;
  }

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
      linkedCells.computeInteractions(buildTraversal.get());
    } else {
      VLCAllCellsGeneratorFunctor<Particle, TraversalOption::vlc_c18> f(
          _aosNeighborList, _particleToCellMap, cutoff + skin, listLengthEstimate,
          linkedCells.getCellBlock().getCellsPerDimensionWithHalo());
      // Build the AoS list using the AoS or SoA functor depending on buildType
      auto buildTraversal = traversalSelector.template generateTraversal<std::remove_reference_t<decltype(f)>>(
          lcBuildTraversalOpt, f, traversalSelectorInfo, dataLayout, useNewton3);
      linkedCells.computeInteractions(buildTraversal.get());
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
