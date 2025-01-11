/**
 * @file VLCCellPairNeighborList.h
 * @author tirgendetwas
 * @date 07.11.20
 */

#pragma once
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
  using listType = typename VerletListsCellsHelpers::PairwiseNeighborListsType<Particle>;

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
    for (auto &neighborCellsLists : _aosNeighborList) {
      for (const auto &cellsLists : neighborCellsLists) {
        for (size_t i{0}; i < cellsLists.size(); ++i) {
          if (cellsLists[i].first == particle) {
            for (const auto &cellsLists2 : neighborCellsLists) {
              if (cellsLists2.size() > i) {
                listSize += cellsLists2[i].second.size();
              }
            }
            return listSize;
          }
        }
      }
    }
    return 0lu;
  }

  /**
   * Returns the neighbor list in AoS layout.
   * @return Neighbor list in AoS layout.
   */
  typename VerletListsCellsHelpers::PairwiseNeighborListsType<Particle> &getAoSNeighborList() {
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

  /**
   * @copydoc VLCNeighborListInterface::buildAoSNeighborList()
   */
  void buildAoSNeighborList(TraversalOption vlcTraversalOpt, LinkedCells<Particle> &linkedCells, bool useNewton3) {
    using namespace utils::ArrayMath::literals;
    // Sanity check.
    if (linkedCells.getCellBlock().getCellsPerInteractionLength() > 1) {
      utils::ExceptionHandler::exception(
          "VLCCellPairNeighborList::rebuildNeighborListsVLC() was called with a CSF < 1 but it only supports CSF>=1.");
    }
    // Define some aliases
    auto &neighborLists = getAoSNeighborList();
    auto &globalToLocalIndex = _globalToLocalIndex;
    auto &cells = linkedCells.getCells();
    const auto interactionLength = linkedCells.getInteractionLength();
    const auto interactionLengthSquared = interactionLength * interactionLength;
    const auto boxSizeWithHalo = linkedCells.getBoxMax() - linkedCells.getBoxMin() +
                                 std::array<double, 3>{interactionLength, interactionLength, interactionLength} * 2.;

    // Helper lambda to compute the relative index from two cells
    auto relIdx = [&](auto cellIndex1, auto cellIndex2) {
      const auto cellsPerDimensionWithHalo = linkedCells.getCellBlock().getCellsPerDimensionWithHalo();
      const auto threeDPosCell1 = utils::ThreeDimensionalMapping::oneToThreeD(static_cast<long unsigned>(cellIndex1),
                                                                              cellsPerDimensionWithHalo);
      const auto threeDPosCell2 = utils::ThreeDimensionalMapping::oneToThreeD(static_cast<long unsigned>(cellIndex2),
                                                                              cellsPerDimensionWithHalo);
      const auto offset = threeDPosCell2 - threeDPosCell1;
      return (offset[0] + 1) * 9 + (offset[1] + 1) * 3 + (offset[2] + 1);
    };

    // Todo reuse old neighbour lists
    neighborLists.clear();
    neighborLists.resize(cells.size());

    globalToLocalIndex.clear();
    globalToLocalIndex.resize(cells.size());

    // Todo: Create an estimate for the average length of a neighbor list.
    // This assumes homogeneous distribution and some overestimation.
    const auto listLengthEstimate = VerletListsCellsHelpers::estimateListLength(
        linkedCells.getNumberOfParticles(IteratorBehavior::ownedOrHalo), boxSizeWithHalo, interactionLength, 1.3);

    // Todo reuse old neighbour lists
    // Reset lists. Don't free any memory, only mark as unused.
    this->setLinkedCellsPointer(&linkedCells);
    for (auto &neighborCellLists : neighborLists) {
      for (auto &cellLists : neighborCellLists) {
        for (auto &[particlePtr, neighbors] : cellLists) {
          particlePtr = nullptr;
          neighbors.clear();
        }
      }
    }
    neighborLists.resize(cells.size());
    for (auto &neighborList : neighborLists) {
      neighborList.resize(27);
    }

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
        if (neighborList.size() <= p1Index) {
          neighborList.resize(p1Index + 1);
        }
        neighborList[p1Index].first = &p1;
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
    const auto offsets = VerletListsCellsHelpers::buildBaseStep(cellsPerDim, vlcTraversalOpt);

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
          auto threadNum = autopas_get_thread_num();

          // Todo: Use esitmate from above for this
          // Allocate memory for ptr-list pairs for this cell.
          // baseCellsLists.resize(27);

          // baseCellsLists.resize(VerletListsCellsHelpers::estimateNumLists(
          //     cellIndexBase, useNewton3, cells, offsets,
          //     utils::ArrayUtils::static_cast_copy_array<size_t>(cellsPerDim)));

          // Re-initialize a neighbor list for all particles in the cell but at most for as many as there are lists
          // const size_t minCellSizeVsNumLists = std::min(baseCell.size(), baseCellsLists.size());
          // for (size_t i = 0; i < minCellSizeVsNumLists; ++i) {
          //   auto &[particlePtr, neighbors] = baseCellsLists[i];
          //   particlePtr = &baseCell[i];
          //   neighbors.reserve(listLengthEstimate);
          // }

          // For any remaining particles create a new list.
          // This case can only happen if estimateNumLists can return values smaller than baseCell.size()
          // for (size_t i = minCellSizeVsNumLists; i < baseCell.size(); ++i) {
          //   baseCellsLists.emplace_back(&baseCell[i], std::vector<Particle *>{});
          //   baseCellsLists.back().second.reserve(listLengthEstimate);
          // }

          // Build c08 lists for this base step according to predefined cell pairs
          for (const auto &[offset1, offset2, _] : offsets) {
            const auto cellIndex1 = cellIndexBase + offset1;
            const auto cellIndex2 = cellIndexBase + offset2;
            auto &cell1List = neighborLists[cellIndex1];
            auto &cell2List = neighborLists[cellIndex2];

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
              if (cellIndex1 == cellIndex2 and vlcTraversalOpt != TraversalOption::vlp_c01) {
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
                  {
                    size_t secondCellIndexInFirst;
                    insert(p1, particleIndexCell1, p2, cellIndex1, cellIndex1,
                           cell1List[relIdx(cellIndex1, cellIndex2)]);
                  }
                  // If the traversal does not use Newton3 the inverse interaction also needs to be stored in p2's list
                  if (not useNewton3 and not(vlcTraversalOpt == TraversalOption::vlp_c01)) {
                    {
                      size_t secondCellIndexInFirst;
                      insert(p2, particleIndexCell2, p1, cellIndex2, cellIndex2,
                             cell2List[relIdx(cellIndex2, cellIndex1)]);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    // Cleanup: Remove any unused ptr-list pairs to avoid accessing nullptr
    for (auto &neighborCellLists : neighborLists) {
      for (auto &cellLists : neighborCellLists) {
        cellLists.erase(std::remove_if(cellLists.begin(), cellLists.end(),
                                       [](const auto &pair) {
                                         const auto &[particlePtr, neighbors] = pair;
                                         return particlePtr == nullptr;
                                       }),
                        cellLists.end());
      }
    }
  }

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
      soaLists.resize(aosLists.size());

      // iterate over each cell's neighboring cells
      for (size_t secondCellIndex = 0; secondCellIndex < aosLists.size(); ++secondCellIndex) {
        const auto &aosCellPairLists = aosLists[secondCellIndex];
        auto &soaCellPairLists = soaLists[secondCellIndex];
        soaCellPairLists.reserve(aosCellPairLists.capacity());

        // iterate over pairs of particle and neighbor list
        for (const auto &[particlePtr, neighbors] : aosCellPairLists) {
          // global index of current particle
          size_t currentParticleGlobalIndex = particlePtrToIndex.at(particlePtr);

          // create SoA neighbor list for current particle
          std::vector<size_t, autopas::AlignedAllocator<size_t>> currentSoANeighborList{};
          currentSoANeighborList.reserve(neighbors.size());

          // fill the SoA neighbor list with the indices of the particles from the corresponding AoS neighbor list
          for (const auto &neighborOfCurrentParticle : neighbors) {
            currentSoANeighborList.emplace_back(particlePtrToIndex.at(neighborOfCurrentParticle));
          }

          // add the newly constructed pair of particle index and SoA particle neighbor list to cell pair
          soaCellPairLists.emplace_back(currentParticleGlobalIndex, currentSoANeighborList);
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
  /**
   * Internal neighbor list structure in AoS format - Verlet lists for each particle for each cell pair.
   */
  typename VerletListsCellsHelpers::PairwiseNeighborListsType<Particle> _aosNeighborList =
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
