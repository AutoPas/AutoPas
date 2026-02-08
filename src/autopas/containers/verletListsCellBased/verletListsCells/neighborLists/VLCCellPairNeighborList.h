/**
 * @file VLCCellPairNeighborList.h
 * @author tirgendetwas
 * @date 07.11.20
 */

#pragma once

#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsHelpers.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/VLCNeighborListInterface.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCCellPairTraversalInterface.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCTraversalInterface.h"

namespace autopas {

/**
 * Neighbor list to be used with VerletListsCells container.
 * Pairwise verlet lists iterates through each pair of neighboring cells
 * and generates a neighbor list for each particle from cell1, which consists of its (potential) partners from cell2.
 * @tparam Particle_T Type of particle to be used for this neighbor list.
 */
template <class Particle_T>
class VLCCellPairNeighborList : public VLCNeighborListInterface<Particle_T> {
 public:
  /**
   * Type of the data structure used to save the neighbor lists.
   */
  using ListType = typename VerletListsCellsHelpers::PairwiseNeighborListsType<Particle_T>;

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

  size_t getNumberOfPartners(const Particle_T *particle) const override {
    size_t listSize{0};
    bool particleFound{false};
    for (auto &neighborListsForOneCell : _aosNeighborList) {
      for (const auto &neighborListsForCellCellPair : neighborListsForOneCell) {
        // Check for the desired particle in the first cell. We check all "partner" cells because the particle might not
        // have any particles from one "partner" cell in its list.
        for (size_t i{0}; i < neighborListsForCellCellPair.size(); ++i) {
          if (neighborListsForCellCellPair[i].first == particle) {
            particleFound = true;
            // We accumulate the number of partners. Partners can be in all lists of neighboring cells (middle for loop)
            listSize += neighborListsForCellCellPair[i].second.size();
            // Since we found the particle, we can skip all the other particles in the current list
            break;
          }
        }
      }
      if (particleFound) {
        // We've found the particle in the cell with neighbor lists neighborListsForOneCell. So we can skip all the
        // other cells (outer for loop)
        return listSize;
      }
    }
    return 0lu;
  }

  /**
   * Returns the neighbor list in AoS layout.
   * @return Neighbor list in AoS layout.
   */
  typename VerletListsCellsHelpers::PairwiseNeighborListsType<Particle_T> &getAoSNeighborList() {
    return _aosNeighborList;
  }

  /**
   * Returns the neighbor list in SoA layout.
   * @return Neighbor list in SoA layout.
   */
  auto &getSoANeighborList() { return _soaNeighborList; }

  /**
   * @copydoc VLCNeighborListInterface::buildAoSNeighborList()
   */
  void buildAoSNeighborList(TraversalOption vlcTraversalOpt, LinkedCells<Particle_T> &linkedCells,
                            bool useNewton3) override {
    using namespace utils::ArrayMath::literals;
    // Sanity check.
    if (linkedCells.getCellBlock().getCellsPerInteractionLength() > 1) {
      utils::ExceptionHandler::exception(
          "VLCCellPairNeighborList::buildAoSNeighborList() was called with a CSF < 1 but it only supports CSF>=1.");
    }
    // Define some aliases
    auto &neighborLists = getAoSNeighborList();
    auto &cells = linkedCells.getCells();
    const auto interactionLength = linkedCells.getInteractionLength();
    const auto interactionLengthSquared = interactionLength * interactionLength;
    const auto boxSizeWithHalo = linkedCells.getBoxMax() - linkedCells.getBoxMin() +
                                 std::array<double, 3>{interactionLength, interactionLength, interactionLength} * 2.;

    // Helper lambda to compute the relative index from two cells within in a 3x3x3 cell-cube
    auto relativeNeighborhoodIndex = [&](auto cellIndex1, auto cellIndex2) {
      const auto cellsPerDimensionWithHalo = linkedCells.getCellBlock().getCellsPerDimensionWithHalo();
      const auto threeDPosCell1 = utils::ThreeDimensionalMapping::oneToThreeD(static_cast<long unsigned>(cellIndex1),
                                                                              cellsPerDimensionWithHalo);
      const auto threeDPosCell2 = utils::ThreeDimensionalMapping::oneToThreeD(static_cast<long unsigned>(cellIndex2),
                                                                              cellsPerDimensionWithHalo);
      const auto offset = threeDPosCell2 - threeDPosCell1;
      return (offset[0] + 1) * 9 + (offset[1] + 1) * 3 + (offset[2] + 1);
    };

    // This assumes homogeneous distribution and some overestimation.
    const auto listLengthEstimate = VerletListsCellsHelpers::estimateListLength(
        linkedCells.getNumberOfParticles(IteratorBehavior::ownedOrHalo), boxSizeWithHalo, interactionLength, 1.3);

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
    // each cell hast 27 neighbors including itself
    for (auto &neighborList : neighborLists) {
      neighborList.resize(27);
    }

    /* This must not be a doc comment (with two **) to not confuse doxygen.
     * Helper function to insert a pointer into a list of the base cell.
     * It considers the cases that neither particle is in the base cell
     * and in that case finds or creates the appropriate list.
     *
     * @param p1 Reference to source particle.
     * @param p2 Reference to target particle.
     * @param neighborList Reference to the list where the particle pair should be stored.
     */
    auto insert = [&](auto &p1, auto &p2, auto &neighborList) {
      // Check if the base cell already has a list for p1
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
          neighborList.emplace_back(&p1, std::vector<Particle_T *>{});
          neighborList.back().second.reserve(listLengthEstimate);
          neighborList.back().second.push_back(&p2);
        }
      }
    };

    const auto &cellsPerDim =
        utils::ArrayUtils::static_cast_copy_array<int>(linkedCells.getCellBlock().getCellsPerDimensionWithHalo());
    // Vector of offsets from the base cell for the given base step
    // and respective factors for the fraction of particles per cell that need neighbor lists in the base cell.
    const auto offsets = VerletListsCellsHelpers::buildBaseStep(cellsPerDim, vlcTraversalOpt);

    int xEnd{};
    int yEnd{};
    int zEnd{};

    switch (vlcTraversalOpt) {
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
    const auto numThreads = autopas_get_preferred_num_threads();
    AUTOPAS_OPENMP(parallel for collapse(3) schedule(dynamic, std::max(cells.size() / (numThreads * 10), 1ul)) \
      num_threads(numThreads) \
    )
    for (int z = 0; z < zEnd; ++z) {
      for (int y = 0; y < yEnd; ++y) {
        for (int x = 0; x < xEnd; ++x) {
          // aliases
          const auto cellIndexBase = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, cellsPerDim);
          auto &baseCell = cells[cellIndexBase];
          auto &baseCellsLists = neighborLists[cellIndexBase];
          auto threadNum = autopas_get_thread_num();

          // Build lists for this base step according to predefined cell pairs
          for (const auto &[offset1, offset2, _] : offsets) {
            const auto cell1Index1D = cellIndexBase + offset1;
            const auto cell2Index1D = cellIndexBase + offset2;
            auto &cell1List = neighborLists[cell1Index1D];
            auto &cell2List = neighborLists[cell2Index1D];

            // For all traversals ensures the partner cell is not outside the boundary
            const auto cell2Index3D = utils::ThreeDimensionalMapping::oneToThreeD(cell2Index1D, cellsPerDim);
            if (cell2Index3D[0] >= cellsPerDim[0] or cell2Index3D[0] < 0 or cell2Index3D[1] >= cellsPerDim[1] or
                cell2Index3D[1] < 0 or cell2Index3D[2] >= cellsPerDim[2] or cell2Index3D[2] < 0) {
              continue;
            }

            // Skip if both cells only contain halos or dummys
            if (not(cells[cell1Index1D].getPossibleParticleOwnerships() == OwnershipState::owned) and
                not(cells[cell2Index1D].getPossibleParticleOwnerships() == OwnershipState::owned)) {
              continue;
            }

            // Go over all particle pairs in the two cells and insert close pairs into their respective lists
            for (size_t particleIndexCell1 = 0; particleIndexCell1 < cells[cell1Index1D].size(); ++particleIndexCell1) {
              auto &p1 = cells[cell1Index1D][particleIndexCell1];

              // Determine the starting index for the second particle in the pair.
              // If both cells are the same and the traversal is newton3 compatible, this is the next particle in the
              // cell. If the cells are different or the traversal is not newton3 compatible, this is index 0. For
              // non-newton 3 compatible traversals (vlp_c01) we have to consider the interaction in both directions, so
              // we always start from index 0.
              size_t startIndexCell2 = 0;
              if (cell1Index1D == cell2Index1D and vlcTraversalOpt != TraversalOption::vlp_c01) {
                startIndexCell2 = particleIndexCell1 + 1;
              }

              for (size_t particleIndexCell2 = startIndexCell2; particleIndexCell2 < cells[cell2Index1D].size();
                   ++particleIndexCell2) {
                auto &p2 = cells[cell2Index1D][particleIndexCell2];
                // Ignore dummies and self interaction
                if (&p1 == &p2 or p1.isDummy() or p2.isDummy()) {
                  continue;
                }

                // If the distance is less than interaction length add the pair to the list
                const auto distVec = p2.getR() - p1.getR();
                const auto distSquared = utils::ArrayMath::dot(distVec, distVec);
                if (distSquared < interactionLengthSquared) {
                  {
                    const size_t secondCellIndexInFirst = relativeNeighborhoodIndex(cell1Index1D, cell2Index1D);
                    insert(p1, p2, cell1List[secondCellIndexInFirst]);
                  }
                  // If the traversal is Newton3 compatible, Newton3 is used for building regardless of if the actual
                  // interactions will use Newton3. In the case that Newton3 will not be used, the inverse interaction
                  // also needs to be stored in p2's list. If the traversal is Newton3 incompatible, the insertion into
                  // p2's list will occur when the p2 particle is p1. This is ensured above by the startIndexCell2.
                  if (not useNewton3 and not(vlcTraversalOpt == TraversalOption::vlp_c01)) {
                    {
                      const size_t secondCellIndexInFirst = relativeNeighborhoodIndex(cell2Index1D, cell1Index1D);
                      insert(p2, p1, cell2List[secondCellIndexInFirst]);
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

  void generateSoAFromAoS(LinkedCells<Particle_T> &linkedCells) override {
    _soaNeighborList.clear();

    // particle pointer to global index of particle
    std::unordered_map<Particle_T *, size_t> particlePtrToIndex;
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
    auto vTraversal = dynamic_cast<VLCCellPairTraversalInterface<Particle_T> *>(traversal);

    if (vTraversal) {
      vTraversal->setVerletList(*this);
    } else {
      auto traversal2 =
          dynamic_cast<VLCTraversalInterface<Particle_T, VLCCellPairNeighborList<Particle_T>> *>(traversal);
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
  typename VerletListsCellsHelpers::PairwiseNeighborListsType<Particle_T> _aosNeighborList =
      std::vector<std::vector<std::vector<std::pair<Particle_T *, std::vector<Particle_T *>>>>>();

  /**
   * Internal neighbor list structure in SoA format - Verlet lists for each particle for each cell pair.
   * Contrary to aosNeighborList it saves global particle indices instead of particle pointers.
   */
  SoAListType _soaNeighborList = std::vector<
      std::vector<std::vector<std::pair<size_t, std::vector<size_t, autopas::AlignedAllocator<size_t>>>>>>();
};
}  // namespace autopas
