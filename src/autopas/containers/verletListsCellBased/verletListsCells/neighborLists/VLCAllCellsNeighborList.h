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

#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsHelpers.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/VLCNeighborListInterface.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCTraversalInterface.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

/**
 * Neighbor list to be used with VerletListsCells container. Classic implementation of verlet lists based on linked
 * cells.
 * @tparam Particle_T Type of particle to be used for this neighbor list.
 */
template <class Particle_T>
class VLCAllCellsNeighborList : public VLCNeighborListInterface<Particle_T> {
 public:
  /**
   * Type of the data structure used to save the neighbor lists.
   */
  using ListType = typename VerletListsCellsHelpers::AllCellsNeighborListsType<Particle_T>;

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
  void buildAoSNeighborList(TraversalOption vlcTraversalOpt, LinkedCells<Particle_T> &linkedCells,
                            bool useNewton3) override {
    using namespace utils::ArrayMath::literals;
    // Sanity check.
    if (linkedCells.getCellBlock().getCellsPerInteractionLength() > 1) {
      utils::ExceptionHandler::exception(
          "VLCAllCellsNeighborLists::buildAoSNeighborList() was called with a CSF < 1 but it only supports CSF>=1.");
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
            neighborList.emplace_back(&p1, std::vector<Particle_T *>{});
            neighborList.back().second.reserve(listLengthEstimate);
            neighborList.back().second.push_back(&p2);
          }
        }
      }
    };

    const auto &cellsPerDim =
        utils::ArrayUtils::static_cast_copy_array<int>(linkedCells.getCellBlock().getCellsPerDimensionWithHalo());
    // Vector of offsets from the base cell for the base step corresponding to the traversal
    // and respective factors for the fraction of particles per cell that need neighbor lists in the base cell.
    const auto offsets = VerletListsCellsHelpers::buildBaseStep(cellsPerDim, vlcTraversalOpt);

    int xEnd{};
    int yEnd{};
    int zEnd{};

    switch (vlcTraversalOpt) {
      case TraversalOption::vlc_c08:
        // cellsPerDim includes halo cells, so outermost cells are halo.
        // Only owned-halo interactions are needed, not halo-halo.
        // c08 is forward-looking: halo base cells only interact with other halo cells -> can skip.
        // Other traversals like C01 and C18 are backward-looking: halo base cells may interact with owned cells -> must
        // keep.
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

          // Allocate memory for ptr-list pairs for this cell.
          baseCellsLists.resize(VerletListsCellsHelpers::estimateNumLists(
              cellIndexBase, useNewton3, cells, offsets,
              utils::ArrayUtils::static_cast_copy_array<size_t>(cellsPerDim)));
          // Re-initialize neighbor lists for all particles in the cell, or for as many lists as currently exist if this
          // is smaller. By "re-initialize", we mean reserving the estimated list length.
          const size_t minNumParticlesVsNumLists = std::min(baseCell.size(), baseCellsLists.size());
          for (size_t i = 0; i < minNumParticlesVsNumLists; ++i) {
            auto &[particlePtr, neighbors] = baseCellsLists[i];
            particlePtr = &baseCell[i];
            neighbors.reserve(listLengthEstimate);
          }
          // For any remaining particles without (size non-zero) neighbor lists, reserve the estimated list length.
          for (size_t i = minNumParticlesVsNumLists; i < baseCell.size(); ++i) {
            baseCellsLists.emplace_back(&baseCell[i], std::vector<Particle_T *>{});
            baseCellsLists.back().second.reserve(listLengthEstimate);
          }

          // Build lists for this base step according to predefined cell pairs
          for (const auto &[offset1, offset2, _] : offsets) {
            const auto cell1Index1D = cellIndexBase + offset1;
            const auto cell2Index1D = cellIndexBase + offset2;

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
              // non-newton 3 compatible traversals (vlc_c01) we have to consider the interaction in both directions, so
              // we always start from index 0.
              size_t startIndexCell2 = 0;
              if (cell1Index1D == cell2Index1D and vlcTraversalOpt != TraversalOption::vlc_c01) {
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
                  insert(p1, particleIndexCell1, p2, cell1Index1D, cellIndexBase, baseCellsLists);
                  // If the traversal is Newton3 compatible, Newton3 is used for building regardless of if the actual
                  // interactions will use Newton3. In the case that Newton3 will not be used, the inverse interaction
                  // also needs to be stored in p2's list. If the traversal is Newton3 incompatible, the insertion into
                  // p2's list will occur when the p2 particle is p1. This is ensured above by the startIndexCell2.
                  if (not useNewton3 and not(vlcTraversalOpt == TraversalOption::vlc_c01)) {
                    insert(p2, particleIndexCell2, p1, cell2Index1D, cellIndexBase, baseCellsLists);
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
   * @copydoc VLCNeighborListInterface::getNumberOfPartners()
   */
  size_t getNumberOfPartners(const Particle_T *particle) const override {
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
  typename VerletListsCellsHelpers::AllCellsNeighborListsType<Particle_T> &getAoSNeighborList() {
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
    auto vTraversal = dynamic_cast<VLCTraversalInterface<Particle_T, VLCAllCellsNeighborList<Particle_T>> *>(traversal);

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
   * Internal neighbor list structure in AoS format - Verlet lists for each particle for each cell.
   */
  typename VerletListsCellsHelpers::AllCellsNeighborListsType<Particle_T> _aosNeighborList{};

  /**
   * Internal neighbor list structure in SoA format - Verlet lists for each particle for each cell.
   * Contrary to aosNeighborList it saves global particle indices instead of particle pointers.
   */
  SoAListType _soaNeighborList{};
};
}  // namespace autopas
