/**
 * @file VerletListsCells.h
 * @author nguyen
 * @date 30.08.18
 */

#pragma once

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/CellBasedParticleContainer.h"
#include "autopas/containers/LoadEstimators.h"
#include "autopas/containers/cellPairTraversals/BalancedTraversal.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/verletListsCellBased/VerletListsLinkedBase.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsHelpers.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/tuning/selectors/TraversalSelector.h"

namespace autopas {

/**
 * Linked Cells with Verlet Lists container.
 * The VerletListsCells class uses neighborhood lists for each cell
 * to calculate pairwise interactions of particles.
 * It is optimized for a constant, i.e. particle independent, cutoff radius of
 * the interaction.
 * Cells are created using a cell size of at least cutoff + skinPerTimestep*rebuildFrequency.
 * @tparam Particle
 * @tparam NeighborList The neighbor list used by this container.
 */

template <class Particle, class NeighborList>
class VerletListsCells : public VerletListsLinkedBase<Particle> {
  using ParticleCell = FullParticleCell<Particle>;

 public:
  /**
   * Constructor of the VerletListsCells class.
   * The neighbor lists are build using a search radius of cutoff + skin*rebuildfrequency.
   * @param boxMin The lower corner of the domain.
   * @param boxMax The upper corner of the domain.
   * @param cutoff The cutoff radius of the interaction.
   * @param rebuildFrequency The rebuild Frequency.
   * @param skinPerTimestep The skin radius per Timestep.
   * @param cellSizeFactor Cell size factor relative to cutoff.
   * @param loadEstimator Load estimation algorithm for balanced traversals.
   * @param dataLayoutDuringListRebuild Data layout during the list generation. Has no influence on list layout.
   */
  VerletListsCells(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, const double cutoff,
                   const double skinPerTimestep = 0, const unsigned int rebuildFrequency = 2,
                   const double cellSizeFactor = 1.0,
                   const LoadEstimatorOption loadEstimator = LoadEstimatorOption::squaredParticlesPerCell,
                   typename VerletListsCellsHelpers::VLCBuildType dataLayoutDuringListRebuild =
                       VerletListsCellsHelpers::VLCBuildType::soaBuild)
      : VerletListsLinkedBase<Particle>(boxMin, boxMax, cutoff, skinPerTimestep, rebuildFrequency,
                                        compatibleTraversals::allVLCCompatibleTraversals(), cellSizeFactor),
        _loadEstimator(loadEstimator),
        _dataLayoutDuringListRebuild(dataLayoutDuringListRebuild) {}

  /**
   * @copydoc ParticleContainerInterface::getContainerType()
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return _neighborList.getContainerType(); }

  /**
   * Generates the load estimation function depending on _loadEstimator.
   * @return load estimator function object.
   */
  BalancedTraversal::EstimatorFunction getLoadEstimatorFunction() {
    switch (this->_loadEstimator) {
      case LoadEstimatorOption::squaredParticlesPerCell: {
        return [&](const std::array<unsigned long, 3> &cellsPerDimension,
                   const std::array<unsigned long, 3> &lowerCorner, const std::array<unsigned long, 3> &upperCorner) {
          return loadEstimators::squaredParticlesPerCell((this->_linkedCells).getCells(), cellsPerDimension,
                                                         lowerCorner, upperCorner);
        };
      }
      case LoadEstimatorOption::neighborListLength: {
        return [&](const std::array<unsigned long, 3> &cellsPerDimension,
                   const std::array<unsigned long, 3> &lowerCorner, const std::array<unsigned long, 3> &upperCorner) {
          return loadEstimators::neighborListLength<Particle, NeighborList>(_neighborList, cellsPerDimension,
                                                                            lowerCorner, upperCorner);
        };
      }

      case LoadEstimatorOption::none:
        [[fallthrough]];
      default: {
        return
            [&](const std::array<unsigned long, 3> &cellsPerDimension, const std::array<unsigned long, 3> &lowerCorner,
                const std::array<unsigned long, 3> &upperCorner) { return 1; };
      }
    }
  }

  void iteratePairwise(TraversalInterface *traversal) override {
    // Check if traversal is allowed for this container and give it the data it needs.
    _neighborList.setUpTraversal(traversal);
    if (auto *balancedTraversal = dynamic_cast<BalancedTraversal *>(traversal)) {
      balancedTraversal->setLoadEstimator(getLoadEstimatorFunction());
    }

    traversal->initTraversal();
    traversal->traverseParticlePairs();
    traversal->endTraversal();
  }

  /**
   * Gets the number of neighbors over all neighbor lists that belong to this particle.
   * @param particle
   * @return the size of the neighbor list(s) of this particle
   */
  size_t getNumberOfPartners(const Particle *particle) const { return _neighborList.getNumberOfPartners(particle); }

  /**
   * Special case of building the neighbor lists in c08 and c18 style where all lists that belong to one base step are
   * stored together.
   * @param traversal The TraversalOption with which the function is called.
   */
  void rebuildNeighborListsC08C18(TraversalOption traversal) {
    using namespace utils::ArrayMath::literals;
    // So far this neighborlist rebuilding only supports c01, c08 c18, sliced traversals.
    if (traversal != TraversalOption::vlc_c08 and traversal != TraversalOption::vlc_c18 and
        traversal != TraversalOption::vlc_c01 and traversal != TraversalOption::vlc_sliced and
        traversal != TraversalOption::vlc_sliced_balanced and traversal != TraversalOption::vlc_sliced_c02) {
      utils::ExceptionHandler::exception(
          "VerletListsCells::rebuildNeighborListsC08C18() was called with an unsupported traversal.");
    }
    // Sanity check.
    if (this->_linkedCells.getCellBlock().getCellsPerInteractionLength() > 1) {
      utils::ExceptionHandler::exception(
          "VerletListsCells::rebuildNeighborListsC08C18() was called with a CSF < 1 but it only supports CSF>=1.");
    }
    // Define some aliases
    auto &neighborLists = _neighborList.getAoSNeighborList();
    auto &particleToCellMap = _neighborList.getParticleToCellMap();
    auto &cells = this->_linkedCells.getCells();
    const auto interactionLength = this->getInteractionLength();
    const auto interactionLengthSquared = interactionLength * interactionLength;
    const auto boxSizeWithHalo = this->_linkedCells.getBoxMax() - this->_linkedCells.getBoxMin() +
                                 std::array<double, 3>{interactionLength, interactionLength, interactionLength} * 2.;
    // Create an estimate for the average length of a neighbor list.
    // This assumes homogeneous distribution and some overestimation.
    const auto listLengthEstimate = VerletListsCellsHelpers::estimateListLength(
        this->_linkedCells.getNumberOfParticles(IteratorBehavior::ownedOrHalo), boxSizeWithHalo, interactionLength,
        1.3);

    // Reset lists. Don't free any memory, only mark as unused.
    _neighborList.setLinkedCellsPointer(&this->_linkedCells);
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

    const auto &cellsPerDim = utils::ArrayUtils::static_cast_copy_array<int>(
        this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo());
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
              cellIndexBase, this->_verletBuiltNewton3, cells, offsets,
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

            if (traversal == TraversalOption::vlc_c18 or traversal == TraversalOption::vlc_sliced or
                traversal == TraversalOption::vlc_sliced_balanced or traversal == TraversalOption::vlc_sliced_c02) {
              const auto cell2Coords = utils::ThreeDimensionalMapping::oneToThreeD(cellIndex2, cellsPerDim);
              if (cell2Coords[0] >= cellsPerDim[0] or cell2Coords[1] >= cellsPerDim[1] or
                  cell2Coords[2] >= cellsPerDim[2]) {
                continue;
              }
            }

            if (traversal == TraversalOption::vlc_c01) {
              const auto cell2Coords = utils::ThreeDimensionalMapping::oneToThreeD(cellIndex2, cellsPerDim);
              if (cell2Coords[0] >= cellsPerDim[0] or cell2Coords[0] < 0 or cell2Coords[1] >= cellsPerDim[1] or
                  cell2Coords[1] < 0 or cell2Coords[2] >= cellsPerDim[2] or cell2Coords[2] < 0) {
                continue;
              }
            }

            // Skip if both cells only contain halos
            if (not(cells[cellIndex1].getPossibleParticleOwnerships() == OwnershipState::owned) and
                not(cells[cellIndex2].getPossibleParticleOwnerships() == OwnershipState::owned)) {
              continue;
            }

            // Go over all particle pairs in the two cells and insert close pairs into their respective lists
            for (size_t particleIndexCell1 = 0; particleIndexCell1 < cells[cellIndex1].size(); ++particleIndexCell1) {
              auto &p1 = cells[cellIndex1][particleIndexCell1];
              for (size_t particleIndexCell2 = (cellIndex1 == cellIndex2 and traversal != TraversalOption::vlc_c01)
                                                   ? particleIndexCell1 + 1
                                                   : 0;
                   particleIndexCell2 < cells[cellIndex2].size(); ++particleIndexCell2) {
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
                  particleToCellMap[&p1] = std::make_pair(cellIndex1, particleIndexCell1);
                  particleToCellMap[&p2] = std::make_pair(cellIndex2, particleIndexCell2);
                  // If the traversal does not use Newton3 the inverse interaction also needs to be stored in p2's list
                  if (not this->_verletBuiltNewton3 and not(traversal == TraversalOption::vlc_c01)) {
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

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    this->_verletBuiltNewton3 = traversal->getUseNewton3();

    // VLP needs constexpr special case because the types and thus interfaces are slightly different
    if constexpr (std::is_same_v<NeighborList, VLCCellPairNeighborList<Particle>>) {
      _neighborList.buildAoSNeighborList(this->_linkedCells, this->_verletBuiltNewton3, this->getCutoff(),
                                         this->getVerletSkin(), this->getInteractionLength(),
                                         traversal->getTraversalType(), _dataLayoutDuringListRebuild);

    } else {
      // Switch to test different implementations for vlc_c08 list generation
      const auto traversalOption = traversal->getTraversalType();
      if (traversalOption == TraversalOption::vlc_c08 or traversalOption == TraversalOption::vlc_c18 or
          traversalOption == TraversalOption::vlc_c01 or traversalOption == TraversalOption::vlc_sliced or
          traversalOption == TraversalOption::vlc_sliced_balanced or
          traversalOption == TraversalOption::vlc_sliced_c02) {
        rebuildNeighborListsC08C18(traversalOption);
      } else {
        _neighborList.buildAoSNeighborList(this->_linkedCells, this->_verletBuiltNewton3, this->getCutoff(),
                                           this->getVerletSkin(), this->getInteractionLength(),
                                           traversal->getTraversalType(), _dataLayoutDuringListRebuild);
      }
    }

    if (traversal->getDataLayout() == DataLayoutOption::soa) {
      _neighborList.generateSoAFromAoS(this->_linkedCells);
    }

    // the neighbor list is now valid
    this->_neighborListIsValid.store(true, std::memory_order_relaxed);
  }

  /**
   * Return the cell length of the underlying linked cells structure, normally needed only for unit tests.
   * @return
   */
  [[nodiscard]] const std::array<double, 3> &getCellLength() const {
    return this->_linkedCells.getCellBlock().getCellLength();
  }

 private:
  /**
   * Neighbor list abstraction for neighbor list used in the container.
   */
  NeighborList _neighborList;

  /**
   * Load estimation algorithm for balanced traversals.
   */
  autopas::LoadEstimatorOption _loadEstimator;

  /**
   * Data layout during the list generation. Has no influence on list layout.
   */
  typename VerletListsCellsHelpers::VLCBuildType _dataLayoutDuringListRebuild;
};
}  // namespace autopas
