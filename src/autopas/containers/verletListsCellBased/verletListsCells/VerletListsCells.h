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

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    this->_verletBuiltNewton3 = traversal->getUseNewton3();

    // VLP needs constexpr special case because the types and thus interfaces are slightly different
    if constexpr (std::is_same_v<NeighborList, VLCCellPairNeighborList<Particle>>) {
      _neighborList.buildAoSNeighborList(this->_linkedCells, this->_verletBuiltNewton3, this->getCutoff(),
                                         this->getVerletSkin(), this->getInteractionLength(),
                                         traversal->getTraversalType(), _dataLayoutDuringListRebuild);

    } else {
// Switch to test two implementations for vlc_c08 list generation
#define BUILD_WITHOUT_FUNCTOR
#ifdef BUILD_WITHOUT_FUNCTOR
      if (traversal->getTraversalType() == TraversalOption::vlc_c08) {
        const auto &cellsPerDim = utils::ArrayUtils::static_cast_copy_array<int>(
            this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo());

        // list of cell interaction pairs relative to the base cell
        const auto offsets = [&]() {  // NOLINT(*-function-cognitive-complexity) TODO fix this
          // cellOffset 1, cellOffset2, list estimation factor
          std::vector<std::tuple<int, int, double>> offsets;
          offsets.reserve(14);
          // This is currently guaranteed by the constructor
          const int interactionCellsPerDim = 1;
          constexpr std::array<double, 4> estimatorFactors{1., 1., 1. / 4. * M_PI, 1. / 6. * M_PI};
//          constexpr std::array<double, 4> estimatorFactors{1., 1. * 0.8, 1. / 4. * M_PI * 0.7, 1. / 6. * M_PI * 0.5};
          for (int z = -interactionCellsPerDim; z <= interactionCellsPerDim; ++z) {
            for (int y = -interactionCellsPerDim; y <= interactionCellsPerDim; ++y) {
              for (int x = -interactionCellsPerDim; x <= interactionCellsPerDim; ++x) {
                int baseCell = 0;
                int partnerCell = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, cellsPerDim);
                if (x < 0) {
                  baseCell -= x;
                  partnerCell -= x;
                }
                if (y < 0) {
                  baseCell -= y * cellsPerDim[0];
                  partnerCell -= y * cellsPerDim[0];
                }
                if (z < 0) {
                  baseCell -= z * cellsPerDim[0] * cellsPerDim[1];
                  partnerCell -= z * cellsPerDim[0] * cellsPerDim[1];
                }
                // Count number of non-aligned dimensions
                const auto factor = estimatorFactors[std::abs(x) + std::abs(y) + std::abs(z)];
                size_t smallerIndex, biggerIndex;
                std::tie(smallerIndex, biggerIndex) = std::minmax(baseCell, partnerCell);
                if (std::find_if(offsets.begin(), offsets.end(), [&](const auto &tuple) {
                      return std::get<0>(tuple) == smallerIndex and std::get<1>(tuple) == biggerIndex;
                    }) == offsets.end()) {
                  offsets.emplace_back(smallerIndex, biggerIndex, factor);
                }
              }
            }
          }
          // sort offsets to group processing of the same cells (-> better cache re-usage)
          std::sort(offsets.begin(), offsets.end(), [](const auto &pair1, const auto &pair2) {
            const auto &[a1, b1, f1] = pair1;
            const auto &[a2, b2, f2] = pair2;

            if (a1 == a2) {
              return b1 < b2;
            } else {
              return a1 < a2;
            }
          });
          return offsets;
        }();

        using namespace utils::ArrayMath::literals;
        // Reset lists and define some aliases
        auto &neighborLists = _neighborList.getAoSNeighborList();
        neighborLists.clear();
        _neighborList.setLinkedCells(&this->_linkedCells);
        auto &cells = this->_linkedCells.getCells();
        const auto interactionLength = this->getInteractionLength();
        const auto interactionLengthSquared = interactionLength * interactionLength;
        const auto boxSizeWithHalo =
            this->_linkedCells.getBoxMax() - this->_linkedCells.getBoxMin() +
            std::array<double, 3>{interactionLength, interactionLength, interactionLength} * 2.;

        const auto listLengthEstimate = VerletListsCellsHelpers::estimateListLength(
            this->_linkedCells.getNumberOfParticles(IteratorBehavior::ownedOrHalo), boxSizeWithHalo, interactionLength,
            1.3);

        neighborLists.resize(cells.size());

        // Helper function to insert a pointer into a list of the base cell.
        // It considers the cases that neither particle is in the base cell and in that case
        // finds or creates the appropriate list
        auto insert = [&](auto &p1, auto i, auto &p2, auto cellIndex1, auto cellIndexBase, auto &currentCellsLists) {
          // Easy case: cell1 is the base cell
          if (cellIndexBase == cellIndex1) {
            currentCellsLists[i].second.push_back(&p2);
          } else {
            // Otherwise, check if the base cell already has a list for p1
            auto iter = std::find_if(currentCellsLists.begin(), currentCellsLists.end(), [&](const auto &pair) {
              const auto &[particlePtr, list] = pair;
              return particlePtr == &p1;
            });
            // If yes, append p2 to it.
            if (iter != currentCellsLists.end()) {
              iter->second.push_back(&p2);
            } else {
              // If no, create one, reserve space and emplace p2
              currentCellsLists.emplace_back(&p1, std::vector<Particle *>{});
              currentCellsLists.back().second.reserve(listLengthEstimate);
              currentCellsLists.back().second.push_back(&p2);
            }
          }
        };

        const auto estimateListSize = [&](size_t baseCellIndex) {
          size_t estimate = 0;
          std::map<int, double> offsetFactors{};
          std::map<int, double> offsetFactorsNoN3{};
          for (const auto [offsetA, offsetB, factor] : offsets) {
            offsetFactors[offsetA] = std::max(offsetFactors[offsetA], factor);
            offsetFactorsNoN3[offsetB] = std::max(offsetFactors[offsetB], factor);
          }
          for (const auto &[offset, factor] : offsetFactors) {
            estimate += cells[baseCellIndex + offset].size() * factor;
          }
          if (not this->_verletBuiltNewton3) {
            for (const auto &[offset, factor] : offsetFactorsNoN3) {
              estimate += cells[baseCellIndex + offset].size() * factor;
            }
          }
          return estimate;
        };

        // Go over all cells except the very last layer
      AUTOPAS_OPENMP(parallel for collapse(3))
      for (int z = 0; z < cellsPerDim[2] - 1; ++z) {
        for (int y = 0; y < cellsPerDim[1] - 1; ++y) {
          for (int x = 0; x < cellsPerDim[0] - 1; ++x) {
            const auto cellIndexBase = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, cellsPerDim);
            auto &baseCell = cells[cellIndexBase];
            auto &baseCellsLists = neighborLists[cellIndexBase];
            // Allocate space for ptr-list pairs for this cell
            baseCellsLists.reserve(estimateListSize(cellIndexBase));
            for (size_t i = 0; i < baseCell.size(); ++i) {
              baseCellsLists.emplace_back(&baseCell[i], std::vector<Particle *>{});
              baseCellsLists.back().second.reserve(listLengthEstimate);
            }

            // Build c08 lists according to predefined cell pairs
            for (const auto &[offset1, offset2, _] : offsets) {
              const auto cellIndex1 = cellIndexBase + offset1;
              const auto cellIndex2 = cellIndexBase + offset2;

              // Skip if both cells only contain halos
              if (not(cells[cellIndex1].getPossibleParticleOwnerships() == OwnershipState::owned) and
                  not(cells[cellIndex2].getPossibleParticleOwnerships() == OwnershipState::owned)) {
                continue;
              }

              // Go over all particle pairs in the two cells
              for (size_t i = 0; i < cells[cellIndex1].size(); ++i) {
                auto &p1 = cells[cellIndex1][i];
                for (size_t j = (cellIndex1 == cellIndex2) ? i + 1 : 0; j < cells[cellIndex2].size(); ++j) {
                  auto &p2 = cells[cellIndex2][j];
                  if (&p1 == &p2 or p1.isDummy() or p2.isDummy()) {
                    continue;
                  }

                  const auto distVec = p2.getR() - p1.getR();
                  const auto distSquared = utils::ArrayMath::dot(distVec, distVec);
                  if (distSquared < interactionLengthSquared) {
                    insert(p1, i, p2, cellIndex1, cellIndexBase, baseCellsLists);

                    if (not this->_verletBuiltNewton3) {
                      insert(p2, j, p1, cellIndex2, cellIndexBase, baseCellsLists);
                    }
                  }
                }
              }
            }
          }
        }
      }
      } else {
#endif
        _neighborList.buildAoSNeighborList(this->_linkedCells, this->_verletBuiltNewton3, this->getCutoff(),
                                           this->getVerletSkin(), this->getInteractionLength(),
                                           traversal->getTraversalType(), _dataLayoutDuringListRebuild);
#ifdef BUILD_WITHOUT_FUNCTOR
      }
#endif
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
