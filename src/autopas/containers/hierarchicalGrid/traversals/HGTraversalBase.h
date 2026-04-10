/**
 * @file HGTraversalBase.h
 * @author atacann
 * @date 09.12.2024
 */

#pragma once

#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {
/**
 * Class for common operations of HGridTraversal classes. It should be a parent class for all HGrid traversals.
 * @tparam ParticleCell_T type of Particle cell
 */
template <class ParticleCell_T>
class HGTraversalBase : public TraversalInterface {
 public:
  using ParticleType = typename ParticleCell_T::ParticleType;

  explicit HGTraversalBase(DataLayoutOption dataLayout, bool useNewton3)
      : TraversalInterface(dataLayout, useNewton3),
        _numLevels(0),
        _levels(nullptr),
        _skin(0),
        _stepsSinceLastRebuild(0),
        _rebuildFrequency(1) {}

  /**
   * Store HGrid data
   * @param levels LinkedCells for each HGrid level
   * @param maxCutoffPerLevel maxCutoffPerLevel of each HGrid level
   * @param skin verlet skin of HGrid container
   * @param stepsSinceLastRebuild number of time-steps since last rebuild
   * @param rebuildFrequency frequency of rebuild
   */
  void setLevels(std::vector<std::unique_ptr<LinkedCells<ParticleType>>> *levels,
                 std::vector<double> &maxCutoffPerLevel, double skin, unsigned int stepsSinceLastRebuild,
                 const unsigned int rebuildFrequency) {
    _numLevels = maxCutoffPerLevel.size();
    _maxCutoffPerLevel = maxCutoffPerLevel;
    _levels = levels;
    _skin = skin;
    _stepsSinceLastRebuild = stepsSinceLastRebuild;
    _rebuildFrequency = rebuildFrequency;
  }

 protected:
  size_t _numLevels;
  std::vector<std::unique_ptr<LinkedCells<ParticleType>>> *_levels;
  std::vector<double> _maxCutoffPerLevel;
  double _skin;
  unsigned int _stepsSinceLastRebuild;
  unsigned int _rebuildFrequency;
  // At this ratio of max Cutoffs between two levels, distance checks are employed to sort out unnecessary interactions.
  // This is needed as the number of interactions can increase drastically if the cell size difference is too large.
  // particles can do particle to level specific instead of level to level specific in the future
  const double distCheckRatio = 0.2;
  // Intralevel traversals used for each level for this iteration
  std::vector<std::unique_ptr<TraversalInterface>> _traversals;

  using CellBlock = internal::CellBlock3D<FullParticleCell<ParticleType>>;

  /**
   * Generate a new Traversal from the given data, needed as each level of HGrid has different cell sizes
   * @param level which HGrid level to generate a traversal for
   * @return A new traversal that is applicable to a specific LinkedCells level
   */
  virtual std::unique_ptr<TraversalInterface> generateNewTraversal(const size_t level) = 0;

  /**
   * Get the interaction length for a given level pair.
   * This is simply the average of the maximum possible cutoff of the two levels plus the verlet skin.
   * @param lowerLevel The lower level index.
   * @param upperLevel The upper level index.
   * @return The interaction length for the given level pair.
   */
  [[nodiscard]] double getInteractionLength(const size_t lowerLevel, const size_t upperLevel) const {
    const double cutoff = (this->_maxCutoffPerLevel[upperLevel] + this->_maxCutoffPerLevel[lowerLevel]) / 2;
    return cutoff + currentSkin();
  }

  /**
   * Check if the ratio of level grid sizes is small enough to that the distance between each cell should be checked
   * before applying functors.
   * @param upperLevel The upper level index.
   * @param lowerLevel The lower level index.
   * @return True if the distance is small enough, false otherwise.
   */
  [[nodiscard]] bool checkDistance(const size_t upperLevel, const size_t lowerLevel) const {
    // size of the lower level cell divided by interactionLength
    // We don't use cellLength because cellSizeFactor might have influenced it
    return (_maxCutoffPerLevel[lowerLevel] + _skin) / getInteractionLength(lowerLevel, upperLevel) <= distCheckRatio;
  }

  [[nodiscard]] double currentSkin() const {
    // We only need to check at most distance cutoff + max displacement of any particle in the container if dynamic
    // containers are used.
    // NOTE: if in the future, if Hgrid will be used as a base container to a verlet list, interactionLength should be
    // always cutoff + _skin.
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
    return this->_skin; //std::min(this->_maxDisplacement * 2 + 1e-9, this->_skin);
#else
    return std::min((this->_skin / _rebuildFrequency) * _stepsSinceLastRebuild + 1e-9, this->_skin);
#endif
  }

  /**
   * Compute the needed stride to avoid race conditions.
   * For SoA consider that two cells interacting on the same lower level cell is problematic, as two threads can operate
   * on the same buffer.
   * @param level The hierarchy level to compute the stride for.
   * @return The stride for each dimension.
   */
  std::array<size_t, 3> computeStride(const size_t level, const bool topDown = true) {
    const std::array<double, 3> levelLength = this->_levels->at(level)->getTraversalSelectorInfo().cellLength;
    std::array<size_t, 3> stride{1, 1, 1};
    for (size_t otherLevel = 0; otherLevel < this->_numLevels; ++otherLevel) {
      if (otherLevel == level) {
        continue;
      }
      if (this->_useNewton3 and ((otherLevel >= level and topDown) or (otherLevel <= level and !topDown))) {
        continue;
      }
      const double interactionLength = this->getInteractionLength(level, otherLevel);
      const std::array<double, 3> otherLevelLength =
          this->_levels->at(otherLevel)->getTraversalSelectorInfo().cellLength;
      std::array<size_t, 3> tempStride{};
      for (size_t i = 0; i < 3; i++) {
        if (this->_useNewton3) {
          // find out the stride so that cells we check on lowerLevel do not intersect
          if (this->_dataLayout == DataLayoutOption::soa) {
            tempStride[i] = 1 + static_cast<size_t>(std::ceil(std::ceil(interactionLength / otherLevelLength[i]) * 2 *
                                                              otherLevelLength[i] / levelLength[i]));
            // the stride calculation below, for AoS, can result in less colors, but two different threads can operate
            // on SoA buffer of the same cell at the same time. They won't update the same value at the same time, but
            // causes race conditions in SoA. (Inbetween loading data into vectors (avx etc.) -> functor calcs -> store
            // again).
          } else {
            tempStride[i] = 1 + static_cast<size_t>(std::ceil(interactionLength * 2 / levelLength[i]));
          }
        } else {
          // do c01 traversal if newton3 is disabled
          tempStride[i] = 1;
        }
        stride[i] = std::max(stride[i], tempStride[i]);
      }
    }
    return stride;
  }

  /**
   *
   * @param CB cell block that contains the cell
   * @param cellIndex 1d cell index of cell belonging to first CellBlock3D
   * @param point
   * @return minimum distance between a cell and a point squared
   */
  double getMinDistBetweenCellAndPointSquared(const CellBlock &CB, const size_t cellIndex,
                                              const std::array<double, 3> &point) {
    const auto cellPos = CB.getCellBoundingBox(cellIndex);
    double totalDist = 0;
    for (size_t i = 0; i < 3; ++i) {
      const auto dist = std::max(0.0, std::max(cellPos.first[i] - point[i], point[i] - cellPos.second[i]));
      totalDist += dist * dist;
    }
    return totalDist;
  }

  /**
   * Returns traversalSelectorInfo for the specific level
   * @param level Which level to get info from
   * @return TraversalSelectorInfo of level
   */
  [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo(int level) const {
    if (level < 0 || level >= _numLevels) {
      autopas::utils::ExceptionHandler::exception("Hierarchical grid level out of range: {}", level);
    }
    // return traversal info of hierarchy with biggest interactionLength
    const auto temp = _levels->at(level)->getTraversalSelectorInfo();
    // adjust interactionLength to actual value
    // _overlap will be smaller than needed, because smaller levels have big halo regions compared to their actual
    // interactionlength so cells past _overlap can also be halo cells, resulting in unnecessary iteration of those halo
    // cells
    const TraversalSelectorInfo ret{temp.cellsPerDim, this->getInteractionLength(level, level), temp.cellLength,
                                    temp.clusterSize};
    return ret;
  }

  /**
   * Compute intra-level interactions for all levels.
   */
  void computeIntraLevelInteractions() {
    for (size_t level = 0; level < this->_numLevels; level++) {
      // We do not simply call computeInteractions() here as we want to store SoA after inter-level traversals are
      // computed. They will be loaded in HGTraversalBase::endTraversal().
      this->_traversals[level]->traverseParticles();
    }
  }

  void initTraversal() override {
    this->_traversals.resize(this->_numLevels);
    for (size_t level = 0; level < this->_numLevels; ++level) {
      this->_traversals[level] = this->generateNewTraversal(level);
      TraversalInterface *intraLevelTraversal = this->_traversals[level].get();
      this->_levels->at(level)->prepareTraversal(intraLevelTraversal);
      this->_traversals[level]->initTraversal();
    }
  }

  void endTraversal() override {
    // stores SoA if SoA is used
    for (size_t level = 0; level < this->_numLevels; level++) {
      this->_traversals[level]->endTraversal();
    }
  }

  // HG traversals only make sense if there are multiple levels, otherwise just use LinkedCells
  [[nodiscard]] bool isApplicable() const override { return not _maxCutoffPerLevel.empty(); }

  /**
   * Traverses a single upper-level cell and lower-level cells that are in the interaction range using AoS.
   * @tparam Functor_T type of the functor
   * @param lowerCB lower level cell block
   * @param upperCB upper level cell block
   * @param upperCellCoords 3d cell index of cell belonging to the upper level
   * @param functor functor to apply
   * @param lowerLevel lower level index
   * @param lowerBound inclusive lower coordinate bound for lower-level owned cells when only interacting with owned
   * particles
   * @param upperBound inclusive upper coordinate bound for lower-level owned cells when only interacting with owned
   * particles
   */
  template <class Functor_T>
  void AoSTraversal(const CellBlock &lowerCB, const CellBlock &upperCB, const std::array<size_t, 3> upperCellCoords,
                    Functor_T *functor, size_t lowerLevel, size_t upperLevel, const std::array<size_t, 3> &lowerBound,
                    const std::array<size_t, 3> &upperBound) {
    // can also use sorted cell optimization? need to be careful to sort only once per upper cell, to not sort for each
    // particle in the upper cell
    using namespace autopas::utils::ArrayMath::literals;
    using utils::ArrayUtils::operator<<;

    const bool enableCellDistanceChecking = checkDistance(upperLevel, lowerLevel);
    const auto &lowerLevelDims = lowerCB.getCellsPerDimensionWithHalo();
    auto &upperCell = upperCB.getCell(upperCellCoords);
    if (upperCell.isEmpty()) {
      return;
    }
    // skip if cell is a halo cell and newton3 is disabled
    auto isHalo = upperCell.getPossibleParticleOwnerships() == OwnershipState::halo;
    if (isHalo && !this->_useNewton3) {
      return;
    }
    const double lowerInteractionLength = currentSkin() + _maxCutoffPerLevel[lowerLevel] / 2;
    const std::array<double, 3> lowerDir{lowerInteractionLength, lowerInteractionLength, lowerInteractionLength};
    // variable to determine if we are only interested in owned particles in the lower level
    const bool containToOwnedOnly = isHalo && this->_useNewton3;
    for (auto p1Ptr = upperCell.begin(); p1Ptr != upperCell.end(); ++p1Ptr) {
      const std::array<double, 3> &pos = p1Ptr->getR();
      const double radius = p1Ptr->getSize() / 2;
      const std::array<double, 3> interactionLength = lowerDir + radius;
      const double interactionLengthSquared = interactionLength[0] * interactionLength[0];
      auto startIndex3D = lowerCB.get3DIndexOfPosition(pos - interactionLength);
      auto stopIndex3D = lowerCB.get3DIndexOfPosition(pos + interactionLength);
      // skip halo cells if we need to consider only owned particles
      if (containToOwnedOnly) {
        startIndex3D = utils::ArrayMath::max(startIndex3D, lowerBound);
        stopIndex3D = utils::ArrayMath::min(stopIndex3D, upperBound);
      }
      const auto particleCellIndex3D =
          enableCellDistanceChecking ? lowerCB.get3DIndexOfPosition(pos) : std::array{0ul, 0ul, 0ul};
      for (size_t zl = startIndex3D[2]; zl <= stopIndex3D[2]; ++zl) {
        for (size_t yl = startIndex3D[1]; yl <= stopIndex3D[1]; ++yl) {
          auto cellIndex1D = autopas::utils::ThreeDimensionalMapping::threeToOneD(
              {static_cast<size_t>(startIndex3D[0]), yl, zl}, lowerLevelDims);
          for (size_t xl = startIndex3D[0]; xl <= stopIndex3D[0]; ++xl, ++cellIndex1D) {
            auto &lowerCell = lowerCB.getCell(cellIndex1D);
            if (lowerCell.isEmpty()) {
              continue;
            }
            // skip if cell is farther than interactionLength
            if (enableCellDistanceChecking and
                this->getMinDistBetweenCellAndPointSquared(lowerCB, cellIndex1D, pos) > interactionLengthSquared) {
              continue;
            }
            for (auto &p : lowerCell) {
              functor->AoSFunctor(*p1Ptr, p, this->_useNewton3);
            }
          }
        }
      }
    }
  }

  /**
   * Traverses a single upper level cell and a lower level cells that are in the interaction range using SoA.
   * SoA is calculated between a single upper level particle and a full lower level cell.
   * @tparam Functor_T type of the functor
   * @param lowerCB lower level cell block
   * @param upperCB upper level cell block
   * @param upperCellCoords 3d cell index of cell belonging to the upper level
   * @param functor functor to apply
   * @param lowerLevel lower level index
   * @param lowerBound lower bound of the coordinates of lower level cells that contain owned particles
   * @param upperBound upper bound of the coordinates of lower level cells that contain owned particles
   * less than the interaction length, otherwise skips the lower level cell
   */
  template <class Functor_T>
  void SoATraversalParticleToCell(const CellBlock &lowerCB, const CellBlock &upperCB,
                                  const std::array<size_t, 3> upperCellCoords, Functor_T *functor, size_t lowerLevel,
                                  size_t upperLevel, const std::array<size_t, 3> &lowerBound,
                                  const std::array<size_t, 3> &upperBound) {
    using namespace autopas::utils::ArrayMath::literals;
    const bool enableCellDistanceChecking = checkDistance(upperLevel, lowerLevel);
    const auto &lowerLevelDims = lowerCB.getCellsPerDimensionWithHalo();
    auto &upperCell = upperCB.getCell(upperCellCoords);
    if (upperCell.isEmpty()) {
      return;
    }
    // skip if cell is a halo cell and newton3 is disabled
    auto isHalo = upperCell.getPossibleParticleOwnerships() == OwnershipState::halo;
    if (isHalo && !this->_useNewton3) {
      return;
    }
    // variable to determine if we are only interested in owned particles in the lower level
    const bool containToOwnedOnly = isHalo && this->_useNewton3;
    auto &soa = upperCell._particleSoABuffer;
    const auto *const __restrict xptr = soa.template begin<ParticleType::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<ParticleType::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<ParticleType::AttributeNames::posZ>();

    const double lowerInteractionLength = currentSkin() + _maxCutoffPerLevel[lowerLevel] / 2;
    const std::array<double, 3> lowerDir{lowerInteractionLength, lowerInteractionLength, lowerInteractionLength};

    for (int idx = 0; idx < upperCell.size(); ++idx) {
      const std::array<double, 3> pos = {xptr[idx], yptr[idx], zptr[idx]};
      const double radius = upperCell[idx].getSize() / 2;
      const std::array<double, 3> interactionLength = lowerDir + radius;
      const double interactionLengthSquared = interactionLength[0] * interactionLength[0];
      auto soaSingleParticle = soa.constructView(idx, idx + 1);
      auto startIndex3D = lowerCB.get3DIndexOfPosition(pos - interactionLength);
      auto stopIndex3D = lowerCB.get3DIndexOfPosition(pos + interactionLength);
      if (containToOwnedOnly) {
        startIndex3D = utils::ArrayMath::max(startIndex3D, lowerBound);
        stopIndex3D = utils::ArrayMath::min(stopIndex3D, upperBound);
      }

      const auto particleCellIndex3D =
          enableCellDistanceChecking ? lowerCB.get3DIndexOfPosition(pos) : std::array{0ul, 0ul, 0ul};
      for (size_t zl = startIndex3D[2]; zl <= stopIndex3D[2]; ++zl) {
        for (size_t yl = startIndex3D[1]; yl <= stopIndex3D[1]; ++yl) {
          auto cellIndex1D = autopas::utils::ThreeDimensionalMapping::threeToOneD(
              {static_cast<size_t>(startIndex3D[0]), yl, zl}, lowerLevelDims);
          for (size_t xl = startIndex3D[0]; xl <= stopIndex3D[0]; ++xl, ++cellIndex1D) {
            auto &lowerCell = lowerCB.getCell(cellIndex1D);
            if (lowerCell.isEmpty()) {
              continue;
            }
            // skip if cell is farther than interactionLength
            if (enableCellDistanceChecking and
                this->getMinDistBetweenCellAndPointSquared(lowerCB, cellIndex1D, pos) > interactionLengthSquared) {
              continue;
            }
            // 1 to n SoAFunctorPair
            functor->SoAFunctorPair(soaSingleParticle, lowerCell._particleSoABuffer, this->_useNewton3);
          }
        }
      }
    }
  }
};
}  // namespace autopas
