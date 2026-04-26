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
class HGTraversalBase {
 public:
  /**
   * Type of particles stored in the traversed cells.
   */
  using ParticleType = typename ParticleCell_T::ParticleType;

  /**
   * Constructor.
   * @param numLevels Number of levels in the hierarchical grid.
   * @param dataLayout Data layout used by this traversal.
   * @param useNewton3 Whether Newton3 optimization is enabled.
   */
  explicit HGTraversalBase(size_t numLevels, DataLayoutOption dataLayout, bool useNewton3)
      : _numLevels(numLevels), _levels(nullptr), _skin(0) {}

  /**
   * Store HGrid data
   * @param levels LinkedCells for each HGrid level
   * @param maxCutoffPerLevel maxCutoffPerLevel of each HGrid level
   * @param skin verlet skin of HGrid container
   */
  void setLevels(std::vector<std::unique_ptr<LinkedCells<ParticleType>>> *levels,
                 std::vector<double> &maxCutoffPerLevel, double skin) {
    _numLevels = maxCutoffPerLevel.size();
    _maxCutoffPerLevel = maxCutoffPerLevel;
    _levels = levels;
    _skin = skin;
  }

 protected:
  /**
   * Number of hierarchy levels.
   */
  size_t _numLevels;
  /**
   * Pointer to all LinkedCells levels of the hierarchical grid.
   */
  std::vector<std::unique_ptr<LinkedCells<ParticleType>>> *_levels;
  /**
   * Maximum cutoff per hierarchy level.
   */
  std::vector<double> _maxCutoffPerLevel;
  /**
   * Configured verlet skin of the container.
   */
  double _skin;
  /**
   * Cutoff ratio threshold for enabling cell distance checks.
   * At this ratio of max Cutoffs between two levels, distance checks are employed to sort out unnecessary interactions.
   * Particles can do particle to level specific instead of level to level specific in the future.
   */
  const double distCheckRatio = 0.2;

  /**
   * Type alias for cell blocks used by hierarchical levels.
   */
  using CellBlock = internal::CellBlock3D<FullParticleCell<ParticleType>>;

  /**
   * Get the interaction length for a given level pair.
   * This is simply the average of the maximum possible cutoff of the two levels plus the verlet skin.
   * @param lowerLevel The lower level index.
   * @param upperLevel The upper level index.
   * @return The interaction length for the given level pair.
   */
  [[nodiscard]] double getInteractionLength(const size_t lowerLevel, const size_t upperLevel) const {
    const double cutoff = (this->_maxCutoffPerLevel[upperLevel] + this->_maxCutoffPerLevel[lowerLevel]) / 2;
    return cutoff + _skin;
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

  // HG traversals only make sense if there are multiple levels, otherwise just use LinkedCells
  [[nodiscard]] bool isApplicable() const { return _numLevels > 1; }
};
}  // namespace autopas
