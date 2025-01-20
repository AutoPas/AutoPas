/**
 * @file HGTraversalBase.h
 * @author atacann
 * @date 09.12.2024
 */

#pragma once

#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {
/**
 * Class for common operations of HGridTraversal classes
 * Also compatible for use in HierarchicalGrid.h as it does not have Functor as a tparam.
 * @tparam ParticleCell type of Particle cell
 */
template <class ParticleCell>
class HGTraversalBase : public TraversalInterface {
 public:
  using Particle = typename ParticleCell::ParticleType;

  explicit HGTraversalBase(DataLayoutOption dataLayout, bool useNewton3)
      : TraversalInterface(dataLayout, useNewton3), _numLevels(0), _levels(nullptr), _skin(0) {}

  /**
   * Store HGrid data
   * @param levels LinkedCells for each HGrid level
   * @param cutoffs cutoffs of each HGrid level
   * @param skin verlet skin of HGrid container
   */
  void setLevels(std::vector<std::unique_ptr<LinkedCells<Particle>>>* levels, std::vector<double> &cutoffs,
                 double skin) {
    _numLevels = cutoffs.size();
    _cutoffs = cutoffs;
    _levels = levels;
    _skin = skin;
  }

 protected:
  size_t _numLevels;
  std::vector<std::unique_ptr<LinkedCells<Particle>>> *_levels;
  std::vector<double> _cutoffs;
  double _skin;

  // double getMinDistBetweenCellsSquared(const internal::CellBlock3D<ParticleCell> &upperCB, std::array<size_t, 3> &upperCell,
  //                                      const internal::CellBlock3D<ParticleCell> &lowerCB, std::array<size_t, 3> &lowerCell) {
  //   const auto posX = upperCB.getCellBoundingBox(upperCell);
  //   const auto posY = lowerCB.getCellBoundingBox(lowerCell);
  //   double sum = 0;
  //   return sum;
  // }

  /**
   *
   * @tparam T type of elements in arrays
   * @tparam SIZE size of arrays
   * @param a first array
   * @param b second array
   * @return max of two arrays element-wise
   */
  template <class T, std::size_t SIZE>
  constexpr std::array<T, SIZE> getMax(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
    std::array<T, SIZE> ret{};
    for (size_t i = 0; i < SIZE; ++i) {
      ret[i] = std::max(a[i], b[i]);
    }
    return ret;
  }

  /**
   *
   * @tparam T type of elements in arrays
   * @tparam SIZE size of arrays
   * @param a first array
   * @param b second array
   * @return min of two arrays element-wise
   */
  template <class T, std::size_t SIZE>
  constexpr std::array<T, SIZE> getMin(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
    std::array<T, SIZE> ret{};
    for (size_t i = 0; i < SIZE; ++i) {
      ret[i] = std::min(a[i], b[i]);
    }
    return ret;
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
    // TODO: _overlap will be smaller than needed, because smaller levels have big halo regions compared to their actual interactionlength
    // TODO: so cells past _overlap can also be halo cells, resulting in unnecessary iteration of those halo cells
    const TraversalSelectorInfo ret{temp.cellsPerDim, _cutoffs[level] + _skin, temp.cellLength, temp.clusterSize};
    return ret;
  }
};
}  // namespace autopas
