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
      : TraversalInterface(dataLayout, useNewton3), _numLevels(0), _levels(nullptr), _skin(0), _maxDisplacement(0) {}

  /**
   * Store HGrid data
   * @param levels LinkedCells for each HGrid level
   * @param cutoffs cutoffs of each HGrid level
   * @param skin verlet skin of HGrid container
   * @param maxDisplacement maximum displacement of any particle in the container, ignored if dynamic containers is not used
   */
  void setLevels(std::vector<std::unique_ptr<LinkedCells<Particle>>> *levels, std::vector<double> &cutoffs,
                 double skin, double maxDisplacement) {
    _numLevels = cutoffs.size();
    _cutoffs = cutoffs;
    _levels = levels;
    _skin = skin;
    _maxDisplacement = maxDisplacement;
  }

 protected:
  size_t _numLevels;
  std::vector<std::unique_ptr<LinkedCells<Particle>>> *_levels;
  std::vector<double> _cutoffs;
  double _skin;
  double _maxDisplacement;

  using CBParticleCell = FullParticleCell<Particle>;

  template <typename T>
  T computeIndex(T x, T y, T z, T L, T W) {
    int term = y * L + ((y % 2 == 0) ? x : (L - 1 - x));
    if (z % 2 == 0) {
      return z * L * W + term;
    } else {
      return (z + 1) * L * W - 1 - term;
    }
  }

  /**
  * Convert a 1d index to a 3d index in a snaky pattern so that between consecutive indexes,
  * the manhattan distance of the 3D coordinates is 1.
  * @tparam T Type of the indices.
  * @param ind The 1d index.
  * @param dims The total dimensions of the index space.
  * @return The 3d index.
  */
  template <typename T>
  constexpr std::vector<std::array<T, 3>> oneToThreeDForStart(const std::array<T, 3> &dims) {
    static_assert(std::is_integral_v<T>, "oneToThreeD requires integral types");
    std::vector<std::array<T, 3> > coords(dims[0] * dims[1] * dims[2]);
    for (T z = 0; z < dims[2]; ++z) {
      for (T y = 0; y < dims[1]; ++y) {
        for (T x = 0; x < dims[0]; ++x) {
          int idx = computeIndex(x, y, z, dims[0], dims[1]);
          coords[idx] = {x, y, z};
        }
      }
    }
    return coords;
  }


  /**
   *
   * @param upperCB first CellBlock3D
   * @param upperCell 3d cell index of cell belonging to first CellBlock3D
   * @param lowerCB second CellBlock3D
   * @param lowerCell 3d cell index of cell belonging to second CellBlock3D
   * @return minimum distance between two cells squared
   */
  double getMinDistBetweenCellsSquared(const internal::CellBlock3D<CBParticleCell> &upperCB,
                                       const std::array<size_t, 3> &upperCell,
                                       const internal::CellBlock3D<CBParticleCell> &lowerCB,
                                       const std::array<size_t, 3> &lowerCell) {
    const auto posX = upperCB.getCellBoundingBox(upperCell);
    const auto posY = lowerCB.getCellBoundingBox(lowerCell);
    double totalDist = 0;
    for (size_t i = 0; i < 3; ++i) {
      const auto dist = std::max(0.0, std::max(posX.first[i] - posY.second[i], posY.first[i] - posX.second[i]));
      totalDist += dist * dist;
    }
    return totalDist;
  }

  /**
   *
   * @param CB cell block that contains the cell
   * @param cellIndex 3d cell index of cell belonging to first CellBlock3D
   * @param point
   * @return minimum distance between a cell and a point squared
   */
  double getMinDistBetweenCellAndPointSquared(const internal::CellBlock3D<CBParticleCell> &CB,
                                              const std::array<size_t, 3> &cellIndex,
                                              const std::array<double, 3> &point) {
    const auto posX = CB.getCellBoundingBox(cellIndex);
    double totalDist = 0;
    for (size_t i = 0; i < 3; ++i) {
      const auto dist = std::max(0.0, std::max(posX.first[i] - point[i], point[i] - posX.second[i]));
      totalDist += dist * dist;
    }
    return totalDist;
  }

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
    // TODO: _overlap will be smaller than needed, because smaller levels have big halo regions compared to their actual
    // interactionlength
    // TODO: so cells past _overlap can also be halo cells, resulting in unnecessary iteration of those halo cells
    const TraversalSelectorInfo ret{temp.cellsPerDim, _cutoffs[level] + _maxDisplacement * 2, temp.cellLength, temp.clusterSize};
    return ret;
  }
};
}  // namespace autopas
