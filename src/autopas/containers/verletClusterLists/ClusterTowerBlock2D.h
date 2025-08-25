/**
 * @file ClusterTowerBlock2D.h
 * @author F. Gratl
 * @date 14.09.23
 */

#pragma once

#include <array>
#include <vector>

#include "autopas/containers/CellBorderAndFlagManager.h"
#include "autopas/containers/verletClusterLists/ClusterTower.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas::internal {

/**
 * Class to manage the grid of towers for the Verlet Cluster Lists container.
 * This also includes the geometric shape of the container.
 *
 * @tparam Particle_T
 */
template <class Particle_T>
class ClusterTowerBlock2D : public CellBorderAndFlagManager {
 public:
  /**
   * Constructor.
   *
   * @note This constructor does not initialize the tower grid! For this, the number of particles is required.
   *
   * @param boxMin
   * @param boxMax
   * @param interactionLength
   */
  ClusterTowerBlock2D(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                      double interactionLength)
      : _interactionLength(interactionLength),
        _boxMin(boxMin),
        _boxMax(boxMax),
        _haloBoxMin{utils::ArrayMath::subScalar(boxMin, interactionLength)},
        _haloBoxMax{utils::ArrayMath::addScalar(boxMax, interactionLength)} {}
  /**
   * Start iterator over towers.
   * @return
   */
  typename std::vector<ClusterTower<Particle_T>>::iterator begin() { return _towers.begin(); }
  /**
   * Start iterator over towers.
   * @return
   */
  typename std::vector<ClusterTower<Particle_T>>::const_iterator begin() const { return _towers.begin(); }
  /**
   * End iterator over towers.
   * @return
   */
  typename std::vector<ClusterTower<Particle_T>>::iterator end() { return _towers.end(); }
  /**
   * End iterator over towers.
   * @return
   */
  typename std::vector<ClusterTower<Particle_T>>::const_iterator end() const { return _towers.end(); }

  /**
   * Access operator for towers.
   * @param i 1D tower index
   * @return
   */
  ClusterTower<Particle_T> &operator[](size_t i) { return _towers[i]; }

  /**
   * Access operator for towers.
   * @param i 1D tower index
   * @return
   */
  const ClusterTower<Particle_T> &operator[](size_t i) const { return _towers[i]; }

  /**
   * Return the number of towers.
   * @return
   */
  size_t size() const { return _towers.size(); }

  /**
   * Resize the internal grid and storage.
   * @param towerSideLength
   * @param towersPerDim
   */
  void resize(const std::array<double, 2> &towerSideLength, const std::array<size_t, 2> &towersPerDim) {
    using namespace autopas::utils::ArrayMath::literals;
    using autopas::utils::ArrayMath::ceil;
    using autopas::utils::ArrayUtils::static_cast_copy_array;
    _towerSideLength = towerSideLength;
    _towerSideLengthReciprocal = 1. / towerSideLength;
    const auto numTowersPerInteractionLength2D =
        static_cast_copy_array<int>(ceil(_interactionLength / _towerSideLength));
    // Sanity check. This is not necessarily a problem with the method but potentially with the way the code is written.
    if (numTowersPerInteractionLength2D[0] != numTowersPerInteractionLength2D[1]) {
      AutoPasLog(WARN,
                 "Number of towers per interaction length differs in X vs Y direction! "
                 "The container is built on the assumption that they are the same, "
                 "so from here on the behavior might be undefined. Values: {}",
                 utils::ArrayUtils::to_string(numTowersPerInteractionLength2D));
    }
    _numTowersPerInteractionLength =
        *std::max_element(numTowersPerInteractionLength2D.begin(), numTowersPerInteractionLength2D.end());
    _towersPerDim = towersPerDim;
    _towers.resize(_towersPerDim[0] * _towersPerDim[1]);

    // converting indices relies on _towersPerDim to already have the new values.
    _firstOwnedTowerIndex = towerIndex2DTo1D(_numTowersPerInteractionLength, _numTowersPerInteractionLength);
    _lastOwnedTowerIndex = towerIndex2DTo1D(_towersPerDim[0] - _numTowersPerInteractionLength,
                                            _towersPerDim[1] - _numTowersPerInteractionLength);
  }

  /**
   * Indicator if the block contains any towers.
   * @return
   */
  bool empty() const { return _towers.empty(); }

  /**
   * Construct a new tower at the end of the storage.
   * @param clusterSize
   */
  void addTower(size_t clusterSize) { _towers.emplace_back(clusterSize); }

  /**
   * Estimates the optimal 2D tower grid side length.
   *
   * sideLength = (clusterSize / (numParticles / volume))^(1/3)
   *
   * The side length is then adjusted such that towers exactly fill the (non-halo) box by increasing
   * the side length such that the number of towers does not change.
   *
   * @param numParticles The number of particles in the container.
   * @param clusterSize the number of particles per cluster.
   * @return tuple{gridSideLength, numberOfTowers}.
   */
  [[nodiscard]] std::tuple<std::array<double, 2>, std::array<size_t, 2>> estimateOptimalGridSideLength(
      size_t numParticles, size_t clusterSize) const {
    using namespace autopas::utils::ArrayMath::literals;
    using utils::ArrayMath::ceil;
    using utils::ArrayUtils::static_cast_copy_array;

    // estimate the grid size for the actual box, then add halo layers as necessary
    const std::array<double, 3> boxSize = _boxMax - _boxMin;
    const std::array<double, 2> boxSize2D{boxSize[0], boxSize[1]};
    if (numParticles == 0) {
      // We always build at least 3 towers per direction (one owned, two halo)
      // This is needed for sliced traversals to work correctly.
      return {boxSize2D, {3, 3}};
    }

    const double volume = boxSize[0] * boxSize[1] * boxSize[2];
    // We use double here to avoid back and forth casting
    // estimate particle density
    const double density = static_cast<double>(numParticles) / volume;
    const auto optimalSideLength = std::cbrt(static_cast<double>(clusterSize) / density);
    // make sure the towers fill the domain exactly
    const auto numTowersOwned(ceil(boxSize2D / optimalSideLength));
    const auto towerSideLengthNew = boxSize2D / numTowersOwned;
    // pad a halo layer of towers in both directions, filling one interaction length each.
    const auto numTowers = numTowersOwned + (ceil(_interactionLength / towerSideLengthNew) * 2.);
    // Do not resize the number of towers here!
    // If we get less towers we need to save the particles in the towers we throw away.
    return {towerSideLengthNew, static_cast_copy_array<size_t>(numTowers)};
  }

  /**
   * Calculates the low and high corner of a tower given by its index.
   *
   * @note If towers are not built yet the corners of the full container are returned.
   *
   * @param index1D The tower's index in _towerBlock.
   * @return tuple<lowCorner, highCorner>
   */
  [[nodiscard]] std::tuple<std::array<double, 3>, std::array<double, 3>> getTowerBoundingBox(size_t index1D) const {
    // case: towers are not built yet.
    if (_towersPerDim[0] == 0) {
      return {_boxMin, _boxMax};
    }
    return getTowerBoundingBox(towerIndex1DTo2D(index1D));
  }

  /**
   * Calculates the low and high corner of a tower given by its 2D grid index.
   *
   * @param index2D
   * @return
   */
  std::tuple<std::array<double, 3>, std::array<double, 3>> getTowerBoundingBox(
      const std::array<size_t, 2> &index2D) const {
    // case: towers are not built yet.
    if (_towersPerDim[0] == 0) {
      return {_boxMin, _boxMax};
    }
    // towerBoxMin[0/1] does NOT start at _haloBoxMin[0/1] because the halo towers might extend beyond the halo box.
    const std::array<double, 3> towerBoxMin{
        _boxMin[0] - _towerSideLength[0] * _numTowersPerInteractionLength +
            _towerSideLength[0] * static_cast<double>(index2D[0]),
        _boxMin[1] - _towerSideLength[1] * _numTowersPerInteractionLength +
            _towerSideLength[1] * static_cast<double>(index2D[1]),
        _haloBoxMin[2],
    };
    const std::array<double, 3> towerBoxMax{
        towerBoxMin[0] + _towerSideLength[0],
        towerBoxMin[1] + _towerSideLength[1],
        _haloBoxMax[2],
    };
    return {towerBoxMin, towerBoxMax};
  }

  /**
   * Returns the 2D index of the tower in the tower grid the given 3D coordinates are in.
   * If the location is outside of the domain, the tower nearest tower is returned.
   * * @param pos The 3D coordinates
   * @return array<X, Y>
   */
  [[nodiscard]] std::array<size_t, 2> getTowerIndex2DAtPosition(const std::array<double, 3> &pos) const {
    // special case: Towers are not yet built, then everything is in this one tower.
    if (_towers.size() == 1) {
      return {0, 0};
    }

    std::array<size_t, 2> towerIndex2D{};

    for (int dim = 0; dim < 2; dim++) {
      const auto towerDimIndex =
          static_cast<long int>(std::floor((pos[dim] - _boxMin[dim]) * _towerSideLengthReciprocal[dim])) +
          _numTowersPerInteractionLength;
      const auto towerDimIndexNonNegative = static_cast<size_t>(std::max(towerDimIndex, 0l));
      const auto towerDimIndexNonLargerValue = std::min(towerDimIndexNonNegative, _towersPerDim[dim] - 1);
      towerIndex2D[dim] = towerDimIndexNonLargerValue;
      /// @todo this is a sanity check to prevent doubling of particles, but could be done better! e.g. by border and
      // flag manager
      if (pos[dim] >= _haloBoxMax[dim]) {
        towerIndex2D[dim] = _towersPerDim[dim] - 1;
      } else if (pos[dim] < _haloBoxMin[dim]) {
        towerIndex2D[dim] = 0;
      }
    }

    return towerIndex2D;
  }

  /**
   * Return the 1D index of the tower at a given position
   * @param pos
   * @return
   */
  [[nodiscard]] size_t getTowerIndex1DAtPosition(const std::array<double, 3> &pos) const {
    const auto [x, y] = getTowerIndex2DAtPosition(pos);
    return towerIndex2DTo1D(x, y);
  }

  /**
   * Return a reference to the tower at a given position in the simulation coordinate system (e.g. particle position).
   * @param pos
   * @return
   */
  ClusterTower<Particle_T> &getTowerAtPosition(const std::array<double, 3> &pos) {
    const auto [x, y] = getTowerIndex2DAtPosition(pos);
    return getTowerByIndex2D(x, y);
  }

  /**
   * Returns a reference to the tower for the given tower grid coordinates.
   * @param x The x-th tower in x direction.
   * @param y The y-th tower in y direction.
   * @return a reference to the tower for the given tower grid coordinates.
   */
  ClusterTower<Particle_T> &getTowerByIndex2D(const size_t x, const size_t y) {
    return _towers[towerIndex2DTo1D(x, y)];
  }

  /**
   * Returns the 1D index for the given 2D-coordinates of a tower.
   *
   * @param x The x-coordinate of the tower.
   * @param y The y-coordinate of the tower.
   * @return the 1D index for the given 2D-coordinates of a tower.
   */
  [[nodiscard]] size_t towerIndex2DTo1D(const size_t x, const size_t y) const { return x + y * _towersPerDim[0]; }

  /**
   * Returns the 2D index for the given 1D index of a tower.
   *
   * @param index
   * @return the 2D index for the given 1D index of a tower.
   */
  [[nodiscard]] std::array<size_t, 2> towerIndex1DTo2D(size_t index) const {
    if (_towersPerDim[0] == 0) {
      return {0, 0};
    } else {
      return {index % _towersPerDim[0], index / _towersPerDim[0]};
    }
  }

  /**
   * Getter
   * @return
   */
  size_t getFirstOwnedTowerIndex() const { return _firstOwnedTowerIndex; }
  /**
   * Getter
   * @return
   */
  size_t getLastOwnedTowerIndex() const { return _lastOwnedTowerIndex; }
  /**
   * Getter
   * @return
   */
  double getInteractionLength() const { return _interactionLength; }
  /**
   * Getter
   * @return
   */
  int getNumTowersPerInteractionLength() const { return _numTowersPerInteractionLength; }
  /**
   * Getter
   * @return
   */
  const std::array<double, 3> &getBoxMin() const { return _boxMin; }
  /**
   * Getter
   * @return
   */
  const std::array<double, 3> &getBoxMax() const { return _boxMax; }
  /**
   * Getter
   * @return
   */
  const std::array<double, 3> &getHaloBoxMin() const { return _haloBoxMin; }
  /**
   * Getter
   * @return
   */
  const std::array<double, 3> &getHaloBoxMax() const { return _haloBoxMax; }
  /**
   * Getter
   * @return
   */
  const std::vector<ClusterTower<Particle_T>> &getTowers() const { return _towers; }
  /**
   * Getter for a mutable reference
   * @return
   */
  std::vector<ClusterTower<Particle_T>> &getTowersRef() { return _towers; }
  /**
   * Getter
   * @return
   */
  const std::array<size_t, 2> &getTowersPerDim() const { return _towersPerDim; }
  /**
   * Getter
   * @return
   */
  const std::array<double, 2> &getTowerSideLength() const { return _towerSideLength; }
  /**
   * Getter
   * @return
   */
  const std::array<double, 2> &getTowerSideLengthReciprocal() const { return _towerSideLengthReciprocal; }

  bool cellCanContainHaloParticles(index_t index1d) const override {
    // Always true because towers cover the whole z dimension
    return true;
  }

  bool cellCanContainOwnedParticles(index_t index1d) const override {
    // check if the tower is strictly in the halo region
    if (index1d < _firstOwnedTowerIndex or index1d > _lastOwnedTowerIndex) {
      return false;
    }
    const auto index2D = towerIndex1DTo2D(index1d);
    bool isHaloTower = false;
    for (size_t i = 0; i < index2D.size(); ++i) {
      if (index2D[i] < _numTowersPerInteractionLength or
          index2D[i] >= _towersPerDim[i] - _numTowersPerInteractionLength) {
        isHaloTower = true;
        break;
      }
    }
    return not isHaloTower;
  }

 private:
  /**
   * Index of the first tower in _towers that can contain owned particles.
   */
  size_t _firstOwnedTowerIndex{};
  /**
   * Index of the last tower in _towers that can contain owned particles.
   */
  size_t _lastOwnedTowerIndex{};
  /**
   * cutoff + skin
   */
  double _interactionLength;
  /**
   * The interaction length in number of towers it reaches.
   * static_cast<int>(std::ceil((this->getInteractionLength()) * _towerSideLengthReciprocal))
   */
  int _numTowersPerInteractionLength{};

  /**
   * Minimum of the container.
   */
  std::array<double, 3> _boxMin;

  /**
   * Maximum of the container.
   */
  std::array<double, 3> _boxMax;

  /**
   * Minimum of the container including halo.
   */
  std::array<double, 3> _haloBoxMin;

  /**
   * Maximum of the container including halo.
   */
  std::array<double, 3> _haloBoxMax;

  /**
   * Internal storage, particles are split into a grid in xy-dimension.
   */
  std::vector<ClusterTower<Particle_T>> _towers;

  /**
   * Dimensions of the 2D xy-grid including halo.
   */
  std::array<size_t, 2> _towersPerDim{};

  /**
   * Side length of xy-grid cells.
   */
  std::array<double, 2> _towerSideLength{0.};

  /**
   * 1/side length of xy-grid cells.
   */
  std::array<double, 2> _towerSideLengthReciprocal{0.};
};
}  // namespace autopas::internal
