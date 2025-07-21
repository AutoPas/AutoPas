/**
 * @file ContainerSelectorInfo.h
 * @author Jan Nguyen
 * @date 05.06.19
 */

#pragma once
#include <array>
#include <memory>

#include "autopas/options/LoadEstimatorOption.h"

namespace autopas {
/**
 * Info to generate a container.
 */
class ContainerSelectorInfo {
 public:
  /**
   * Default Constructor.
   */
  ContainerSelectorInfo()
      : cellSizeFactor(1.),
        verletSkin(0.),
        verletRebuildFrequency(0),
        verletClusterSize(64),
        loadEstimator(autopas::LoadEstimatorOption::none),
        cao(1),
        alpha(0.5),
        grid_dims({0,0,0}) {}

  /**
   * Constructor.
   * @param cellSizeFactor Cell size factor to be used in this container (only relevant for LinkedCells, VerletLists and
   * VerletListsCells).
   * @param verletSkin Length added to the cutoff for the verlet lists' skin per timestep inbetween
   * rebuilding lists.
   * @param verletRebuildFrequency rebuild frequency.
   * @param verletClusterSize Size of verlet Clusters
   * @param loadEstimator load estimation algorithm for balanced traversals.
   */
  explicit ContainerSelectorInfo(double cellSizeFactor, double verletSkin, unsigned int verletRebuildFrequency,
                                 unsigned int verletClusterSize, autopas::LoadEstimatorOption loadEstimator, unsigned int cao, double alpha, std::array<unsigned int, 3> grid_dims)
      : cellSizeFactor(cellSizeFactor),
        verletSkin(verletSkin),
        verletRebuildFrequency(verletRebuildFrequency),
        verletClusterSize(verletClusterSize),
        loadEstimator(loadEstimator),
        cao(cao),
        alpha(alpha),
        grid_dims(grid_dims) {}

  /**
   * Equality between ContainerSelectorInfo
   * @param other
   * @return True iff all member equal
   */
  bool operator==(const ContainerSelectorInfo &other) const {
    return cellSizeFactor == other.cellSizeFactor and verletSkin == other.verletSkin and
           verletClusterSize == other.verletClusterSize and loadEstimator == other.loadEstimator and cao == other.cao and alpha == other.alpha and grid_dims == other.grid_dims;
  }

  /**
   * Inequality between ContainerSelectorInfo
   * @param other
   * @return False iff all member euqal
   */
  bool operator!=(const ContainerSelectorInfo &other) const { return !(*this == other); }

  /**
   * Comparison operator for ContainerSelectorInfo objects.
   * Configurations are compared member wise in the order: _cellSizeFactor, _verletSkin,
   * _verlerRebuildFrequency, loadEstimator, cao, grid_dims
   *
   * @param other
   * @return
   */
  bool operator<(const ContainerSelectorInfo &other) {
    return std::tie(cellSizeFactor, verletSkin, verletRebuildFrequency, verletClusterSize, loadEstimator, cao, alpha, grid_dims) <
           std::tie(other.cellSizeFactor, other.verletSkin, other.verletRebuildFrequency, other.verletClusterSize,
                    other.loadEstimator, other.cao, other.alpha, other.grid_dims);
  }

  /**
   * cellSizeFactor Cell size factor to be used in this container (only relevant for LinkedCells)
   */
  double cellSizeFactor;
  /**
   * Length added to the cutoff for the verlet lists' skin inbetween rebuilding lists.
   */
  double verletSkin;
  /**
   * The rebuild frequency.
   */
  unsigned int verletRebuildFrequency;
  /**
   * Size of Verlet Clusters
   */
  unsigned int verletClusterSize;
  /**
   * Load estimator for balanced sliced traversals.
   */
  autopas::LoadEstimatorOption loadEstimator;
  /**
   * Charge assignment order for P3M
   */
  unsigned int cao;
  /**
   * Ewald parameter for P3M
   */
  double alpha;
  /**
   * Number of gridpoints per dimension for P3M
   */
  std::array<unsigned int, 3> grid_dims;
};

}  // namespace autopas
