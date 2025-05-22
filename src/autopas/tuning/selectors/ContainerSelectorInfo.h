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
      : boxMin({0., 0., 0.}),
        boxMax({1., 1., 1.}),
        cutoff(1.),
        cellSizeFactor(1.),
        verletSkin(0.),
        verletRebuildFrequency(0),
        verletClusterSize(64),
        loadEstimator(LoadEstimatorOption::none) {}

  /**
   * Constructor.
   * @param boxMin Lower corner of the container.
   * @param boxMax Upper corner of the container.
   * @param cutoff Cutoff radius to be used in this container.
   * @param cellSizeFactor Cell size factor to be used in this container (only relevant for LinkedCells, VerletLists and
   * VerletListsCells).
   * @param verletSkin Length added to the cutoff for the verlet lists' skin per timestep inbetween
   * rebuilding lists.
   * @param verletRebuildFrequency rebuild frequency.
   * @param verletClusterSize Size of verlet Clusters
   * @param loadEstimator load estimation algorithm for balanced traversals.
   */
  explicit ContainerSelectorInfo(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, double cutoff, double cellSizeFactor, double verletSkin, unsigned int verletRebuildFrequency,
                                 unsigned int verletClusterSize, LoadEstimatorOption loadEstimator)
      : boxMin(boxMin),
        boxMax(boxMax),
        cutoff(cutoff),
        cellSizeFactor(cellSizeFactor),
        verletSkin(verletSkin),
        verletRebuildFrequency(verletRebuildFrequency),
        verletClusterSize(verletClusterSize),
        loadEstimator(loadEstimator) {}

  /**
   * Equality between ContainerSelectorInfo
   * @param other
   * @return True iff all members equal
   */
  bool operator==(const ContainerSelectorInfo &other) const {
    return cellSizeFactor == other.cellSizeFactor and verletSkin == other.verletSkin and
           verletClusterSize == other.verletClusterSize and loadEstimator == other.loadEstimator;
  }

  /**
   * Inequality between ContainerSelectorInfo
   * @param other
   * @return False iff all members equal
   */
  bool operator!=(const ContainerSelectorInfo &other) const { return !(*this == other); }

  /**
   * Comparison operator for ContainerSelectorInfo objects.
   * Configurations are compared member wise in the order: _cellSizeFactor, _verletSkin,
   * _verletRebuildFrequency, loadEstimator
   *
   * @param other
   * @return
   */
  bool operator<(const ContainerSelectorInfo &other) {
    return std::tie(cellSizeFactor, verletSkin, verletRebuildFrequency, verletClusterSize, loadEstimator) <
           std::tie(other.cellSizeFactor, other.verletSkin, other.verletRebuildFrequency, other.verletClusterSize,
                    other.loadEstimator);
  }

  /**
   * Lower corner of the container.
   */
  std::array<double, 3> boxMin;

  /**
   * Upper corner of the container.
   */
  std::array<double, 3> boxMax;

  /**
   * Cutoff radius to be used in this container.
   */
  const double cutoff;

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
  LoadEstimatorOption loadEstimator;
};

}  // namespace autopas
