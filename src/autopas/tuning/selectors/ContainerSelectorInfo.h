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
        verletSkinPerTimestep(0.),
        verletRebuildFrequency(0),
        verletClusterSize(64),
        loadEstimator(autopas::LoadEstimatorOption::none),
        numberOfHGLevels(1) {}

  /**
   * Constructor.
   * @param cellSizeFactor Cell size factor to be used in this container (only relevant for LinkedCells, VerletLists and
   * VerletListsCells).
   * @param verletSkinPerTimestep Length added to the cutoff for the verlet lists' skin per timestep inbetween
   * rebuilding lists.
   * @param verletRebuildFrequency rebuild frequency.
   * @param verletClusterSize Size of verlet Clusters
   * @param loadEstimator load estimation algorithm for balanced traversals.
   */
  explicit ContainerSelectorInfo(double cellSizeFactor, double verletSkinPerTimestep,
                                 unsigned int verletRebuildFrequency, unsigned int verletClusterSize,
                                 autopas::LoadEstimatorOption loadEstimator,
                                 unsigned int numberOfHGLevels)
      : cellSizeFactor(cellSizeFactor),
        verletSkinPerTimestep(verletSkinPerTimestep),
        verletRebuildFrequency(verletRebuildFrequency),
        verletClusterSize(verletClusterSize),
        loadEstimator(loadEstimator),
        numberOfHGLevels(numberOfHGLevels) {}

  /**
   * Equality between ContainerSelectorInfo
   * @param other
   * @return True iff all member equal
   */
  bool operator==(const ContainerSelectorInfo &other) const {
    return cellSizeFactor == other.cellSizeFactor and verletSkinPerTimestep == other.verletSkinPerTimestep and
           verletClusterSize == other.verletClusterSize and loadEstimator == other.loadEstimator and
           numberOfHGLevels == other.numberOfHGLevels;
  }

  /**
   * Inequality between ContainerSelectorInfo
   * @param other
   * @return False iff all member euqal
   */
  bool operator!=(const ContainerSelectorInfo &other) const { return !(*this == other); }

  /**
   * Comparison operator for ContainerSelectorInfo objects.
   * Configurations are compared member wise in the order: _cellSizeFactor, _verletSkinPerTimestep,
   * _verlerRebuildFrequency, loadEstimator
   *
   * @param other
   * @return
   */
  bool operator<(const ContainerSelectorInfo &other) {
    return std::tie(cellSizeFactor, verletSkinPerTimestep, verletRebuildFrequency, verletClusterSize, loadEstimator,
                    numberOfHGLevels) <
           std::tie(other.cellSizeFactor, other.verletSkinPerTimestep, other.verletRebuildFrequency,
                    other.verletClusterSize, other.loadEstimator, other.numberOfHGLevels);
  }

  /**
   * cellSizeFactor Cell size factor to be used in this container (only relevant for LinkedCells)
   */
  double cellSizeFactor;
  /**
   * Length added to the cutoff for the verlet lists' skin per timestep inbetween rebuilding lists.
   */
  double verletSkinPerTimestep;
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
   * Number of the Hierarchical Grid (H-Grid) levels (DEM only).
   */
  unsigned int numberOfHGLevels;
};

}  // namespace autopas
