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
      : cellSizeFactor(1.), verletSkinPerTimestep(0.),verletRebuildFrequency(0.), verletClusterSize(64), loadEstimator(autopas::LoadEstimatorOption::none) {}

  /**
   * Constructor.
   * @param cellSizeFactor Cell size factor to be used in this container (only relevant for LinkedCells, VerletLists and
   * VerletListsCells).
   * @param verletSkinPerTimestep Length added to the cutoff for the verlet lists' skin per timestep.
   * @param verletClusterSize Size of verlet Clusters
   * @param loadEstimator load estimation algorithm for balanced traversals.
   */
  explicit ContainerSelectorInfo(double cellSizeFactor, double verletSkinPerTimestep, unsigned int verletClusterSize,
                                 autopas::LoadEstimatorOption loadEstimator)
      : cellSizeFactor(cellSizeFactor),
        verletSkinPerTimestep(verletSkinPerTimestep),
        verletRebuildFrequency(verletRebuildFrequency),
        verletClusterSize(verletClusterSize),
        loadEstimator(loadEstimator) {}

  /**
   * Equality between ContainerSelectorInfo
   * @param other
   * @return True iff all member equal
   */
  bool operator==(const ContainerSelectorInfo &other) const {
    return cellSizeFactor == other.cellSizeFactor and verletSkinPerTimestep == other.verletSkinperTimestep and
           verletClusterSize == other.verletClusterSize and loadEstimator == other.loadEstimator;
  }

  /**
   * Inequality between ContainerSelectorInfo
   * @param other
   * @return False iff all member euqal
   */
  bool operator!=(const ContainerSelectorInfo &other) const { return !(*this == other); }

  /**
   * Comparison operator for ContainerSelectorInfo objects.
   * Configurations are compared member wise in the order: _cellSizeFactor, _verletSkin, _verlerRebuildFrequency,
   * loadEstimator
   *
   * @param other
   * @return
   */
  bool operator<(const ContainerSelectorInfo &other) {
    return std::tie(cellSizeFactor, verletSkinPerTimestep, verletClusterSize, loadEstimator) <
           std::tie(other.cellSizeFactor, other.verletSkinPerTimestep, other.verletClusterSize, other.loadEstimator);
  }

  /**
   * cellSizeFactor Cell size factor to be used in this container (only relevant for LinkedCells)
   */
  double cellSizeFactor;
  /**
   * Length added to the cutoff for the verlet lists' skin.
   */
  double verletSkinPerTimestep;
  /**
   * Length added to the cutoff for the verlet lists' skin.
   */
  double verletRebuildFrequency;
  /**
   * Size of Verlet Clusters
   */
  unsigned int verletClusterSize;
  /**
   * Load estimator for balanced sliced traversals.
   */
  autopas::LoadEstimatorOption loadEstimator;
};

}  // namespace autopas
