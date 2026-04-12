/**
 * @file IterationMeasurements.h
 * @authors muehlhaeusser
 * @date 05.06.2024
 */

#pragma once

namespace autopas {
/**
 * Struct to collect all sorts of measurements taken during a computeInteractions iteration.
 */
struct IterationMeasurements {
  /**
   * Time it takes for the LogicHandler's computeInteractions() function.
   */
  long timeComputeInteractions{};

  /**
   * Time it takes for the Remainder Traversal.
   */
  long timeRemainderTraversal{};

  /**
   * Time it takes for rebuilding neighbor lists.
   */
  long timeRebuild{};

  /**
   * Time it takes for the complete iteratePairwise pipeline.
   */
  long timeTotal{};

  /**
   * Bool whether energy measurements are currently possible.
   */
  bool energyMeasurementsPossible{false};

  /**
   * Average energy consumed per time in Watts.
   */
  double energyWatts{};

  /**
   * Total energy consumed in Joules.
   */
  double energyJoules{};

  /**
   * Time in seconds during which energy was consumed.
   */
  double energyDeltaT{};

  /**
   * Total energy consumed during rebuilding.
   */
  long energyTotalRebuild{};

  /**
   * Total energy consumed during compute interactions and remainder traversal.
   */
  long energyTotalNonRebuild{};

  /**
   * Total energy consumed so far
   */
  long energyTotal{};
};
}  // namespace autopas
