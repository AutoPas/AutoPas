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
   * Time it takes for the LogicHandler's iteratePairwise() function.
   */
  long timeIteratePairwise{};

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
   * Total energy consumed in watts.
   */
  double energyWatts{};

  /**
   * Total energy consumed in joules.
   */
  double energyJoules{};

  /**
   * Time in which energy was consumed.
   */
  double energySeconds{};

  /**
   * Total energy consumed so far
   */
  long energyTotal{};

};
}  // namespace autopas
