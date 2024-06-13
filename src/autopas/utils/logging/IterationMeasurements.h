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
   * Energy consumed by the entire system.
   */
  double energyPsys{};

  /**
   * Energy consumed by the CPU package.
   */
  double energyPkg{};

  /**
   * Energy consumed by the RAM.
   */
  double energyRam{};

  /**
   * Total energy consumed. This is identical to energyPsys if it is available, otherwise it is the sum of all other
   * energy sources.
   */
  long energyTotal{};
};
}  // namespace autopas
