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
   * Number of rebuild per iteration
   */
  long numberRebuild{};

  /**
   * Number of fast particles per iteration
   */
  long numberFastParticles{};

  /**
   * Size of particle buffer at the end of updateContainer()
   */
  size_t particleBufferSize{};

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
