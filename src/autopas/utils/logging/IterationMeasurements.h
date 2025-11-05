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
   * Total energy consumed so far
   */
  long energyTotal{};

  /**
   * Size of particle buffer at the end of updateContainer()
   */
  size_t particleBufferSize{};

  /**
   * Number of owned particles per iteration
   */
  size_t numParticlesOwned{};

  /**
   * Number of halo particles per iteration
   */
  size_t numParticlesHalo{};

  /**
   * Number of fast particles per iteration
   */
  size_t numParticlesFast{};
};
}  // namespace autopas
