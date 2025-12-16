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
   * Number of particles in the buffer during remainder Traversal
   */
  size_t numParticlesBuffer{};

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

  /**
   * Estimate of particles in the buffer used for predicting estimateRemainderTraversalTime
   */
  size_t estimateNumParticlesBuffer{};

  /**
   * Linear regression estimate of remainderTraversalTime given estimated amount of particles in the buffer
   */
  double estimateRemainderTraversalTime{};

  /**
   * Mean of rebuildNeighborTime in the current iteration excluding tuning rebuilds
   */
  double estimateRebuildNeighborTime{};

  /**
   * True if dynamic rebuild would be initiated in the current iteration otherwise false
   */
  bool doDynamicRebuild{false};
};
}  // namespace autopas
