/**
 * @file LogicHandlerInfo.h
 * @author F. Gratl
 * @date 19/06/2023
 */

#pragma once

#include "array"
#include "string"

namespace autopas {
/**
 * Class that wraps all arguments for the logic handler to provide a more stable API.
 */
class LogicHandlerInfo {
 public:
  /**
   * Lower corner of the container.
   */
  std::array<double, 3> boxMin;
  /**
   * Upper corner of the container.
   */
  std::array<double, 3> boxMax;
  /**
   * Cutoff radius to be used in this simulation.
   */
  double cutoff;
  /**
   * Length added to the cutoff for the Verlet lists' skin.
   */
  double verletSkinPerTimestep;
  /**
   * The rebuild frequency for all containers.
   */
  unsigned int rebuildFrequency;
  /**
   * Number of particles in a cluster to use in VCL.
   */
  unsigned int verletClusterSize;
  /**
   * Suffix for all output (log) files produced by the LogicHandler.
   */
  std::string outputSuffix;
};
}  // namespace autopas
