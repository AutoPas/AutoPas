/**
 * @file TuningDataLogger.h
 * @author F. Gratl
 * @date 25/01/2021
 */

#pragma once

#include "autopas/tuning/Configuration.h"

namespace autopas {

/**
 * Helper to log results of the tuning process to a csv file for easier analysis.
 *
 * It uses an asynchronous spd logger to write a csv file named "AutoPas_iterationPerformance_<dateStamp>.csv".
 *
 * By default logging the data is disabled. It can be enabled by setting the cmake variable AUTOPAS_LOG_TUNINGDATA
 * to ON.
 */
class TuningDataLogger {
 public:
  /**
   * Constructor initializes the logger and sets the output file name.
   * @param numSamples Number of samples that are taken per configuration.
   * @param rebuildFrequency Number of iterations at which the neighbor lists are updated. This is used here to evaluate
   * how many rebuild samples we will collect.
   * @param outputSuffix Suffix for all output files produced by this class.
   */
  explicit TuningDataLogger(size_t numSamples, size_t rebuildFrequency, const std::string &outputSuffix = "");

  /**
   * Destructor drops the logger from the spd registry.
   */
  ~TuningDataLogger();

  /**
   * Log the result of a tuning phase.
   * @param configuration
   * @param samplesRebuildingNeighborLists
   * @param samplesNotRebuildingNeighborLists
   * @param iteration
   * @param reducedValue
   * @param smoothedVale
   * @param meanRebuildFrequency
   */
  void logTuningData(const autopas::Configuration &configuration,
                     const std::vector<long> &samplesRebuildingNeighborLists,
                     const std::vector<long> &samplesTraverseInteractions, size_t iteration, long reducedValue,
                     long smoothedVale, double meanRebuildFrequency) const;

 private:
  std::string _loggerName;
};

}  // namespace autopas