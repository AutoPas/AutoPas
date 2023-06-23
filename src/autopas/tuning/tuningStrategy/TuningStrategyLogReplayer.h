/**
 * @file TuningStrategyLogReplayer.h
 * @author humig
 * @date 24.09.2021
 */

#pragma once

#include <memory>
#include <optional>
#include <string>

#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyInterface.h"

namespace autopas {
/**
 * This class is able to replay a log file to a tuning strategy to observe its behavior in the logged scenario. Can only
 * provide evidence for configurations that has actually been logged, so make sure the necessary data is in the log.
 * Collecting data with FullSearch ensures this.
 */
class TuningStrategyLogReplayer {
 public:
  /**
   * Creates a log replayer from a log file and a tuning strategy to replay to.
   * @param filename The name of the log file.
   * @param tuningStrategy The tuning strategy to replay to.
   */
  TuningStrategyLogReplayer(std::string filename, std::shared_ptr<TuningStrategyInterface> tuningStrategy);

  /**
   * Replays the log to the tuning strategy.
   * @return The best configuration found by the tuning strategy, if there is any
   */
  std::optional<Configuration> replay();

 private:
  /**
   * The filename of the log file to replay.
   */
  std::string _filename;
  /**
   * The tuning strategy to replay the log to.
   */
  std::shared_ptr<TuningStrategyInterface> _tuningStrategy;
};
}  // namespace autopas