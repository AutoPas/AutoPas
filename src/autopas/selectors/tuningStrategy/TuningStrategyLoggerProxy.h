/**
 * @file TuningStrategyLoggerProxy.h
 * @author humig
 * @date 24.09.2021
 */

#pragma once

#include <any>
#include <fstream>
#include <optional>

#include "TuningStrategyInterface.h"

namespace autopas {

/**
 * A proxy for a tuning strategy that logs all calls and their parameters. Just hands over all calls to a wrapped real
 * tuning strategy.
 *
 * The log can be replayed using the TuningStrategyLogReplayer.
 */
class TuningStrategyLoggerProxy : public TuningStrategyInterface {
 public:
  /**
   * Creates a tuning strategy logger proxy.
   * @param actualTuningStrategy The tuning strategy to hand all calls over.
   * @param outputSuffix The output suffix for all written log files.
   */
  explicit TuningStrategyLoggerProxy(std::unique_ptr<TuningStrategyInterface> actualTuningStrategy,
                                     const std::string &outputSuffix);

  ~TuningStrategyLoggerProxy() override;

  /**
   * Logs evidence and hands it to wrapped.
   * @param time
   * @param iteration
   */
  void addEvidence(long time, size_t iteration) override;
  /**
   * Only hands over call to wrapped tuning strategy.
   * @param configuration
   * @return
   */
  [[nodiscard]] long getEvidence(Configuration configuration) const override;
  /**
   * Only hands over call to wrapped tuning strategy.
   * @return
   */
  [[nodiscard]] const Configuration &getCurrentConfiguration() const override;

  /**
   * Logs call to tune, then hands it over to wrapped.
   * @param currentInvalid
   * @return
   */
  bool tune(bool currentInvalid = false) override;

  /**
   * Logs call to reset, then hands it over to wrapped.
   * @param iteration
   */
  void reset(size_t iteration) override;

  /**
   * Returns true. The live info is always needed for logging.
   * @return true
   */
  [[nodiscard]] bool needsLiveInfo() const override;

  /**
   * Logs the live info and hands it over to wrapped.
   * @param info
   */
  void receiveLiveInfo(const LiveInfo& info) override;

  /**
   * Hands call over to wrapped.
   * @return
   */
  [[nodiscard]] std::set<ContainerOption> getAllowedContainerOptions() const override;

  /**
   * Hands over call to wrapped.
   */
  void removeN3Option(Newton3Option) override;

  /**
   * Hands over call to wrapped.
   * @return
   */
  [[nodiscard]] bool searchSpaceIsTrivial() const override;

  /**
   * Hands over call to wrapped.
   * @return
   */
  [[nodiscard]] bool searchSpaceIsEmpty() const override;

 private:
  /**
   * The wrapped real tuning strategy that is used to in almost all method implementations.
   */
  std::unique_ptr<TuningStrategyInterface> _actualTuningStrategy;
  /**
   * The ofstream of the log file.
   */
  std::ofstream _logOut;
};

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
