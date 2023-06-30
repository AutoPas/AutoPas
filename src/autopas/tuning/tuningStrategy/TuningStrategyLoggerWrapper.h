/**
 * @file TuningStrategyLoggerWrapper.h
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
 * This wrapper can wrap any tuning strategy and produce a log of its activity.
 * All methods first invoke a log call and are then passed on to the wrapped strategy.
 *
 * The log can be replayed using the TuningStrategyLogReplayer.
 */
class TuningStrategyLoggerWrapper final : public TuningStrategyInterface {
 public:
  /**
   * Creates a wrapper logger for a tuning strategy.
   * @param actualTuningStrategy The tuning strategy to hand all calls over.
   * @param outputSuffix The output suffix for all written log files.
   */
  explicit TuningStrategyLoggerWrapper(std::unique_ptr<TuningStrategyInterface> actualTuningStrategy,
                                       const std::string &outputSuffix);

  ~TuningStrategyLoggerWrapper() override;

  /**
   * Logs evidence and hands it to wrapped.
   * @param configuration
   * @param evidence
   */
  void addEvidence(const Configuration &configuration, const Evidence &evidence) override;

  /**
   * Logs call to tuneConfiguration, then hands it over to wrapped.
   * @param currentInvalid
   * @return
   */
  void optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;

  /**
   * Logs call to reset, then hands it over to wrapped.
   * @param iteration
   */
  void reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const autopas::EvidenceCollection &evidenceCollection) override;

  /**
   * Returns true. The live info is always needed for logging.
   * @return true
   */
  [[nodiscard]] bool needsLiveInfo() const override;

  /**
   * Logs the live info and hands it over to wrapped.
   * @param info
   */
  void receiveLiveInfo(const LiveInfo &info) override;

  /**
   * Hands over call to wrapped
   * @return
   */
  bool smoothedHomogeneityAndMaxDensityNeeded() const override;

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

}  // namespace autopas
