/**
 * @file TuningStrategyLogger.h
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
 * This dummy strategy logs any call made to the list of strategies.
 * All methods only invoke a logger and do nothing else.
 *
 * The log can be replayed using the TuningStrategyLogReplayer.
 */
class TuningStrategyLogger final : public TuningStrategyInterface {
 public:
  /**
   * Creates a wrapper logger for a tuning strategy.
   * @param outputSuffix The output suffix for all written log files.
   */
  explicit TuningStrategyLogger(const std::string &outputSuffix);

  ~TuningStrategyLogger() override;

  void addEvidence(const Configuration &configuration, const Evidence &evidence) override;

  void optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;

  void reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const autopas::EvidenceCollection &evidenceCollection) override;

  [[nodiscard]] bool needsLiveInfo() const override;

  void receiveLiveInfo(const LiveInfo &info) override;

 private:
  /**
   * The ofstream of the log file.
   */
  std::ofstream _logOut;
};

}  // namespace autopas
