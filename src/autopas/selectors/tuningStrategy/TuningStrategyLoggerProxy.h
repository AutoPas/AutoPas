/**
 * @file TuningStrategyLoggerProxy.h
 * @author humig
 * @date 24.09.2021
 */

#pragma once

#include <any>
#include <fstream>

#include "TuningStrategyInterface.h"

namespace autopas {

class TuningStrategyLoggerProxy : public TuningStrategyInterface {
 public:
  explicit TuningStrategyLoggerProxy(std::unique_ptr<TuningStrategyInterface> actualTuningStrategy,
                                     const std::string &outputSuffix);

  ~TuningStrategyLoggerProxy() override;

  void addEvidence(long time, size_t iteration) override;
  [[nodiscard]] long getEvidence(Configuration configuration) const override;
  [[nodiscard]] const Configuration &getCurrentConfiguration() const override;

  bool tune(bool currentInvalid = false) override;

  void reset(size_t iteration) override;

  [[nodiscard]] bool needsLiveInfo() const override;

  void receiveLiveInfo(LiveInfo info) override;

  [[nodiscard]] std::set<ContainerOption> getAllowedContainerOptions() const override;

  void removeN3Option(Newton3Option) override;

  [[nodiscard]] bool searchSpaceIsTrivial() const override;

  [[nodiscard]] bool searchSpaceIsEmpty() const override;

 private:
  std::unique_ptr<TuningStrategyInterface> _actualTuningStrategy;
  std::ofstream _logOut;
};

class TuningStrategyLogReplayer {
 public:
  TuningStrategyLogReplayer(std::string filename, std::unique_ptr<TuningStrategyInterface> tuningStrategy);

  void replay() const;

 private:
  std::string _filename;
  std::unique_ptr<TuningStrategyInterface> _tuningStrategy;
};
}  // namespace autopas
