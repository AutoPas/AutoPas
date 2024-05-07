/**
 * @file FuzzyTuning.h
 * @author Manuel Lerchner
 * @date 17.04.24
 */

#pragma once

#include <string>

#include "autopas/tuning/tuningStrategy/TuningStrategyInterface.h"

namespace autopas {

class FuzzyTuning : public TuningStrategyInterface {
 public:
  FuzzyTuning(const std::string &fuzzyRuleFileName);

  TuningStrategyOption getOptionType() override;

  bool needsLiveInfo() const override;

  void receiveLiveInfo(const LiveInfo &info) override;

  void addEvidence(const Configuration &configuration, const Evidence &evidence) override;

  void reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const autopas::EvidenceCollection &evidenceCollection) override;

  void optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;

 private:
  std::string _fuzzyRuleFileName;

  LiveInfo _currentLiveInfo;

  int _tuningPhase;
};
};  // namespace autopas
