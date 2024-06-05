/**
 * @file FuzzyTuning.h
 * @author Manuel Lerchner
 * @date 17.04.24
 */

#pragma once

#include <string>

#include "autopas/tuning/tuningStrategy/TuningStrategyInterface.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/FuzzyControlSystem.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/LinguisticVariable.h"
#include "tuning/tuningStrategy/fuzzyTuning/fuzzyController/LinguisticVariable.h"

namespace autopas {

using namespace autopas::fuzzy_logic;

class FuzzyTuning : public TuningStrategyInterface {
 public:
  explicit FuzzyTuning(std::string fuzzyRuleFileName);

  TuningStrategyOption getOptionType() const override;

  bool needsLiveInfo() const override;

  void receiveLiveInfo(const LiveInfo &info) override;

  void addEvidence(const Configuration &configuration, const Evidence &evidence) override;

  void reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const autopas::EvidenceCollection &evidenceCollection) override;

  void optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;

 private:
  /**
   * Parses the given fuzzy rule file and returns a map of the form {consequent_domain: FuzzyControlSystem}.
   * @param fuzzyRuleFilename The name of the fuzzy rule file.
   * @return A map of the form {consequent_domain: FuzzyControlSystem}.
   */
  [[nodiscard]] static std::pair<std::vector<std::shared_ptr<LinguisticVariable>>,
                                 std::map<std::string, std::shared_ptr<FuzzyControlSystem>>>
  parse(const std::string &fuzzyRuleFilename);

  /**
   * The name of the fuzzy rule file.
   */
  std::string _fuzzyRuleFileName;

  /**
   * The current live info used to make predictions.
   */
  std::map<std::string, double> _currentLiveInfo;

  /**
   * The fuzzy control systems parsed from the fuzzy rule file.
   */
  std::map<std::string, std::shared_ptr<FuzzyControlSystem>> _fuzzyControlSystems;
};

};  // namespace autopas
