/**
 * @file RuleBasedTuning.h
 * @author humig
 * @date 28.06.2021
 */

#pragma once

#include <fstream>
#include <iostream>
#include <list>
#include <variant>
#include <unordered_set>

#include "RuleBasedProgramParser.h"
#include "RuleBasedProgramTree.h"
#include "RuleVM.h"
#include "autopas/selectors/tuningStrategy/FullSearch.h"
#include "autopas/utils/logging/Logger.h"

namespace autopas {

class RuleBasedTuning : public FullSearch {
 public:
  using PrintTuningErrorFunType = std::function<void(const rule_syntax::ConfigurationOrder& order,
                                                     const Configuration& actualBetterConfig, unsigned long betterRuntime,
                                                     const Configuration& shouldBeBetterConfig, unsigned long shouldBeBetterRuntime,
                                                     const LiveInfo& liveInfo)>;

  RuleBasedTuning(const std::set<ContainerOption> &allowedContainerOptions,
                  const std::set<double> &allowedCellSizeFactors,
                  const std::set<TraversalOption> &allowedTraversalOptions,
                  const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
                  const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                  const std::set<Newton3Option> &allowedNewton3Options, bool verifyModeEnabled = false,
                  std::string ruleFileName = "tuningRules.rule",
                  PrintTuningErrorFunType tuningErrorPrinter = {})
      : FullSearch(allowedContainerOptions, allowedCellSizeFactors, allowedTraversalOptions,
                   allowedLoadEstimatorOptions, allowedDataLayoutOptions, allowedNewton3Options),
        _originalSearchSpace(this->_searchSpace.begin(), this->_searchSpace.end()),
        _verifyModeEnabled(verifyModeEnabled),
        _ruleFileName(std::move(ruleFileName)), _tuningErrorPrinter(std::move(tuningErrorPrinter)) {}

  bool needsLiveInfo() const override { return true; }

  void receiveLiveInfo(const LiveInfo& info) override {
    _currentLiveInfo = info;
    _lastApplicableConfigurationOrders = applyRules();
  }

  void addEvidence(long time, size_t iteration) override {
    FullSearch::addEvidence(time, iteration);
    tuningTime += time;
    tuningTimeLifetime += time;
    if (_verifyModeEnabled) {
      verifyCurrentConfigTime();
      if(_removedConfigurations.find(this->getCurrentConfiguration()) != _removedConfigurations.end()) {
        wouldHaveSkippedTuningTime += time;
        wouldHaveSkippedTuningTimeLifetime += time;
      }
    }
  }

  void reset(size_t iteration) override {
    FullSearch::reset(iteration);

    if(_verifyModeEnabled) {
      if(tuningTime > 0) {
        auto percent = static_cast<double>(wouldHaveSkippedTuningTime) / static_cast<double>(tuningTime) * 100;
        auto percentRounded = std::round(percent * 100) / 100;
        AutoPasLog(info, "Rules would have saved {} ns removing {}/{} configurations. ({}% of total tuning time)",
                   wouldHaveSkippedTuningTime, _removedConfigurations.size(), _originalSearchSpace.size() ,percentRounded);
      }
    }

    tuningTime = 0;
    wouldHaveSkippedTuningTime = 0;
  }

  [[nodiscard]] auto getLifetimeWouldHaveSkippedTuningTime() const {
    return wouldHaveSkippedTuningTimeLifetime;
  }

  [[nodiscard]] auto getLifetimeTuningTime() const {
    return tuningTimeLifetime;
  }

 private:
  void verifyCurrentConfigTime() const {
    for (const auto &order : _lastApplicableConfigurationOrders) {
      bool shouldBeBetter;
      if (order.smaller.matches(*_currentConfig)) {
        shouldBeBetter = true;
      } else if (order.greater.matches(*_currentConfig)) {
        shouldBeBetter = false;
      } else {
        continue;
      }

      auto currentConfigTime = _traversalTimes.at(*_currentConfig);
      const auto &comparePattern = shouldBeBetter ? order.greater : order.smaller;
      for (const auto &[otherConfig, time] : _traversalTimes) {
        bool error = false;
        if (comparePattern.matches(otherConfig) and order.haveEqualSameProperties(*_currentConfig, otherConfig)) {
          if ((shouldBeBetter and time > currentConfigTime) or
              (not shouldBeBetter and time < currentConfigTime)) {
            error = true;
          }
        }
        if (error) {
          if(shouldBeBetter) {
            _tuningErrorPrinter(order, *_currentConfig, currentConfigTime, otherConfig, time, _currentLiveInfo);
          } else {
            _tuningErrorPrinter(order, otherConfig, time, *_currentConfig, currentConfigTime, _currentLiveInfo);
          }
        }
      }
    }
  }

  std::vector<rule_syntax::ConfigurationOrder> applyRules() {
    AutoPasLog(debug, _currentLiveInfo.toString());

    std::vector<RuleVM::MemoryCell> initialStack;
    std::vector<std::pair<std::string, rule_syntax::Define>> defines{};
    for (const auto &[name, value] : _currentLiveInfo.get()) {
      initialStack.emplace_back(value);
      defines.push_back({name, {name, std::make_shared<rule_syntax::Literal>(value)}});
    }

    std::ifstream ifs{_ruleFileName};
    std::string content((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));

    rule_syntax::RuleBasedProgramParser parser{defines};
    auto [programTree, context] = parser.parse(content);

    auto generatedProgram = programTree.generateCode(context);

    RuleVM vm{};
    auto removePatterns = vm.execute(generatedProgram, initialStack);

    AutoPasLog(debug, "Remove patterns (Count {}):", removePatterns.size());
    std::vector<ConfigurationPattern> toRemovePatterns{};
    std::vector<rule_syntax::ConfigurationOrder> applicableConfigurationOrders{};
    for (const auto &patternIdx : removePatterns) {
      auto pattern = context.smallerConfigurationPatternByIndex(patternIdx);
      toRemovePatterns.push_back(pattern);
      auto str = pattern.toString();
      AutoPasLog(debug, "Remove {}", str);

      applicableConfigurationOrders.push_back(context.getConfigurationOrders().at(patternIdx));
    }

    auto newSearchSpace = _originalSearchSpace;
    _removedConfigurations.clear();
    auto& removedConfigsLocal = _removedConfigurations;
    newSearchSpace.remove_if([&toRemovePatterns, &removedConfigsLocal](const Configuration &configuration) {
      bool remove = std::any_of(toRemovePatterns.begin(), toRemovePatterns.end(),
                         [&configuration](const auto &pattern) { return pattern.matches(configuration); });
      if(remove) {
        removedConfigsLocal.insert(configuration);
      }
      return remove;
    });
    AutoPasLog(debug, "Rules remove {} out of {} configurations", _originalSearchSpace.size() - newSearchSpace.size(),
               _originalSearchSpace.size());
    if (not _verifyModeEnabled) {
      this->_searchSpace = {newSearchSpace.begin(), newSearchSpace.end()};
    }

    return applicableConfigurationOrders;
  }

 private:
  const std::list<Configuration> _originalSearchSpace;
  std::unordered_set<Configuration, ConfigHash> _removedConfigurations;
  std::vector<rule_syntax::ConfigurationOrder> _lastApplicableConfigurationOrders;
  bool _verifyModeEnabled;
  std::string _ruleFileName;

  long tuningTime = 0;
  long wouldHaveSkippedTuningTime = 0;
  long tuningTimeLifetime = 0;
  long wouldHaveSkippedTuningTimeLifetime = 0;

  PrintTuningErrorFunType _tuningErrorPrinter;

  LiveInfo _currentLiveInfo;
};
}  // namespace autopas
